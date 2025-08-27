import sys,os
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, f'{project_root}/../')
import time
from pathlib import Path
import torch
from tqdm import tqdm
import pandas as pd
import warnings
from rxnrecer.config import config as cfg
# Suppress noisy warnings and logging
warnings.filterwarnings("ignore", category=FutureWarning)
from transformers import T5EncoderModel, T5Tokenizer, T5Config
from transformers.utils import logging
logging.set_verbosity_error()


def get_device():
    """Return the best available device."""
    if torch.cuda.is_available():
        return torch.device('cuda')
    elif torch.backends.mps.is_available():
        return torch.device('mps')
    return torch.device('cpu')


def read_fasta(fasta_path, split_char='!', id_field=0, is_3Di=False):
    """Read FASTA file and return a dict of {id: sequence}."""
    sequences = {}
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                uniprot_id = line[1:].strip().split(split_char)[id_field].replace("/", "_").replace(".", "_")
                sequences[uniprot_id] = ''
            else:
                seq_line = ''.join(line.split()).replace("-", "")
                sequences[uniprot_id] += seq_line.lower() if is_3Di else seq_line
    return sequences


def load_model_and_prepare(model_dir):
    """Load and configure T5 model + tokenizer."""
    device = get_device()
    config = T5Config.from_pretrained(model_dir)
    config.is_encoder_decoder = False
    config.architectures = ['T5EncoderModel']
    config.use_cache = False

    model = T5EncoderModel.from_pretrained(model_dir, config=config).to(device).eval()
    if device.type == 'cuda':
        model = model.half()

    tokenizer = T5Tokenizer.from_pretrained(model_dir, do_lower_case=False, legacy=True)
    return model, tokenizer


def process_batches(sequences, model, tokenizer, prefix, per_protein, max_residues, max_seq_len, max_batch, device):
    emb_dict = {}
    batch = []
    start = time.time()
    for seq_idx, (identifier, seq) in tqdm(enumerate(sequences, 1), total=len(sequences), desc="Processing Sequences"):
        seq = seq.replace('U', 'X').replace('Z', 'X').replace('O', 'X')
        seq_len = len(seq)
        seq = f"{prefix} {' '.join(seq)}"
        batch.append((identifier, seq, seq_len))
        n_res_batch = sum(s_len for _, _, s_len in batch)
        if len(batch) >= max_batch or n_res_batch >= max_residues or seq_idx == len(sequences) or seq_len > max_seq_len:
            pdb_ids, seqs, seq_lens = zip(*batch)
            batch = []
            tokenized = tokenizer.batch_encode_plus(seqs, add_special_tokens=True, padding="longest", return_tensors='pt').to(device)
            try:
                with torch.no_grad():
                    outputs = model(tokenized.input_ids, attention_mask=tokenized.attention_mask)
            except RuntimeError:
                print(f"RuntimeError on {identifier} (Length={seq_len})")
                continue

            for i, pid in enumerate(pdb_ids):
                emb = outputs.last_hidden_state[i, 1:seq_lens[i]+1]
                if per_protein:
                    emb = emb.mean(dim=0)
                emb_dict[pid] = emb.detach().cpu().numpy().squeeze()

    print(f"Processed {len(emb_dict)} proteins in {time.time() - start:.2f} seconds")
    return emb_dict

def get_embeddings(seq_path, id_field=0, per_protein=True, is_3Di=False, split_char='!',
                   max_residues=4000, max_seq_len=2701, max_batch=40, device='cuda'
                   ,model_dir=cfg.CKPT_PROSTT5):
    seq_dict = read_fasta(seq_path, split_char, id_field, is_3Di)
    sequences = sorted(seq_dict.items(), key=lambda kv: len(kv[1]), reverse=True)
    prefix = "<fold2AA>" if is_3Di else "<AA2fold>"
    model, tokenizer = load_model_and_prepare(model_dir=model_dir)
    return process_batches(sequences, model, tokenizer, prefix, per_protein, max_residues, max_seq_len, max_batch, torch.device(device))

def save(emb_dict, emb_path):
    df = pd.DataFrame(emb_dict).T
    df.columns = [f'f{i+1}' for i in range(df.shape[1])]
    df.insert(0, 'uniprot_id', df.index)
    df.reset_index(drop=True, inplace=True)
    df.to_feather(emb_path)
    print(f"Embeddings saved to: {emb_path}")

def get_embd_seq(seqdfwithid, batch_size=40, max_seq_len=2500, model_dir=cfg.CKPT_PROSTT5):
    seq_dict = dict(zip(seqdfwithid['id'], seqdfwithid['seq'].str.upper().str[:max_seq_len]))
    sequences = sorted(seq_dict.items(), key=lambda kv: len(kv[1]), reverse=True)
    prefix = "<AA2fold>"
    model, tokenizer = load_model_and_prepare(model_dir=model_dir)
    result = process_batches(sequences, model, tokenizer, prefix, per_protein=True, max_residues=4000, max_seq_len=2700, max_batch=batch_size, device=get_device())
    return pd.DataFrame(result.items(), columns=['id', 't5'])

def test():
    seq_path = Path(f'{cfg.SAMPLE_DIR}/sample10.fasta')
    emb_path = Path(f'{cfg.TEMP_DIR}/3di_embeddings_sample10.feather')
    result = get_embeddings(seq_path=seq_path, id_field=0, per_protein=True, is_3Di=False)
    save(emb_dict=result, emb_path=emb_path)

if __name__ == '__main__':
    test()