import time
from pathlib import Path
import torch
from tqdm import tqdm
import pandas as pd
from transformers import T5EncoderModel, T5Tokenizer


def get_device():
    """Returns the appropriate device (cuda, mps, or cpu)."""
    if torch.cuda.is_available():
        return torch.device('cuda')
    elif torch.backends.mps.is_available():
        return torch.device('mps')
    else:
        return torch.device('cpu')


# Original function to read sequences from FASTA remains the same
def read_fasta(fasta_path, split_char, id_field, is_3Di):
    """
    Reads a FASTA file containing multiple sequences.
    Returns a dictionary of sequences.
    """
    sequences = {}
    with open(fasta_path, 'r') as fasta_f:
        for line in fasta_f:
            if line.startswith('>'):
                uniprot_id = line.replace('>', '').strip().split(split_char)[id_field]
                uniprot_id = uniprot_id.replace("/", "_").replace(".", "_")
                sequences[uniprot_id] = ''
            else:
                if is_3Di:
                    sequences[uniprot_id] += ''.join(line.split()).replace("-", "").lower()
                else:
                    sequences[uniprot_id] += ''.join(line.split()).replace("-", "")
    return sequences


def load_model_and_prepare(model_dir="Rostlab/ProstT5"):
    """
    Load the T5 model and tokenizer, optionally set to half-precision.

    Parameters:
        model_dir (str): Path to the model directory.
        half_precision (bool): Whether to use half-precision (float16).
        device (torch.device): The device to load the model onto (e.g., 'cpu' or 'cuda').

    Returns:
        model: The T5 model.
        vocab: The tokenizer.
    """
    
    # print("Loading T5 model and tokenizer from: {}".format(model_dir))
    
    device = get_device()
    if str(device) =='cuda':
        half_precision = True
    else:
        half_precision = False
    
    # Load the model
    model = T5EncoderModel.from_pretrained(model_dir).to(device)
    model = model.eval()  # Set to evaluation mode

    # Use half-precision if specified
    if half_precision:
        model = model.half()
        # print("Using model in half-precision!")

    # Load the tokenizer
    vocab = T5Tokenizer.from_pretrained(model_dir, do_lower_case=False)

    return model, vocab


def process_batches(sequences, model, vocab, prefix, per_protein, max_residues, max_seq_len, max_batch, device):
    """
    Process sequences in batches and generate embeddings.

    Parameters:
        sequences (list): List of tuples (id, sequence).
        model: The T5 model.
        vocab: The tokenizer.
        prefix (str): Prefix to prepend to sequences.
        per_protein (bool): Whether to average embeddings per protein.
        max_residues (int): Maximum residues in a batch.
        max_seq_len (int): Maximum length of a sequence.
        max_batch (int): Maximum number of sequences per batch.
        device (torch.device): The device to run inference on.

    Returns:
        dict: A dictionary mapping protein IDs to embeddings.
    """
    emb_dict = {}
    batch = []
    start = time.time()

    # Wrap the sequence processing loop with tqdm to show progress
    for seq_idx, (identifier, seq) in tqdm(enumerate(sequences, 1), total=len(sequences), desc="Processing Sequences"):
        # Replace non-standard amino acids with 'X'
        seq = seq.replace('U', 'X').replace('Z', 'X').replace('O', 'X')
        seq_len = len(seq)
        seq = f"{prefix} {' '.join(seq)}"
        batch.append((identifier, seq, seq_len))

        # Calculate total residues in the current batch
        n_res_batch = sum(s_len for _, _, s_len in batch) + seq_len

        # Process the batch if limits are reached
        if len(batch) >= max_batch or n_res_batch >= max_residues or seq_idx == len(sequences) or seq_len > max_seq_len:
            pdb_ids, seqs, seq_lens = zip(*batch)
            batch = []

            # Tokenize sequences
            token_encoding = vocab.batch_encode_plus(
                seqs, add_special_tokens=True, padding="longest", return_tensors='pt'
            ).to(device)

            # Generate embeddings
            try:
                with torch.no_grad():
                    embedding_repr = model(
                        token_encoding.input_ids,
                        attention_mask=token_encoding.attention_mask
                    )
            except RuntimeError:
                print(f"RuntimeError during embedding for {identifier} (Length={seq_len})")
                continue

            # Extract embeddings
            for batch_idx, identifier in enumerate(pdb_ids):
                s_len = seq_lens[batch_idx]
                emb = embedding_repr.last_hidden_state[batch_idx, 1:s_len + 1]

                # Optionally average embeddings
                if per_protein:
                    emb = emb.mean(dim=0)

                # Save embedding
                emb_dict[identifier] = emb.detach().cpu().numpy().squeeze()

                # Log example embedding
                # if len(emb_dict) == 1:
                #     print(f"Example: Embedded protein {identifier} (Length={s_len}) to embedding of shape {emb.shape}")

    end = time.time()
    # print('########################################')
    # print(f"Total time: {end - start:.2f} seconds")
    # print(f"Time per protein: {(end - start) / len(emb_dict):.4f} seconds")
    # print('########################################')

    return emb_dict


def get_embeddings(seq_path, id_field, per_protein, is_3Di=False, split_char='!', max_residues=4000, max_seq_len=2700, max_batch=40, device='cuda'):
    """
    Generate embeddings from a FASTA file.

    Parameters:
        seq_path (str): Path to the FASTA file.
        model_dir (str): Path to the model directory.
        split_char (str): Character used to split sequences in FASTA headers.
        id_field (int): Field in the header to use as the identifier.
        per_protein (bool): Whether to average embeddings per protein.
        half_precision (bool): Whether to use half-precision.
        is_3Di (bool): Whether the input is 3Di sequences.
        max_residues (int): Maximum residues in a batch.
        max_seq_len (int): Maximum sequence length.
        max_batch (int): Maximum sequences per batch.
        device (str): Device to use for computation ('cuda' or 'cpu').

    Returns:
        dict: A dictionary mapping protein IDs to embeddings.
    """
    # Read sequences from FASTA
    seq_dict = read_fasta(seq_path, split_char, id_field, is_3Di)
    sequences = sorted(seq_dict.items(), key=lambda kv: len(kv[1]), reverse=True)
    prefix = "<fold2AA>" if is_3Di else "<AA2fold>"

    # Load model and vocab
    model, vocab = load_model_and_prepare(model_dir='/hpcfs/fhome/shizhenkun/.cache/huggingface/hub/models--Rostlab--ProstT5/snapshots/d7d097d5bf9a993ab8f68488b4681d6ca70db9e5')

    # Process batches
    return process_batches(sequences, model, vocab, prefix, per_protein, max_residues, max_seq_len, max_batch, torch.device(device))


def save(emb_dict, emb_path):
    res = pd.DataFrame(emb_dict).T
    res.columns=[f'f{item}' for item in range(1, 1025)]
    res.insert(0, 'uniprot_id', res.index)
    res.reset_index(drop=True, inplace=True)
    res.to_feather(emb_path)
    print(f'File saved to:{emb_path} successfully')


def get_embd_seq(seqdfwithid, batch_szise=40):
    seq_dict = dict(zip(seqdfwithid['id'], seqdfwithid['seq'].str.upper()))
    sequences = sorted(seq_dict.items(), key=lambda kv: len(kv[1]), reverse=True)
    prefix = '<AA2fold>'
    model, vocab = load_model_and_prepare(model_dir='/hpcfs/fhome/shizhenkun/.cache/huggingface/hub/models--Rostlab--ProstT5/snapshots/d7d097d5bf9a993ab8f68488b4681d6ca70db9e5')
    rpr = process_batches(sequences, model, vocab, prefix, per_protein=True, max_residues=4000, max_seq_len=2700, max_batch=batch_szise, device =get_device())
    rpr = pd.DataFrame(list(rpr.items()), columns=['id', 't5']) 
    return rpr
    

def main():
    # Define variables directly instead of using argparse
    seq_path = Path('/hpcfs/fhome/shizhenkun/sample100.fasta')  # path to input FASTAS
    emb_path = Path('/tmp/3di_embeddings.feather')  # path where embeddings should be stored
    
    id_field = 0  # field index for the uniprot identifier
    per_protein = True  # Whether to average embeddings per protein
    is_3Di = False  # Whether the input is 3Di sequences

    # Generate embeddings
    rpr = get_embeddings(seq_path=seq_path, id_field=id_field, per_protein=per_protein, is_3Di=is_3Di)
    save(emb_dict=rpr, emb_path=emb_path)


if __name__ == '__main__':
    main()
