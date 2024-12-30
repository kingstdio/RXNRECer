import sys,os
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1,'../../')
from config import conf as cfg

from transformers import T5Tokenizer, T5EncoderModel
import torch
import re
import pandas as pd
from tqdm import tqdm


def load_model(device):
    """
    Load the T5 model and tokenizer with exception handling.
    
    Parameters:
        device (torch.device): The device to load the model onto.
        
    Returns:
        model (T5EncoderModel): The loaded T5 encoder model.
        tokenizer (T5Tokenizer): The loaded T5 tokenizer.
    """

    tokenizer = T5Tokenizer.from_pretrained('Rostlab/ProstT5', do_lower_case=False)


    try:
        # Load the model
        model = T5EncoderModel.from_pretrained("Rostlab/ProstT5").to(device)
        if device.type == 'cpu':
            model.float()
        else:
            model.half()
    except Exception as e:
        print(f"Error loading model: {e}")
        print("Possible reasons: Network issues, incorrect model name, or insufficient GPU memory.")
        raise SystemExit("Failed to load the model. Check your internet connection, model name, or device configuration.")
    
    return model, tokenizer


def get_embd_using_3di(sequence_3di):
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    
    model,tokenizer = load_model(device=device)
    
    # prepare your protein sequences/structures as a list.
    # Amino acid sequences are expected to be upper-case ("PRTEINO" below)
    # while 3Di-sequences need to be lower-case ("strctr" below).
    sequence_3di = list(map(str.lower, sequence_3di))
    # replace all rare/ambiguous amino acids by X (3Di sequences do not have those) and introduce white-space between all sequences (AAs and 3Di)
    sequence_3di = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequence_3di]
    
    # The direction of the translation is indicated by two special tokens:
    # if you want to embed AAs, you need to prepend "<AA2fold>"
    # if you want to embed 3Di, you need to prepend "<fold2AA>"
    sequence_3di = [ "<AA2fold>" + " " + s if s.isupper() else "<fold2AA>" + " " + s  for s in sequence_3di] # this expects 3Di sequences to be already lower-case
    
    # tokenize sequences and pad up to the longest sequence in the batch
    ids = tokenizer.batch_encode_plus(sequence_3di, add_special_tokens=True, padding="longest", return_tensors='pt').to(device)
        # generate embeddings
    with torch.no_grad():
        embedding_rpr = model( ids.input_ids,  attention_mask=ids.attention_mask)

    
    embd_protein = embedding_rpr.last_hidden_state.mean(dim=1)    
    
    
    return embd_protein
    
    

def get_embd_using_3di_batch(sequence_3di, batch_size):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    sequence_3di = list(map(str.lower, sequence_3di))
    # Replace all rare/ambiguous amino acids by X (3Di sequences do not have those)
    # and introduce white-space between all sequences (AAs and 3Di)
    sequence_3di = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequence_3di]
    
    # Prepend direction tokens
    sequence_3di = ["<AA2fold> " + s if s.isupper() else "<fold2AA> " + s for s in sequence_3di]
    model, tokenizer = load_model(device=device)  # Initially load model on CPU
    
    # Tokenize sequences and pad up to the longest sequence
    ids = tokenizer.batch_encode_plus(sequence_3di, add_special_tokens=True, padding="longest", return_tensors='pt')
    # print(ids)
    input_ids = ids.input_ids
    attention_mask = ids.attention_mask
    
    embd_list = []
    
    # Process each batch to avoid OOM errors
    num_sequences = input_ids.size(0)
    model.eval()
    with torch.no_grad():
        for start_idx in tqdm(range(0, num_sequences, batch_size), desc="Processing Batches", unit="batch"):
            end_idx = start_idx + batch_size
            input_ids_batch = input_ids[start_idx:end_idx].to(device)
            attention_mask_batch = attention_mask[start_idx:end_idx].to(device)
            
            # Generate embeddings for the batch
            embedding_rpr = model(input_ids_batch, attention_mask=attention_mask_batch)


            # Take the mean of the last hidden state for each sequence in the batch
            # emb = embedding_rpr.last_hidden_state[batch_idx,1:s_len+1]
            
            embd_batch = embedding_rpr.last_hidden_state.mean(dim=1)
            embd_list.append(embd_batch.cpu())  # Move embeddings back to CPU to save GPU memory
    
    # Concatenate all the embeddings from each batch
    embd_protein = torch.cat(embd_list, dim=0)
    embd_protein_np = embd_protein.numpy()
    
    return embd_protein_np




if __name__ == '__main__':
    data = pd.read_feather(cfg.FILE_DS_3DI_LIST)
    res = get_embd_using_3di_batch(sequence_3di=data.head(10).token_3di.to_list(), batch_size=10)
    print(res)
    # res = pd.DataFrame(res)
    # target_file = '/hpcfs/fhome/shizhenkun/codebase/RXNRECer/data/featurebank/t5/embd/protein3di_embd.feater'
    # # res.to_feather(f'{cfg.TEMP_DIR}/protein3di_embd.feater')
    # res.to_feather(target_file)
    


        
    
    