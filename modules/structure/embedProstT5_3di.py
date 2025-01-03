#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 14:27:44 2023

@author: mheinzinger
"""

import argparse
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


def load_model_and_prepare(half_precision, model_dir="Rostlab/ProstT5"):
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
    print("Loading T5 model and tokenizer from: {}".format(model_dir))
    
    device = get_device()
    # Load the model
    model = T5EncoderModel.from_pretrained(model_dir).to(device)
    model = model.eval()  # Set to evaluation mode

    # Use half-precision if specified
    if half_precision:
        model = model.half()
        print("Using model in half-precision!")

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
                if len(emb_dict) == 1:
                    print(f"Example: Embedded protein {identifier} (Length={s_len}) to embedding of shape {emb.shape}")

    end = time.time()
    print('########################################')
    print(f"Total time: {end - start:.2f} seconds")
    print(f"Time per protein: {(end - start) / len(emb_dict):.4f} seconds")
    print('########################################')

    return emb_dict



def get_embeddings(seq_path, split_char, id_field, per_protein, half_precision, is_3Di,
                   max_residues=4000, max_seq_len=2700, max_batch=40, device='cuda'):
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
    model, vocab = load_model_and_prepare(half_precision, model_dir='/hpcfs/fhome/shizhenkun/.cache/huggingface/hub/models--Rostlab--ProstT5/snapshots/d7d097d5bf9a993ab8f68488b4681d6ca70db9e5')

    # Process batches
    return process_batches(sequences, model, vocab, prefix, per_protein, max_residues, max_seq_len, max_batch, torch.device(device))


def get_embeddings_with_df(df_token_with_id, per_protein, half_precision, is_3Di,
                           max_residues=4000, max_seq_len=2700, max_batch=40, device='cuda'):
    """
    Generate embeddings from a DataFrame.

    Parameters:
        df_token_with_id (pd.DataFrame): A DataFrame with columns ['uniprot_id', 'sequence'].
        model_dir (str): Path to the model directory.
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
    # Preprocess sequences
    if is_3Di:
        df_token_with_id['sequence'] = df_token_with_id['sequence'].str.lower()
    sequences = df_token_with_id.sort_values(by='sequence', key=lambda col: col.str.len(), ascending=False).to_records(index=False)
    prefix = "<fold2AA>" if is_3Di else "<AA2fold>"

    # Load model and vocab
    model, vocab = load_model_and_prepare(half_precision)

    # Process batches
    return process_batches(sequences, model, vocab, prefix, per_protein, max_residues, max_seq_len, max_batch, torch.device(device))




def save(emb_dict, emb_path):
    res = pd.DataFrame(emb_dict).T
    res.columns=[f'f{item}' for item in range(1, 1025)]
    res.insert(0, 'uniprot_id', res.index)
    res.reset_index(drop=True, inplace=True)
    res.to_feather(emb_path)
    print(f'File saved to:{emb_path} successfully')
    # with h5py.File(str(emb_path), "w") as hf:
    #     for sequence_id, embedding in emb_dict.items():
    #         # noinspection PyUnboundLocalVariable
    #         hf.create_dataset(sequence_id, data=embedding)

    print('\n############# STATS #############')
    print('Total number of embeddings: {}'.format(len(emb_dict)))



def create_arg_parser():
    """"Creates and returns the ArgumentParser object."""

    # Instantiate the parser
    parser = argparse.ArgumentParser(description=( 
            'embed.py creates ProstT5-Encoder embeddings for a given text '+
            ' file containing sequence(s) in FASTA-format.' +
            'Example: python embed.py --input /path/to/some_sequences.fasta --output /path/to/some_embeddings.feather --half 1 --is_3Di 0 --per_protein 1' ) )
    
    # Required positional argument
    parser.add_argument( '-i', '--input', required=False, type=str,default='/tmp/3di.fasta',
                    help='A path to a fasta-formatted text file containing protein sequence(s).')

    # Optional positional argument
    parser.add_argument( '-o', '--output', required=False, type=str, default='/tmp/3di_embeddings.feather',
                    help='A path for saving the created embeddings as NumPy npz file.')


    # Optional argument
    parser.add_argument('--split_char', type=str, 
                    default='!',
                    help='The character for splitting the FASTA header in order to retrieve ' +
                        "the protein identifier. Should be used in conjunction with --id." +
                        "Default: '!' ")
    
    # Optional argument
    parser.add_argument('--id', type=int, 
                    default=0,
                    help='The index for the uniprot identifier field after splitting the ' +
                        "FASTA header after each symbole in ['|', '#', ':', ' ']." +
                        'Default: 0')
    # Optional argument
    parser.add_argument('--per_protein', type=int, 
                    default=1,
                    help="Whether to return per-residue embeddings (0: default) or the mean-pooled per-protein representation (1).")
        
    parser.add_argument('--half', type=int, 
                    default=1,
                    help="Whether to use half_precision or not. Default: 1 (full-precision)")
    
    parser.add_argument('--is_3Di', type=int, 
                    default=1,
                    help="Whether to create embeddings for 3Di or AA file. Default: 0 (generate AA-embeddings)")
    
    return parser

def main():
    parser     = create_arg_parser()
    args       = parser.parse_args()
    
    seq_path   = Path( args.input ) # path to input FASTAS
    emb_path   = Path( args.output) # path where embeddings should be stored
    model_dir  = "Rostlab/ProstT5"
    
    split_char = args.split_char
    id_field   = args.id

    per_protein    = False if int(args.per_protein) == 0 else True
    half_precision = False if int(args.half)        == 0 else True
    is_3Di         = False if int(args.is_3Di)      == 0 else True

    rpr = get_embeddings(seq_path,  split_char,  id_field,  per_protein=per_protein, half_precision=half_precision,  is_3Di=is_3Di)
    save(emb_dict=rpr, emb_path=emb_path)


if __name__ == '__main__':
    main()
