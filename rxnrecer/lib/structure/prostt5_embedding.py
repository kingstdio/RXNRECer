"""
ProstT5 embedding generation for protein sequences.

Generate embeddings from FASTA files or DataFrames using ProstT5 model.
Supports both amino acid sequences and 3Di sequences.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Dict, Tuple, Optional

import pandas as pd
import torch
from tqdm import tqdm
from transformers import T5EncoderModel, T5Tokenizer


def get_device() -> torch.device:
    """Returns the appropriate device (cuda, mps, or cpu)."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    elif torch.backends.mps.is_available():
        return torch.device("mps")
    else:
        return torch.device("cpu")


def read_fasta(
    fasta_path: str, split_char: str = "|", id_field: int = 0, is_3di: bool = False
) -> Dict[str, str]:
    """
    Read sequences from a FASTA file.

    Args:
        fasta_path: Path to FASTA file
        split_char: Character to split header fields
        id_field: Index of ID field in header
        is_3di: Whether sequences are 3Di (lowercase)

    Returns:
        Dictionary mapping protein IDs to sequences
    """
    sequences = {}
    with open(fasta_path, "r") as fasta_f:
        for line in fasta_f:
            if line.startswith(">"):
                uniprot_id = line.replace(">", "").strip().split(split_char)[id_field]
                uniprot_id = uniprot_id.replace("/", "_").replace(".", "_")
                sequences[uniprot_id] = ""
            else:
                if is_3di:
                    sequences[uniprot_id] += "".join(line.split()).replace("-", "").lower()
                else:
                    sequences[uniprot_id] += "".join(line.split()).replace("-", "")
    return sequences


def load_model_and_prepare(
    half_precision: bool = False, model_dir: str = "Rostlab/ProstT5"
) -> Tuple[T5EncoderModel, T5Tokenizer]:
    """
    Load ProstT5 model and tokenizer.

    Args:
        half_precision: Whether to use float16
        model_dir: Model directory or HuggingFace model name

    Returns:
        Tuple of (model, tokenizer)
    """
    device = get_device()
    print(f"Loading T5 model and tokenizer from: {model_dir}")

    model = T5EncoderModel.from_pretrained(model_dir).to(device)
    model = model.eval()

    if half_precision:
        model = model.half()
        print("Using model in half-precision!")

    vocab = T5Tokenizer.from_pretrained(model_dir, do_lower_case=False)

    return model, vocab


def process_batches(
    sequences,
    model: T5EncoderModel,
    vocab: T5Tokenizer,
    prefix: str,
    per_protein: bool = True,
    max_residues: int = 4000,
    max_seq_len: int = 2700,
    max_batch: int = 40,
    device: Optional[torch.device] = None,
) -> Dict[str, numpy.ndarray]:
    """
    Process sequences in batches and generate embeddings.

    Args:
        sequences: List of tuples (id, sequence)
        model: ProstT5 model
        vocab: Tokenizer
        prefix: Prefix to prepend to sequences
        per_protein: Whether to average embeddings per protein
        max_residues: Maximum residues per batch
        max_seq_len: Maximum sequence length
        max_batch: Maximum sequences per batch
        device: torch device

    Returns:
        Dictionary mapping protein IDs to embeddings
    """
    import numpy

    if device is None:
        device = get_device()

    emb_dict = {}
    batch = []
    start = time.time()

    for seq_idx, (identifier, seq) in tqdm(enumerate(sequences, 1), total=len(sequences), desc="Processing Sequences"):
        seq = seq.replace("U", "X").replace("Z", "X").replace("O", "X")
        seq_len = len(seq)
        seq = f"{prefix} {' '.join(seq)}"
        batch.append((identifier, seq, seq_len))

        n_res_batch = sum(s_len for _, _, s_len in batch) + seq_len

        if len(batch) >= max_batch or n_res_batch >= max_residues or seq_idx == len(sequences) or seq_len > max_seq_len:
            pdb_ids, seqs, seq_lens = zip(*batch)
            batch = []

            token_encoding = vocab.batch_encode_plus(
                seqs, add_special_tokens=True, padding="longest", return_tensors="pt"
            ).to(device)

            try:
                with torch.no_grad():
                    embedding_repr = model(
                        token_encoding.input_ids, attention_mask=token_encoding.attention_mask
                    )
            except RuntimeError:
                print(f"RuntimeError during embedding for {identifier} (Length={seq_len})")
                continue

            for batch_idx, identifier in enumerate(pdb_ids):
                s_len = seq_lens[batch_idx]
                emb = embedding_repr.last_hidden_state[batch_idx, 1 : s_len + 1]

                if per_protein:
                    emb = emb.mean(dim=0)

                emb_dict[identifier] = emb.detach().cpu().numpy().squeeze()

                if len(emb_dict) == 1:
                    print(f"Example: Embedded protein {identifier} (Length={s_len}) to embedding of shape {emb.shape}")

    end = time.time()
    print("########################################")
    print(f"Total time: {end - start:.2f} seconds")
    print(f"Time per protein: {(end - start) / len(emb_dict):.4f} seconds")
    print("########################################")

    return emb_dict


def get_embeddings(
    seq_path: str,
    split_char: str = "|",
    id_field: int = 0,
    per_protein: bool = True,
    half_precision: bool = False,
    is_3di: bool = False,
    max_residues: int = 4000,
    max_seq_len: int = 2700,
    max_batch: int = 40,
    device: str = "cuda",
    model_dir: str = "Rostlab/ProstT5",
) -> Dict[str, numpy.ndarray]:
    """
    Generate embeddings from a FASTA file.

    Args:
        seq_path: Path to FASTA file
        split_char: Character to split header fields
        id_field: Index of ID field in header
        per_protein: Whether to average embeddings per protein
        half_precision: Whether to use float16
        is_3di: Whether input is 3Di sequences
        max_residues: Maximum residues per batch
        max_seq_len: Maximum sequence length
        max_batch: Maximum sequences per batch
        device: Device string ('cuda' or 'cpu')
        model_dir: Model directory

    Returns:
        Dictionary mapping protein IDs to embeddings
    """
    import numpy

    seq_dict = read_fasta(seq_path, split_char, id_field, is_3di)
    sequences = sorted(seq_dict.items(), key=lambda kv: len(kv[1]), reverse=True)
    prefix = "<fold2AA>" if is_3di else "<AA2fold>"

    model, vocab = load_model_and_prepare(half_precision, model_dir)

    return process_batches(
        sequences, model, vocab, prefix, per_protein, max_residues, max_seq_len, max_batch, torch.device(device)
    )


def get_embeddings_with_df(
    df_token_with_id: pd.DataFrame,
    per_protein: bool = True,
    half_precision: bool = False,
    is_3di: bool = False,
    max_residues: int = 4000,
    max_seq_len: int = 2700,
    max_batch: int = 40,
    device: str = "cuda",
    model_dir: str = "Rostlab/ProstT5",
) -> Dict[str, numpy.ndarray]:
    """
    Generate embeddings from a DataFrame.

    Args:
        df_token_with_id: DataFrame with ['uniprot_id', 'sequence'] columns
        per_protein: Whether to average embeddings per protein
        half_precision: Whether to use float16
        is_3di: Whether input is 3Di sequences
        max_residues: Maximum residues per batch
        max_seq_len: Maximum sequence length
        max_batch: Maximum sequences per batch
        device: Device string
        model_dir: Model directory

    Returns:
        Dictionary mapping protein IDs to embeddings
    """
    import numpy

    if is_3di:
        df_token_with_id = df_token_with_id.copy()
        df_token_with_id["sequence"] = df_token_with_id["sequence"].str.lower()

    sequences = df_token_with_id.sort_values(
        by="sequence", key=lambda col: col.str.len(), ascending=False
    ).to_records(index=False)
    prefix = "<fold2AA>" if is_3di else "<AA2fold>"

    model, vocab = load_model_and_prepare(half_precision, model_dir)

    return process_batches(
        sequences, model, vocab, prefix, per_protein, max_residues, max_seq_len, max_batch, torch.device(device)
    )


def save_embeddings(emb_dict: Dict, emb_path: str) -> None:
    """
    Save embeddings to a Feather file.

    Args:
        emb_dict: Dictionary mapping IDs to embedding arrays
        emb_path: Output file path
    """
    import numpy

    res = pd.DataFrame(emb_dict).T
    res.columns = [f"f{item}" for item in range(1, 1025)]
    res.insert(0, "uniprot_id", res.index)
    res.reset_index(drop=True, inplace=True)
    res.to_feather(emb_path)
    print(f"File saved to: {emb_path} successfully")

    print("\n############# STATS #############")
    print(f"Total number of embeddings: {len(emb_dict)}")


def create_arg_parser() -> argparse.ArgumentParser:
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description="embed.py creates ProstT5-Encoder embeddings for a given text "
        "file containing sequence(s) in FASTA-format."
        "Example: python embed.py --input /path/to/some_sequences.fasta --output /path/to/some_embeddings.feather --half 1 --is_3Di 0 --per_protein 1"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=False,
        type=str,
        default="/tmp/3di.fasta",
        help="A path to a fasta-formatted text file containing protein sequence(s).",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=False,
        type=str,
        default="/tmp/3di_embeddings.feather",
        help="A path for saving the created embeddings as NumPy npz file.",
    )

    parser.add_argument(
        "--split_char",
        type=str,
        default="|",
        help="The character for splitting the FASTA header in order to retrieve "
        "the protein identifier. Should be used in conjunction with --id. Default: '|'",
    )

    parser.add_argument(
        "--id",
        type=int,
        default=0,
        help="The index for the uniprot identifier field after splitting the "
        "FASTA header after each symbole in ['|', '#', ':', ' '].",
    )

    parser.add_argument(
        "--half",
        type=int,
        default=0,
        help="Whether to use half precision (1) or not (0). Default: 0",
    )

    parser.add_argument(
        "--is_3Di",
        type=int,
        default=0,
        help="Whether the input sequences are 3Di (1) or AA (0). Default: 0",
    )

    parser.add_argument(
        "--per_protein",
        type=int,
        default=1,
        help="Whether to average the per-residue embeddings for each protein (1) or not (0). Default: 1",
    )

    parser.add_argument(
        "--device",
        type=str,
        default="cuda",
        help="Device to use for computation ('cuda' or 'cpu'). Default: cuda",
    )

    return parser


def main():
    """Main entry point for command line usage."""
    parser = create_arg_parser()
    args = parser.parse_args()

    half_precision = bool(args.half)
    is_3di = bool(args.is_3Di)
    per_protein = bool(args.per_protein)

    emb_dict = get_embeddings(
        seq_path=args.input,
        split_char=args.split_char,
        id_field=args.id,
        per_protein=per_protein,
        half_precision=half_precision,
        is_3di=is_3di,
        device=args.device,
    )

    save_embeddings(emb_dict, args.output)


if __name__ == "__main__":
    main()
