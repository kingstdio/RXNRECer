"""
ESMfold batch processing for protein structure prediction.

Functions for processing FASTA files, chunking sequences by length,
and generating SLURM scripts for batch ESMfold prediction on HPC.
"""

from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd
from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def read_fasta2df(fasta_path: str) -> pd.DataFrame:
    """
    Read FASTA file and return DataFrame with id, description and sequence.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        DataFrame with columns: id, description, sequence
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    data = [{
        "id": r.id,
        "description": r.description[len(r.id):].strip(),
        "sequence": str(r.seq)
    } for r in records]
    return pd.DataFrame(data)


def move_pdb_files(src: str, target: str) -> None:
    """
    Move PDB file from source to target.

    Args:
        src: Source file path
        target: Target file path
    """
    if os.path.exists(src):
        Path(target).parent.mkdir(parents=True, exist_ok=True)
        shutil.move(src, target)
    else:
        logging.warning(f"File does not exist, skipping: {src}")


def chunk_fasta(
    df: pd.DataFrame,
    sequences_per_chunk_default: int = 400,
    partition_chunk_sizes: Optional[Dict[str, int]] = None,
    outdir: Path = Path("./chunks")
) -> Dict[str, Dict]:
    """
    Split FASTA file by sequence length partitions:
      - qgpu_a40: sequence length <= 900 bp
      - qgpu_a800: sequence length in (900, 2600] bp
    Sequences exceeding 2600 bp are not processed.

    Args:
        df: DataFrame with 'id', 'sequence' columns (should have 'seqlen' column added beforehand)
        sequences_per_chunk_default: Default number of sequences per chunk (for a40 partition)
        partition_chunk_sizes: dict, specifying sequences per chunk for each partition
                              e.g., {"qgpu_a40": 400, "qgpu_a800": 100}
        outdir: Output directory

    Returns:
        Dictionary with partition info including chunks count, total sequences, files, options
    """
    outdir.mkdir(parents=True, exist_ok=True)

    partitions = {
        "qgpu_a40": df[df["seqlen"] <= 900],
        "qgpu_a800": df[(df["seqlen"] > 900) & (df["seqlen"] <= 2600)]
    }
    summary = {}

    if partition_chunk_sizes is None:
        partition_chunk_sizes = {
            "qgpu_a40": sequences_per_chunk_default,
            "qgpu_a800": 100
        }

    for partition, subset in partitions.items():
        subset = subset.reset_index(drop=True)
        seqs_per_chunk = partition_chunk_sizes.get(partition, sequences_per_chunk_default)
        chunk_files = []
        chunk_options = {}
        num_chunks = 0

        for i in range(0, len(subset), seqs_per_chunk):
            chunk_df = subset.iloc[i:i+seqs_per_chunk].copy()
            oversized = chunk_df[chunk_df['seqlen'] > 5000]
            chunk_df = chunk_df[chunk_df['seqlen'] <= 5000]
            if chunk_df.empty:
                continue

            num_chunks += 1
            chunk_path = outdir / f"{partition}_chunk_{num_chunks}.fasta"
            with open(chunk_path, "w") as f:
                lines = [f">{row['id']}\n{row['sequence']}\n" for _, row in chunk_df.iterrows()]
                f.writelines(lines)
            chunk_files.append(str(chunk_path))

            min_len = chunk_df['seqlen'].min()
            max_len = chunk_df['seqlen'].max()
            avg_len = chunk_df['seqlen'].mean()
            extra_arg = ""
            if partition == "qgpu_a800":
                if max_len <= 2000:
                    extra_arg = "--chunk-size=128"
                elif max_len <= 2600:
                    extra_arg = "--chunk-size=64"
                else:
                    extra_arg = "⚠️ Too long, consider splitting manually"
            chunk_options[num_chunks] = extra_arg
            logging.info(f"Chunk {num_chunks}: {len(chunk_df)} seqs, min = {min_len}, max = {max_len}, "
                         f"avg = {avg_len:.1f}, partition = {partition}" +
                         (f", extra = {extra_arg}" if extra_arg else ""))

            if not oversized.empty:
                logging.warning("The following sequences exceed length 5000 and were not written to chunk files:")
                for _, row in oversized.iterrows():
                    logging.warning(f"  ID: {row['id']}, Length: {row['seqlen']}")

        summary[partition] = {
            "chunks": num_chunks,
            "total_seqs": len(subset),
            "files": chunk_files,
            "options": chunk_options
        }
    return summary


def generate_all_slurm_scripts(
    chunk_info: Dict,
    input_dir: str,
    output_dir: str,
    log_dir: str,
    script_dir: str,
    container_path: str = "/hpcfs/fpublic/container/singularity/app/esmfold/esmfold.sif",
    esmfold_script: str = "/esmfold.sh"
) -> None:
    """
    Generate SLURM submission scripts for all partitions.

    Args:
        chunk_info: Dictionary with partition info from chunk_fasta
        input_dir: Input directory containing chunk FASTA files
        output_dir: Output directory for PDB files
        log_dir: Directory for log files
        script_dir: Directory to save SLURM scripts
        container_path: Path to Singularity container
        esmfold_script: Path to ESMfold shell script inside container
    """
    Path(script_dir).mkdir(parents=True, exist_ok=True)
    Path(log_dir).mkdir(parents=True, exist_ok=True)

    for partition, info in chunk_info.items():
        files = info["files"]
        options = info.get("options", {})
        if not files:
            logging.warning(f"No chunks for {partition}, skipping.")
            continue

        script_path = Path(script_dir) / f"submit_{partition}.slurm"
        with open(script_path, "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH --job-name=esmfold-{partition}
#SBATCH --partition={partition}
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --gres=gpu:1
#SBATCH --output={log_dir}/esmfold-{partition}-%A_%a.out
#SBATCH --error={log_dir}/esmfold-{partition}-%A_%a.err
#SBATCH --array=1-{len(files)}

INPUT_DIR="{input_dir}"
OUTPUT_DIR="{output_dir}"
INPUT_FILE="{partition}_chunk_${{SLURM_ARRAY_TASK_ID}}.fasta"
EXTRA_ARGS=""
""")
            f.write("case $SLURM_ARRAY_TASK_ID in\n")
            for i, arg in options.items():
                if arg and not arg.startswith("⚠️"):
                    f.write(f"    {i}) EXTRA_ARGS='{arg}' ;;\n")
            f.write("esac\n\n")
            f.write(f"""mkdir -p $OUTPUT_DIR

echo "[`date`] Running on node: `hostname`"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running command:"
echo "singularity exec --nv --bind /hpcfs:/hpcfs {container_path} \\
    bash {esmfold_script} -i $INPUT_DIR/$INPUT_FILE -o $OUTPUT_DIR $EXTRA_ARGS"

singularity exec --nv --bind /hpcfs:/hpcfs {container_path} \\
    bash {esmfold_script} -i $INPUT_DIR/$INPUT_FILE -o $OUTPUT_DIR $EXTRA_ARGS
""")
        logging.info(f"[{partition}] SLURM script generated: {script_path}")
        logging.info(f"Submit command: sbatch {script_path}")


def load_existing_pdb_files(pdb_dir: str) -> list:
    """
    Return list of stems for all PDB files in target directory.

    Args:
        pdb_dir: Directory containing PDB files

    Returns:
        List of PDB file stems (without extension)
    """
    return [f.stem for f in Path(pdb_dir).rglob("*.pdb")]


def prepare_fasta_for_esmfold(
    fasta_path: str,
    target_pdb_dir: str,
    max_seq_len: int = 2600
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read FASTA file and filter out sequences that already have PDB files.

    Args:
        fasta_path: Path to input FASTA file
        target_pdb_dir: Directory with existing PDB files
        max_seq_len: Maximum sequence length to process

    Returns:
        Tuple of (filtered_df, too_long_df)
    """
    records_df = read_fasta2df(fasta_path)
    existing_pdb_files = load_existing_pdb_files(target_pdb_dir)
    records_df = records_df[~records_df.id.isin(existing_pdb_files)].reset_index(drop=True)
    records_df['seqlen'] = records_df.sequence.str.len()
    records_df['sequence'] = records_df.sequence.str.upper().str.strip()
    records_df = records_df.sort_values("seqlen").reset_index(drop=True)

    too_long = records_df[records_df["seqlen"] > max_seq_len]
    if not too_long.empty:
        logging.warning("The following sequences exceed length 2600 and were not processed:")
        for _, row in too_long.iterrows():
            logging.warning(f"  ID: {row['id']}, Length: {row['seqlen']}")
    records_df = records_df[records_df["seqlen"] <= max_seq_len].reset_index(drop=True)

    return records_df, too_long


def step_by_step_run(
    pdb_output_dir: str,
    fasta_path: str,
    tempdir: str,
    slurm_dir: str,
    log_dir: str,
    target_pdb_dir: str,
    sequences_per_chunk_default: int = 400,
    partition_chunk_sizes: Optional[Dict[str, int]] = None
) -> Dict[str, Dict]:
    """
    Run the complete ESMfold batch processing workflow.

    Args:
        pdb_output_dir: Output directory for PDB files
        fasta_path: Path to input FASTA file
        tempdir: Temporary directory for chunk files
        slurm_dir: Directory for SLURM scripts
        log_dir: Directory for log files
        target_pdb_dir: Target directory for final PDB files
        sequences_per_chunk_default: Default chunk size
        partition_chunk_sizes: Partition-specific chunk sizes

    Returns:
        Dictionary with chunk information
    """
    for path in [tempdir, slurm_dir, log_dir, pdb_output_dir, target_pdb_dir]:
        logging.info(f"Created directory: {path}")
        Path(path).mkdir(parents=True, exist_ok=True)

    records_df, too_long = prepare_fasta_for_esmfold(fasta_path, target_pdb_dir)

    chunk_info = chunk_fasta(
        records_df,
        sequences_per_chunk_default=sequences_per_chunk_default,
        partition_chunk_sizes=partition_chunk_sizes,
        outdir=Path(tempdir)
    )
    generate_all_slurm_scripts(
        chunk_info,
        input_dir=tempdir,
        output_dir=pdb_output_dir,
        log_dir=log_dir,
        script_dir=slurm_dir
    )

    return chunk_info
