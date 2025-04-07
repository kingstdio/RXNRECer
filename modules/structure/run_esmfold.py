import os
from pathlib import Path
import shutil
import pandas as pd
from Bio import SeqIO
import logging

# 配置日志记录
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def read_fasta2df(fasta_path):
    """读取FASTA文件，并返回包含id、描述和序列的DataFrame"""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    data = [{
        "id": r.id,
        "description": r.description[len(r.id):].strip(),  # 去掉前面的id
        "sequence": str(r.seq)
    } for r in records]
    return pd.DataFrame(data)

def move_pdb_files(src, target):
    """移动文件，如果源文件不存在则输出提示信息"""
    if os.path.exists(src):
        Path(target).parent.mkdir(parents=True, exist_ok=True)
        shutil.move(src, target)
    else:
        logging.warning(f"文件不存在，跳过：{src}")

def chunk_fasta(df, sequences_per_chunk_default=400, partition_chunk_sizes=None, outdir=Path("./chunks")):
    """
    拆分FASTA文件，并按照序列长度分区：
      - qgpu_a40: 序列长度 ≤ 900 bp
      - qgpu_a800: 序列长度在 (900, 2600] bp
    超过 2600 bp 的序列不处理，会在日志中提示。

    参数:
      sequences_per_chunk_default: 默认每个块的序列数（用于 a40 分区）
      partition_chunk_sizes: dict，用于指定各分区的每块序列数，例如
                             {"qgpu_a40": 400, "qgpu_a800": 100}
      outdir: 输出目录
    """
    outdir.mkdir(parents=True, exist_ok=True)
    
    # 根据长度划分分区
    partitions = {
        "qgpu_a40": df[df["seqlen"] <= 900],
        "qgpu_a800": df[(df["seqlen"] > 900) & (df["seqlen"] <= 2600)]
    }
    summary = {}
    
    if partition_chunk_sizes is None:
        partition_chunk_sizes = {
            "qgpu_a40": sequences_per_chunk_default,
            "qgpu_a800": 100  # 针对长序列使用较小的块大小
        }
    
    for partition, subset in partitions.items():
        subset = subset.reset_index(drop=True)
        seqs_per_chunk = partition_chunk_sizes.get(partition, sequences_per_chunk_default)
        chunk_files = []
        chunk_options = {}
        num_chunks = 0

        for i in range(0, len(subset), seqs_per_chunk):
            chunk_df = subset.iloc[i:i+seqs_per_chunk].copy()
            # 如果chunk内有极长序列（超过5000 bp），则不写入并提示（一般此处不会发生，因为上面已过滤）
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
                logging.warning("以下序列长度超过5000，未写入chunk文件，请手动处理：")
                for _, row in oversized.iterrows():
                    logging.warning(f"  ID: {row['id']}, Length: {row['seqlen']}")
                    
        summary[partition] = {
            "chunks": num_chunks,
            "total_seqs": len(subset),
            "files": chunk_files,
            "options": chunk_options
        }
    return summary

def generate_all_slurm_scripts(chunk_info, input_dir, output_dir, log_dir, script_dir):
    """生成所有分区的SLURM提交脚本"""
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
            f.write("""mkdir -p $OUTPUT_DIR

echo "[`date`] Running on node: `hostname`"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running command:"
echo "singularity exec --nv --bind /hpcfs:/hpcfs /hpcfs/fpublic/container/singularity/app/esmfold/esmfold.sif \\
    bash /esmfold.sh -i $INPUT_DIR/$INPUT_FILE -o $OUTPUT_DIR $EXTRA_ARGS"

singularity exec --nv --bind /hpcfs:/hpcfs /hpcfs/fpublic/container/singularity/app/esmfold/esmfold.sif \\
    bash /esmfold.sh -i $INPUT_DIR/$INPUT_FILE -o $OUTPUT_DIR $EXTRA_ARGS
""")
        logging.info(f"[{partition}] SLURM 脚本生成成功: {script_path}")
        logging.info(f"提交命令：sbatch {script_path}")

def load_existing_pdb_files(pdb_dir):
    """返回目标目录下所有PDB文件的stem列表"""
    return [f.stem for f in Path(pdb_dir).rglob("*.pdb")]

def step_by_step_run():
    # 定义各文件夹路径
    pdb_output_dir = Path("/hpcfs/fhome/shizhenkun/codebase/RXNRECer/data/structure/pdb/ncbi/")
    base_dir = Path("/hpcfs/fhome/shizhenkun/codebase/RXNRECer/case/fusarium_venenatum/data/")
    fasta_path = base_dir / "ncbi_protein.fasta"
    tempdir = base_dir / "slurm_temp"
    slurm_dir = base_dir / "slurm_script"
    log_dir = base_dir / "slurm_log"
    target_pdb_dir = Path("/hpcfs/fhome/shizhenkun/codebase/RXNRECer/data/structure/pdb/ncbi/")

    for path in [tempdir, slurm_dir, log_dir, pdb_output_dir, target_pdb_dir]:
        logging.info(f"创建文件夹: {path}")
        path.mkdir(parents=True, exist_ok=True)

    # 读取FASTA文件并过滤已存在PDB文件对应的序列
    records_df = read_fasta2df(fasta_path)
    existing_pdb_files = load_existing_pdb_files(target_pdb_dir)
    records_df = records_df[~records_df.id.isin(existing_pdb_files)].reset_index(drop=True)
    records_df['seqlen'] = records_df.sequence.str.len()
    records_df['sequence'] = records_df.sequence.str.upper().str.strip()
    records_df = records_df.sort_values("seqlen").reset_index(drop=True)

    # 筛选出可处理的序列（≤2600 bp），超过 2600 bp 的给出警告
    too_long = records_df[records_df["seqlen"] > 2600]
    if not too_long.empty:
        logging.warning("以下序列长度超过2600，未处理，请手动处理：")
        for _, row in too_long.iterrows():
            logging.warning(f"  ID: {row['id']}, Length: {row['seqlen']}")
    records_df = records_df[records_df["seqlen"] <= 2600].reset_index(drop=True)

    # 拆分FASTA文件并生成SLURM脚本
    chunk_info = chunk_fasta(records_df, sequences_per_chunk_default=400, outdir=tempdir)
    generate_all_slurm_scripts(chunk_info,
                               input_dir=tempdir,
                               output_dir=pdb_output_dir,
                               log_dir=log_dir,
                               script_dir=slurm_dir)

if __name__ == "__main__":
    step_by_step_run()