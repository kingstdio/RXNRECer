"""
数据摄取（ingestion）编排层：组合 data.sources 与本地 I/O。

负责组织从远端拉取数据并持久化为 feather/tsv 等格式。
"""

from __future__ import annotations

import os
import pandas as pd
from .sources import uniprot as uniprot_src
from .sources import rhea as rhea_src
from rxnrecer.utils.file_utils import ensure_dir


def build_uniprot_rhea_relation(out_feather: str) -> pd.DataFrame:
    """获取 UniProt→Rhea 关系并写入 feather。"""
    df = uniprot_src.fetch_uniprot_rhea_relation()
    ensure_dir(os.path.dirname(out_feather))
    df.to_feather(out_feather)
    return df


def build_rhea_reactions_smiles(out_feather: str) -> pd.DataFrame:
    """获取 Rhea SMILES 表并写入 feather。"""
    df = rhea_src.fetch_rhea_reactions_smiles()
    ensure_dir(os.path.dirname(out_feather))
    df.to_feather(out_feather)
    return df


