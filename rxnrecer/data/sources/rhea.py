"""
Rhea 数据源客户端。

提供获取 Rhea 反应等表格数据的轻量封装。
"""

from __future__ import annotations

import requests
import pandas as pd
from io import StringIO
from requests.adapters import HTTPAdapter, Retry

DEFAULT_TIMEOUT: int = 60
DEFAULT_RETRIES: int = 5
BACKOFF_FACTOR: float = 0.25


#region: HTTP 会话
def get_http_session() -> requests.Session:
    session = requests.Session()
    retries = Retry(total=DEFAULT_RETRIES, backoff_factor=BACKOFF_FACTOR, status_forcelist=[500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session
#endregion

DEFAULT_TIMEOUT = 60


def fetch_rhea_reactions_smiles(timeout: int = DEFAULT_TIMEOUT) -> pd.DataFrame:
    """
    获取 Rhea reaction SMILES 表（TSV）并解析为 DataFrame。
    """
    url = "https://ftp.expasy.org/databases/rhea/tsv/rhea-reaction-smiles.tsv"
    session = get_http_session()
    resp = session.get(url, timeout=timeout)
    resp.raise_for_status()
    return pd.read_csv(StringIO(resp.text), sep="\t")


def get_rhea_reactions():
    #rhea_web_reactions
    rhea_reactions=pd.read_csv('https://www.rhea-db.org/rhea/?query=&columns=rhea-id,equation,chebi-id,ec&format=tsv', sep='\t')
    return rhea_reactions

