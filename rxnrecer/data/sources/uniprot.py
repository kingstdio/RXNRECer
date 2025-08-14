"""
UniProt 数据源客户端。

本模块提供访问 UniProt REST API 的轻量封装，并将响应解析为 Pandas DataFrame。
仅负责“取数与解析”，不涉及本地落盘路径等工程细节。
"""

from __future__ import annotations

from typing import List
import requests
import pandas as pd
from io import StringIO
from requests.adapters import HTTPAdapter, Retry
from typing import Generator, Optional, Tuple
import re
from tqdm import tqdm

# 默认请求参数
DEFAULT_TIMEOUT: int = 30
DEFAULT_RETRIES: int = 5
BACKOFF_FACTOR: float = 0.25


#region: HTTP 会话
def get_http_session() -> requests.Session:
    """构建带重试策略的 requests.Session。"""
    session = requests.Session()
    retries = Retry(total=DEFAULT_RETRIES, backoff_factor=BACKOFF_FACTOR, status_forcelist=[500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session
#endregion



def get_next_link(headers) -> Optional[str]:
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    link_val = headers.get("Link")
    if link_val:
        match = re_next_link.match(link_val)
        if match:
            return match.group(1)
    return None
        
def get_batch(batch_url: str, session: requests.Session) -> Generator[Tuple[requests.Response, str], None, None]:
    while batch_url:
        response = session.get(batch_url, timeout=DEFAULT_TIMEOUT)
        response.raise_for_status()
        total = response.headers.get("x-total-results", "0")
        yield response, total
        batch_url = get_next_link(response.headers)

def get_batch_data_from_uniprot_rest_api(url: str) -> List[List[str]]:
    session = requests.Session()
    retries = Retry(total=DEFAULT_RETRIES, backoff_factor=BACKOFF_FACTOR, status_forcelist=[500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    try:
        rows: List[str] = []
        for batch, _ in tqdm(get_batch(url, session)):
            batch_records = batch.text.splitlines()[1:]
            rows.extend(batch_records)
        return [item.split('\t') for item in rows if item.strip()]
    except requests.RequestException as e:
        raise RuntimeError(f"Failed to fetch data from UniProt: {e}")


def get_one_record_from_uniprot_rest_api(url: str) -> str:
    resp = requests.get(url, timeout=DEFAULT_TIMEOUT)
    resp.raise_for_status()
    return resp.text




# Function to fetch UniProt data for multiple IDs at once with batch processing
def get_uniprot_records_by_ids(ids: List[str], batch_size: int = 40) -> pd.DataFrame:
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    session = requests.Session()
    retries = Retry(total=DEFAULT_RETRIES, backoff_factor=BACKOFF_FACTOR, status_forcelist=[500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    all_results: List[str] = []
    for i in tqdm(range(0, len(ids), batch_size)):
        batch_ids = ids[i:i + batch_size]
        ids_query = " OR accession:".join(batch_ids)
        query = f"accession:{ids_query}"
        api_url = f"{base_url}?query={query}&format=tsv&fields=accession,reviewed,protein_name,gene_names,gene_oln,organism_name,length,ec,rhea,sequence&compressed=false"
        try:
            response = session.get(api_url, timeout=DEFAULT_TIMEOUT)
            response.raise_for_status()
            batch_result = response.text.strip().split('\n')[1:]
            all_results.extend(batch_result)
        except requests.RequestException as e:
            print(f"Error fetching batch {i // batch_size + 1}: {e}")
    rows = [item.split('\t') for item in all_results if item.strip()]
    columns = ['uniprot_id', 'reviewed', 'protein_name', 'gene_names', 'gene_oln', 'organism_name', 'length', 'ec', 'reaction_id', 'seq']
    all_results_df = pd.DataFrame(rows, columns=columns)
    all_results_df = all_results_df.dropna(how='all')
    return all_results_df


def fetch_uniprot_rhea_relation(size: int = 500, timeout: int = DEFAULT_TIMEOUT) -> pd.DataFrame:
    """
    拉取 UniProt → Rhea 关系表（TSV），返回 DataFrame：['uniprot_id','ec','reaction_id']。
    """
    api_url = (
        "https://rest.uniprot.org/uniprotkb/search?"
        "query=reviewed=true&format=tsv&fields=accession,ec,rhea&size=" + str(size)
    )
    session = get_http_session()
    response = session.get(api_url, timeout=timeout)
    response.raise_for_status()
    df = pd.read_csv(StringIO(response.text), sep="\t")
    df.columns = ["uniprot_id", "ec", "reaction_id"]
    df = df.fillna("-")

    def _strip_comp(s: str) -> str:
        if not isinstance(s, str) or not s:
            return "-"
        return ";".join([item for item in s.split(" ") if "RHEA-COMP" not in item]) or "-"

    df["reaction_id"] = df["reaction_id"].astype(str).map(_strip_comp)
    return df

def read_snapshot_tsv(path: str) -> pd.DataFrame:
    """Read a UniProt snapshot TSV produced by upstream extraction scripts."""
    return pd.read_csv(path, sep="\t", header=0)


