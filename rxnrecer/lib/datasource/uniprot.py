"""
UniProt 数据源客户端。

本模块提供访问 UniProt REST API 的轻量封装，并将响应解析为 Pandas DataFrame。
仅负责“取数与解析”，不涉及本地落盘路径等工程细节。
"""

from __future__ import annotations

from io import StringIO
from typing import List
import requests
import pandas as pd
from requests.adapters import HTTPAdapter, Retry
from typing import Generator, Optional, Tuple
import re
from tqdm import tqdm

# 默认请求参数
DEFAULT_TIMEOUT: int = 30
DEFAULT_RETRIES: int = 5
BACKOFF_FACTOR: float = 0.25
DEFAULT_FIELDS = (
    "accession,reviewed,id,protein_name,gene_names,organism_name,length,"
    "ec,rhea,protein_evidence,cc_function,sequence"
)


#region: HTTP 会话
def get_http_session() -> requests.Session:
    """构建带重试策略的 requests.Session。"""
    session = requests.Session()
    retries = Retry(
        total=DEFAULT_RETRIES,
        backoff_factor=BACKOFF_FACTOR,
        status_forcelist=[429, 500, 502, 503, 504],
    )
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


def parse_tsv_response(tsv_text: str) -> pd.DataFrame:
    """将 UniProt TSV 响应解析为 DataFrame。"""
    if not tsv_text.strip():
        return pd.DataFrame()
    return pd.read_csv(StringIO(tsv_text), sep="\t")


class UniProtClient:
    """带分页和重试的 UniProt REST 客户端。"""

    def __init__(self, session: Optional[requests.Session] = None, timeout: int = DEFAULT_TIMEOUT):
        self.session = session or get_http_session()
        self.timeout = timeout
        self.base_url = "https://rest.uniprot.org/uniprotkb/search"

    def fetch_records(
        self,
        query: Optional[str] = None,
        ids: Optional[List[str]] = None,
        fields: Optional[str] = None,
        format: str = "tsv",
        batch_size: int = 500,
    ) -> pd.DataFrame:
        """
        获取 UniProt 记录。

        - `ids` 模式下会按批次组装 accession 查询
        - `query` 模式下直接执行查询
        - 自动处理分页，并去重
        """
        if not query and not ids:
            raise ValueError("Either query or ids must be provided")

        if format != "tsv":
            raise ValueError("Only TSV format is supported")

        fields = fields or DEFAULT_FIELDS
        result_frames: List[pd.DataFrame] = []

        if ids:
            for start in tqdm(range(0, len(ids), batch_size), desc="Fetching UniProt records by IDs"):
                batch_ids = ids[start:start + batch_size]
                batch_query = " OR ".join(f"accession:{item}" for item in batch_ids)
                result_frames.extend(self._request_query_frames(batch_query, fields=fields, format=format))
        else:
            result_frames.extend(self._request_query_frames(query=query, fields=fields, format=format))

        if not result_frames:
            return pd.DataFrame()

        final_df = pd.concat(result_frames, ignore_index=True)
        entry_col = "Entry" if "Entry" in final_df.columns else "accession" if "accession" in final_df.columns else None
        if entry_col:
            final_df = final_df.drop_duplicates(subset=[entry_col])
            if entry_col == "Entry":
                final_df = final_df.rename(columns={"Entry": "uniprot_id"})
            elif entry_col == "accession":
                final_df = final_df.rename(columns={"accession": "uniprot_id"})
        return final_df

    def _request_query_frames(self, query: str, fields: str, format: str) -> List[pd.DataFrame]:
        params = {"query": query, "fields": fields, "format": format, "size": 500}
        next_url = self.base_url
        frames: List[pd.DataFrame] = []
        while next_url:
            response = self.session.get(
                next_url,
                params=params if next_url == self.base_url else None,
                timeout=self.timeout,
            )
            response.raise_for_status()
            frame = parse_tsv_response(response.text)
            if not frame.empty:
                frames.append(frame)
            next_url = get_next_link(response.headers)
        return frames
        
def get_batch(batch_url: str, session: requests.Session) -> Generator[Tuple[requests.Response, str], None, None]:
    while batch_url:
        response = session.get(batch_url, timeout=DEFAULT_TIMEOUT)
        response.raise_for_status()
        total = response.headers.get("x-total-results", "0")
        yield response, total
        batch_url = get_next_link(response.headers)

def get_batch_data_from_uniprot_rest_api(url: str) -> List[List[str]]:
    session = get_http_session()
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
    client = UniProtClient()
    fields = "accession,reviewed,protein_name,gene_names,gene_oln,organism_name,length,ec,rhea,sequence"
    result = client.fetch_records(ids=ids, batch_size=batch_size, fields=fields)
    rename_map = {
        "Reviewed": "reviewed",
        "Protein names": "protein_name",
        "Gene Names": "gene_names",
        "Gene Names (ordered locus)": "gene_oln",
        "Organism": "organism_name",
        "Length": "length",
        "EC number": "ec",
        "Rhea": "reaction_id",
        "Sequence": "seq",
    }
    return result.rename(columns=rename_map)


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

