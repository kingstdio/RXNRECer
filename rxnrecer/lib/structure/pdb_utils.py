"""
PDB structure selection and download utilities.

Functions for selecting the best PDB structure for a given UniProt ID
based on resolution, experimental method, and chain ID.
"""

from __future__ import annotations

import requests
from typing import Optional, Tuple, Dict, Any

import pandas as pd


def json_to_dataframe(json_data: Dict) -> pd.DataFrame:
    """
    Convert PDB JSON data to DataFrame with UniProt ID column.

    Args:
        json_data: JSON data with UniProt ID as key

    Returns:
        DataFrame with 'uniprot_id' as first column
    """
    uniprot_id = list(json_data.keys())[0]

    for entry in json_data[uniprot_id]:
        entry["uniprot_id"] = uniprot_id

    df = pd.DataFrame(json_data[uniprot_id])

    cols = ["uniprot_id"] + [col for col in df.columns if col != "uniprot_id"]
    return df[cols]


def select_best_pdb(df: pd.DataFrame) -> Tuple[str, str]:
    """
    Select the best PDB structure based on selection criteria.

    Selection rules:
    1. Lower resolution is better
    2. X-ray method is preferred over other methods
    3. Earlier chain ID (A before B, etc.)

    Args:
        df: DataFrame with PDB information including resolution, experimental_method, chain_id

    Returns:
        Tuple of (pdb_id, chain_id)
    """
    df_sorted = df.sort_values(
        by=["resolution", "experimental_method"],
        ascending=[True, False],
    )

    df_sorted = df_sorted.sort_values("chain_id", ascending=True)

    best_row = df_sorted.iloc[0]
    return best_row["pdb_id"], best_row["chain_id"]


def download_pdb(pdb_id: str, save_path: Optional[str] = None) -> Optional[str]:
    """
    Download PDB file from RCSB.

    Args:
        pdb_id: 4-character PDB ID
        save_path: Optional path to save the PDB file

    Returns:
        PDB file content or None if download fails
    """
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url)
        response.raise_for_status()

        if save_path:
            with open(save_path, "w") as f:
                f.write(response.text)

        return response.text
    except Exception as e:
        print(f"Download failed: {e}")
        return None


def parse_pdb_info(pdb_json: Dict) -> Tuple[Dict[str, Any], list]:
    """
    Parse PDB JSON data and extract key information.

    Args:
        pdb_json: JSON data from RCSB REST API

    Returns:
        Tuple of (main_info dict, ligands list)
    """
    main_info = {
        "pdb_id": pdb_json.get("entry", {}).get("id", ""),
        "resolution": pdb_json.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
        "method": pdb_json.get("exptl", [{}])[0].get("method", "").title(),
        "journal": f"{pdb_json.get('citation', [{}])[0].get('journal_abbrev', '')} "
        f"({pdb_json.get('citation', [{}])[0].get('year', '')})",
        "ref_doi": pdb_json.get("citation", [{}])[0].get("pdbx_database_id_doi", ""),
    }

    ligands = []
    if "rcsb_binding_affinity" in pdb_json:
        ligands = [
            {
                "chemical_id": lig["chemical_id"],
                "chemical_name": lig.get("chemical_name", ""),
                "formula": lig.get("formula", ""),
            }
            for lig in pdb_json["rcsb_binding_affinity"]
        ]

    return main_info, ligands


def get_best_pdb(uniprot_id: str, target_dir: str = "") -> Optional[str]:
    """
    Get the best PDB structure for a given UniProt ID.

    Args:
        uniprot_id: UniProt accession ID
        target_dir: Directory to save the downloaded PDB file

    Returns:
        PDB ID if successful, None otherwise
    """
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
    response = requests.get(url).json()

    if not response:
        raise ValueError(f"No PDB found for UniProt ID: {uniprot_id}")

    response = json_to_dataframe(response)

    best_pdb, best_chain = select_best_pdb(response)
    print(f"Best PDB ID: {uniprot_id} -> {best_pdb}")

    try:
        pdb_api = f"https://data.rcsb.org/rest/v1/core/entry/{best_pdb}"
        pdb_info = requests.get(pdb_api).json()
        main_info, ligands = parse_pdb_info(pdb_info)
        print(main_info, ligands)
    except Exception as e:
        print(f"Error processing PDB RCSB_{uniprot_id}_{best_pdb}: {e}")

    download_pdb(best_pdb, f"{target_dir}/RCSB_{uniprot_id}_{best_pdb}.pdb")

    return best_pdb
