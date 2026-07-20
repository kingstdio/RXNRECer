"""Utilities for categorizing evidence annotations from UniProt text fields."""

from __future__ import annotations

import re

import pandas as pd


EVIDENCE_PRIORITY = {"EXP": 4, "CUR": 3, "SIM": 2, "AUT": 1, "UNK": 0}


def categorize_evidence(text: str | float | None) -> str:
    """Map UniProt evidence text to a compact evidence level."""
    if not text or pd.isna(text):
        return "UNK"

    tags = [f"ECO:{eco}" for eco in re.findall(r"ECO:(\d{7})", str(text))]
    if "PubMed" in str(text):
        tags.append("PubMed")
    if "Author statement" in str(text):
        tags.append("Author statement")
    if not tags:
        return "UNK"

    categories = set()
    for tag in tags:
        if tag in {"ECO:0000269", "ECO:0000314", "ECO:0000353", "ECO:0000315", "PubMed"}:
            categories.add("EXP")
        elif tag in {"ECO:0000305", "ECO:0000303", "Author statement"}:
            categories.add("CUR")
        elif tag in {"ECO:0000250", "ECO:0000255", "ECO:0000266"}:
            categories.add("SIM")
        elif tag in {"ECO:0000256", "ECO:0000259", "ECO:0007744"}:
            categories.add("AUT")

    for level in ("EXP", "CUR", "SIM", "AUT"):
        if level in categories:
            return level
    return "UNK"


def merge_evidence_levels(*levels: str) -> str:
    """Return the strongest evidence level among inputs."""
    best = "UNK"
    for level in levels:
        if EVIDENCE_PRIORITY.get(level, 0) > EVIDENCE_PRIORITY[best]:
            best = level
    return best


def annotate_evidence_levels(
    df: pd.DataFrame,
    function_col: str = "Function [CC]",
    catalytic_col: str = "Catalytic activity",
    output_col: str = "evidence_level",
) -> pd.DataFrame:
    """Annotate a UniProt-like table with a merged evidence level."""
    if function_col not in df.columns:
        raise KeyError(f"Missing function column: {function_col}")
    if catalytic_col not in df.columns:
        raise KeyError(f"Missing catalytic column: {catalytic_col}")

    result = df.copy()
    result["_evidence_func"] = result[function_col].apply(categorize_evidence)
    result["_evidence_cata"] = result[catalytic_col].apply(categorize_evidence)
    result[output_col] = result.apply(
        lambda row: merge_evidence_levels(row["_evidence_func"], row["_evidence_cata"]),
        axis=1,
    )
    return result.drop(columns=["_evidence_func", "_evidence_cata"])
