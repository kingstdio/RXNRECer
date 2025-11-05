import os
import sys
import argparse
from typing import Optional
import pandas as pd
from tqdm import tqdm

from rxnrecer.config import config as cfg
from rxnrecer.lib.rxn.Reaction import Reaction


def generate_svgs_from_reactions(limit: Optional[int] = None) -> int:
    """
    Generate molecule SVGs for all reactants and products referenced by reactions
    in cfg.FILE_RHEA_REACTION. Files are saved under cfg.DIR_CPD_SVG.

    Returns the total number of SVG write attempts (unique molecules are deterministically
    named and will be skipped if already exist).
    """
    # Ensure output directory exists
    os.makedirs(cfg.DIR_CPD_SVG, exist_ok=True)

    df = pd.read_feather(cfg.FILE_RHEA_REACTION)

    if limit is not None:
        df = df.head(limit)

    total = 0
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Generating molecule SVGs"):
        try:
            rxn = Reaction(
                rxn_smiles=row.equation_smiles,
                rxn_equation=row.equation,
                rxn_equation_ref_chebi=row.equation_chebi,
                rxn_id=row.reaction_id,
                rxn_ec=row.ec_number if 'ec_number' in row.index else ''
            )
            # Accessing reactants/products triggers SVG generation via Molecule.draw_mol_simple
            total += len(rxn.reactants) + len(rxn.products)
        except Exception as e:
            print(f"Warning: failed to process reaction {getattr(row, 'reaction_id', '?')}: {e}")

    return total


def main():
    parser = argparse.ArgumentParser(
        description="Generate molecule SVGs for all reactions in the reaction table",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process the first N reactions (for testing)",
    )
    args = parser.parse_args()

    print(f"Project root: {cfg.DIR_PROJECT_ROOT}")
    print(f"Reaction table: {cfg.FILE_RHEA_REACTION}")
    print(f"Output SVG dir: {cfg.DIR_CPD_SVG}")

    total = generate_svgs_from_reactions(limit=args.limit)
    print(f"Done. Touched {total} molecules across reactions.")


if __name__ == "__main__":
    sys.exit(main())


