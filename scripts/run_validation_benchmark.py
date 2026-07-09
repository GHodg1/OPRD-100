"""Benchmark the OPRD-100 validation study under both reaction-SMILES similarity methods.

Scores the human re-extraction dataset (``data/validation_reactions.json``) against the
OPRD-100 ground truth using both the default ``inchikey`` metric and the graded
``tanimoto`` metric, and prints a side-by-side comparison. Use this to confirm the
Tanimoto option still gives strong scores on trusted human re-extraction (and to see how
much it relaxes near-miss molecules).

Run from anywhere:

    python scripts/run_validation_benchmark.py
    python scripts/run_validation_benchmark.py --submission data/validation_reactions.json
"""

import argparse
import os
import sys

import pandas as pd

# DataComparer hardcodes a relative "../data/OPRD-100.json"; run from a repo subdir.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_SCRIPT_DIR)
sys.path.insert(0, os.path.join(_REPO_ROOT, "src"))

from validation import DataComparer  # noqa: E402

_METRICS = [
    "reaction_smiles_similarity",
    "reaction_steps_similarity",
    "yield_similarity",
    "reagent_name_similarity",
    "solvent_similarity",
    "time_similarity",
    "temperature_similarity",
    "total_similarity",
]


def _combined_scores(submission_abs: str, method: str) -> pd.DataFrame:
    dc = DataComparer(submission_abs, similarity_method=method)
    exp = dc.compute_comparison_scores(dc.val_experimental_data, dc.oprd_experimental_data)
    scheme = dc.compute_comparison_scores(dc.val_scheme_data, dc.oprd_scheme_data)
    table = dc.compute_comparison_scores(dc.val_table_data, dc.oprd_table_data)
    return pd.concat([exp, scheme, table], ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--submission",
        default=os.path.join(_REPO_ROOT, "data", "validation_reactions.json"),
        help="Path to the validation/submission JSON (default: data/validation_reactions.json)",
    )
    args = parser.parse_args()
    submission_abs = os.path.abspath(args.submission)

    # DataComparer resolves ../data/OPRD-100.json relative to CWD → run from scripts/.
    os.chdir(_SCRIPT_DIR)

    results: dict[str, pd.DataFrame] = {}
    for method in ("inchikey", "tanimoto"):
        results[method] = _combined_scores(submission_abs, method)

    n = len(results["inchikey"])
    print(f"Validation study: {n} matched reactions\n")
    header = f"{'metric':<28}{'inchikey':>12}{'tanimoto':>12}{'delta':>10}"
    print(header)
    print("-" * len(header))
    for metric in _METRICS:
        a = results["inchikey"][metric].mean()
        b = results["tanimoto"][metric].mean()
        print(f"{metric:<28}{a:>12.4f}{b:>12.4f}{b - a:>+10.4f}")


if __name__ == "__main__":
    main()
