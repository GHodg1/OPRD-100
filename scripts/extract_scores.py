"""
Script to extract key metrics from scoring results.
"""
import argparse
import json
import pandas as pd
import numpy as np
import os


def extract_lenient_scores(results_dir, output_file):
    """Flatten a lenient scoring run's ``lenient_scores.json`` into a leaderboard record.

    The lenient scorer (``run_scoring.py --lenient``) already writes
    ``lenient_scores.json`` with all aggregate metrics plus a ``per_type_breakdown``.
    This lifts the per-type ``combined_mean`` values to top-level ``scheme_mean`` /
    ``table_mean`` / ``exp_mean`` so the lenient leaderboard can render them the same
    way the strict leaderboard does.
    """
    print(f"Extracting lenient scores from: {results_dir}")
    with open(os.path.join(results_dir, "lenient_scores.json")) as f:
        raw = json.load(f)

    breakdown = raw.get("per_type_breakdown", {})

    def _type_combined(t):
        return float(breakdown.get(t, {}).get("combined_mean", 0.0))

    scores = {
        "mode": "lenient",
        "similarity_method": raw.get("similarity_method", "tanimoto"),
        # Headline
        "combined_mean": float(raw.get("combined_mean", 0.0)),           # per-match quality
        "combined_mean_covered": float(raw.get("combined_mean_covered", 0.0)),  # x coverage
        "coverage": float(raw.get("coverage", 0.0)),
        # Per-type combined (from the breakdown)
        "exp_mean": _type_combined("Experimental"),
        "table_mean": _type_combined("Table"),
        "scheme_mean": _type_combined("Scheme"),
        # Individual metric means (over matched pairs)
        "reagent_name_mean": float(raw.get("reagent_name_mean", 0.0)),
        "reagent_amount_mean": float(raw.get("reagent_amount_mean", 0.0)),
        "reaction_smiles_mean": float(raw.get("reaction_smiles_mean", 0.0)),
        "solvent_mean": float(raw.get("solvent_mean", 0.0)),
        "time_mean": float(raw.get("time_mean", 0.0)),
        "temperature_mean": float(raw.get("temperature_mean", 0.0)),
        "yield_mean": float(raw.get("yield_mean", 0.0)),
        "reaction_steps_mean": float(raw.get("reaction_steps_mean", 0.0)),
        # Counts
        "num_matched": int(raw.get("num_matched", 0)),
        "num_gold": int(raw.get("num_gold", 0)),
        "num_extracted": int(raw.get("num_extracted", 0)),
    }

    with open(output_file, "w") as f:
        json.dump(scores, f, indent=2)

    print(f"Scores saved to: {output_file}")
    print(f"Combined mean (matched)   : {scores['combined_mean']:.3f}")
    print(f"Combined x coverage       : {scores['combined_mean_covered']:.3f}")
    print(f"Coverage                  : {scores['num_matched']}/{scores['num_gold']} "
          f"({100 * scores['coverage']:.1f}%)")
    return scores


def extract_scores(results_dir, output_file):
    """Extract key metrics from scoring results."""
    
    print(f"Extracting scores from: {results_dir}")
    
    # Load results CSV files
    exp_results = pd.read_csv(os.path.join(results_dir, "experimental_results.csv"))
    table_results = pd.read_csv(os.path.join(results_dir, "table_results.csv"))
    scheme_results = pd.read_csv(os.path.join(results_dir, "scheme_results.csv"))
    combined_results = pd.read_csv(os.path.join(results_dir, "combined_results.csv"))
    
    # Extract key statistics
    scores = {
        # Overall combined scores
        "combined_mean": float(combined_results['total_similarity'].mean()),
        "combined_median": float(combined_results['total_similarity'].median()),
        "combined_std": float(combined_results['total_similarity'].std()),
        
        # Experimental scores
        "exp_mean": float(exp_results['total_similarity'].mean()) if len(exp_results) > 0 else 0.0,
        "exp_median": float(exp_results['total_similarity'].median()) if len(exp_results) > 0 else 0.0,
        "exp_std": float(exp_results['total_similarity'].std()) if len(exp_results) > 0 else 0.0,
        
        # Table scores
        "table_mean": float(table_results['total_similarity'].mean()) if len(table_results) > 0 else 0.0,
        "table_median": float(table_results['total_similarity'].median()) if len(table_results) > 0 else 0.0,
        "table_std": float(table_results['total_similarity'].std()) if len(table_results) > 0 else 0.0,
        
        # Scheme scores
        "scheme_mean": float(scheme_results['total_similarity'].mean()) if len(scheme_results) > 0 else 0.0,
        "scheme_median": float(scheme_results['total_similarity'].median()) if len(scheme_results) > 0 else 0.0,
        "scheme_std": float(scheme_results['total_similarity'].std()) if len(scheme_results) > 0 else 0.0,
        
        # Individual metric means (combined)
        "reagent_name_mean": float(combined_results['reagent_name_similarity'].mean()),
        "reagent_amount_mean": float(combined_results['reagent_amount_similarity'].mean()),
        "reaction_smiles_mean": float(combined_results['reaction_smiles_similarity'].mean()),
        "solvent_mean": float(combined_results['solvent_similarity'].mean()),
        "time_mean": float(combined_results['time_similarity'].mean()),
        "temperature_mean": float(combined_results['temperature_similarity'].mean()),
        "yield_mean": float(combined_results['yield_similarity'].mean()),
        "reaction_steps_mean": float(combined_results['reaction_steps_similarity'].mean()),
        
        # Counts
        "total_reactions": len(combined_results),
        "num_experimental": len(exp_results),
        "num_table": len(table_results),
        "num_scheme": len(scheme_results)
    }
    
    # Save scores
    with open(output_file, 'w') as f:
        json.dump(scores, f, indent=2)
    
    print(f"Scores saved to: {output_file}")
    print(f"Combined mean similarity: {scores['combined_mean']:.3f}")
    print(f"Combined median similarity: {scores['combined_median']:.3f}")
    
    return scores


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract scores from results')
    parser.add_argument('--results-dir', required=True, help='Directory containing results CSV files')
    parser.add_argument('--output', required=True, help='Output JSON file path')
    parser.add_argument(
        '--mode', default='strict', choices=['strict', 'lenient'],
        help="Which scoring run to extract: 'strict' (reads *_results.csv) or "
             "'lenient' (reads lenient_scores.json).",
    )
    args = parser.parse_args()

    if args.mode == 'lenient':
        extract_lenient_scores(args.results_dir, args.output)
    else:
        extract_scores(args.results_dir, args.output)
