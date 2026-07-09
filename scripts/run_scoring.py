"""
Script to run scoring analysis on submission data.
This replicates the analysis from example_scoring.ipynb programmatically.
"""
import argparse
import json
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

import pandas as pd
import numpy as np
from validation import DataComparer
from plotting import Plot
from security import validate_submission


def run_scoring(submission_file, output_dir, similarity_method="inchikey"):
    """Run scoring analysis on submission data."""
    
    print(f"Loading submission from: {submission_file}")
    print(f"Output directory: {output_dir}")
    print(f"Similarity method: {similarity_method}")
    
    # SECURITY: Validate submission file before processing
    print("Validating submission for security threats...")
    is_valid, error_msg, submission = validate_submission(submission_file)
    
    if not is_valid:
        print(f"❌ SECURITY VALIDATION FAILED: {error_msg}")
        print("Submission rejected due to security concerns.")
        sys.exit(1)
    
    print("✅ Security validation passed")
    
    # Load submission (already validated and loaded)
    submitter_name = submission.get('submitter_name', 'anonymous')
    method_description = submission.get('method_description', 'N/A')
    
    print(f"Submitter: {submitter_name}")
    print(f"Method: {method_description}")
    
    # Create temporary file with submission reactions
    submission_data_path = os.path.join(output_dir, "submission_data.json")
    with open(submission_data_path, 'w') as f:
        json.dump(submission['reactions'], f, indent=2)
    
    # Initialize DataComparer with submission data
    # Need to change to the directory where validation.py expects to be run from
    print("Initializing DataComparer...")
    original_dir = os.getcwd()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(script_dir)
    notebooks_dir = os.path.join(repo_root, 'notebooks')
    
    # Change to notebooks directory temporarily (DataComparer uses relative paths)
    os.chdir(notebooks_dir)
    
    # Make submission path absolute before changing directory
    submission_data_path_abs = os.path.abspath(os.path.join(original_dir, submission_data_path))
    
    dc = DataComparer(submission_data_path_abs, similarity_method=similarity_method)
    
    # Run all the same analyses as in example_scoring.ipynb
    print("Computing comparison scores...")
    exp_results = dc.compute_comparison_scores(dc.val_experimental_data, dc.oprd_experimental_data)
    scheme_results = dc.compute_comparison_scores(dc.val_scheme_data, dc.oprd_scheme_data)
    table_results = dc.compute_comparison_scores(dc.val_table_data, dc.oprd_table_data)
    combined_results = pd.concat([exp_results, scheme_results, table_results])
    
    print(f"Total reactions scored: {len(combined_results)}")
    print(f"  Experimental: {len(exp_results)}")
    print(f"  Scheme: {len(scheme_results)}")
    print(f"  Table: {len(table_results)}")
    
    # Change back to original directory for output
    os.chdir(original_dir)
    output_dir_abs = os.path.abspath(output_dir)
    
    # Generate all plots
    print("Generating plots...")
    generate_all_plots(dc, exp_results, scheme_results, table_results, combined_results, output_dir_abs)
    
    # Save results DataFrames
    print("Saving results...")
    exp_results.to_csv(os.path.join(output_dir_abs, "experimental_results.csv"), index=False)
    scheme_results.to_csv(os.path.join(output_dir_abs, "scheme_results.csv"), index=False)
    table_results.to_csv(os.path.join(output_dir_abs, "table_results.csv"), index=False)
    combined_results.to_csv(os.path.join(output_dir_abs, "combined_results.csv"), index=False)
    
    # Save submission metadata
    metadata = {
        'submitter_name': submitter_name,
        'method_description': method_description,
        'repository_url': submission.get('repository_url', ''),
        'paper_url': submission.get('paper_url', ''),
        'contact_email': submission.get('contact_email', ''),
        'num_reactions': len(combined_results),
        'num_experimental': len(exp_results),
        'num_scheme': len(scheme_results),
        'num_table': len(table_results)
    }
    with open(os.path.join(output_dir_abs, "metadata.json"), 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print("Scoring complete!")
    return exp_results, scheme_results, table_results, combined_results


def generate_all_plots(dc, exp_results, scheme_results, table_results, combined_results, output_dir):
    """Generate all plots from example_scoring.ipynb"""
    
    # Plot 1: Entry count comparison
    print("  - Entry count comparison...")
    oprd_counts_table, val_counts_table = dc.get_entry_counts(dc.oprd_table_data, dc.val_table_data)
    oprd_counts_scheme, val_counts_scheme = dc.get_entry_counts(dc.oprd_scheme_data, dc.val_scheme_data)
    
    oprd_counts_table_values = list(oprd_counts_table.values())
    oprd_counts_table_keys = [i[0] + f" ({i[-1]})" for i in oprd_counts_table.keys()]
    val_counts_table_values = list(val_counts_table.values())
    oprd_counts_scheme_values = list(oprd_counts_scheme.values())
    oprd_counts_scheme_keys = [i[0] + f" ({i[-1]})" for i in oprd_counts_scheme.keys()]
    val_counts_scheme_values = list(val_counts_scheme.values())
    
    shared_title_fontdict = {'fontsize': 16, 'fontweight': 'bold'}
    label_fontdict = {'fontsize': 12}
    tick_fontdict = {'fontsize': 12}
    legend_fontdict = {'size': 15}
    
    plot = Plot(n_rows=1, n_cols=2, figsize=(8, 4), 
                shared_title="Reaction Entry Count Comparison\n", 
                shared_title_fontdict=shared_title_fontdict, 
                label_fontdict=label_fontdict, 
                tick_fontdict=tick_fontdict, 
                legend_fontdict=legend_fontdict)
    
    plot.add_bar(row=0, col=0, categories=oprd_counts_scheme_keys, 
                 values=[oprd_counts_scheme_values, val_counts_scheme_values], 
                 side_by_side=True, title="Scheme", bar_labels=["", ""], 
                 colour=["blue", "orange"], x_label="DOI (Scheme number)", 
                 y_label="Count", alpha=0.4, rotation=45)
    
    plot.add_bar(row=0, col=1, categories=oprd_counts_table_keys, 
                 values=[oprd_counts_table_values, val_counts_table_values], 
                 side_by_side=True, title="Table", bar_labels=["", ""], 
                 colour=["blue", "orange"], x_label="DOI (Table number)", 
                 y_label="", alpha=0.4, rotation=45)
    
    entries = [
        {'label': 'OPRD-100', 'color': 'blue', 'alpha': 0.4},
        {'label': 'Submission', 'color': 'orange', 'alpha': 0.4},
    ]
    
    plot.add_figure_legend(entries, title=None, anchor=(0.5, 0.775), 
                          loc='lower center', frame=True, 
                          frame_facecolor='white', frame_alpha=0.9, 
                          frame_edgecolor='black', ncol=1, 
                          label_fontsize=11, label_fontweight='bold')
    
    plot.plot(savefig=True, filepath=os.path.join(output_dir, "ReactionEntryCountComparison.png"))
    
    # Plot 2: Combined validation score distributions
    print("  - Combined score distributions...")
    plot_score_distributions(combined_results, len(combined_results), 
                            "Combined", output_dir, "CombinedValidationSimilarityScoreDistributions.png",
                            )
    
    # Plot 3: Experimental score distributions
    if len(exp_results) > 0:
        print("  - Experimental score distributions...")
        plot_score_distributions(exp_results, len(exp_results), 
                                "Experimental", output_dir, 
                                "ExperimentalValidationSimilarityScoreDistributions.png",
                                )
    
    # Plot 4: Table score distributions
    if len(table_results) > 0:
        print("  - Table score distributions...")
        plot_score_distributions(table_results, len(table_results), 
                                "Table", output_dir, 
                                "TableValidationSimilarityScoreDistributions.png",
                                )
    
    # Plot 5: Scheme score distributions
    if len(scheme_results) > 0:
        print("  - Scheme score distributions...")
        plot_score_distributions(scheme_results, len(scheme_results), 
                                "Scheme", output_dir, 
                                "SchemeValidationSimilarityScoreDistributions.png",
                                )


def plot_score_distributions(results, num_reactions, category, output_dir, filename,):
    """Plot score distribution histograms for a category of results."""
    
    shared_title_fontdict = {'fontsize': 16, 'fontweight': 'bold'}
    label_fontdict = {'fontsize': 14}
    tick_fontdict = {'fontsize': 13.5}
    
    title = f"Validation Score Distributions ({category} - {num_reactions} Reactions)"
    
    plot = Plot(n_rows=2, n_cols=4, figsize=(10, 6), 
                shared_title=title, 
                shared_title_fontdict=shared_title_fontdict, 
                label_fontdict=label_fontdict, 
                tick_fontdict=tick_fontdict)
    
    # Count binary metrics
    equal_step_counts = results['reaction_steps_similarity'].value_counts().to_dict()
    equal_step_count_values = [equal_step_counts.get(1.0, 0), equal_step_counts.get(0.0, 0)]
    
    equal_temp_counts = results['temperature_similarity'].value_counts().to_dict()
    equal_temp_count_values = [equal_temp_counts.get(1.0, 0), equal_temp_counts.get(0.0, 0)]
    
    # Plot histograms and bar charts
    plot.add_hist(row=0, col=0, data=results['reagent_name_similarity'], bins=20, 
                  title="", x_label="Reagent Name Score", y_label="Frequency", 
                  alpha=0.7, colour='blue', )
    
    plot.add_hist(row=0, col=1, data=results['reagent_amount_similarity'], bins=20, 
                  title="", x_label="Reagent Amount\nScore", y_label="", 
                  alpha=0.7, colour='blue', )
    
    plot.add_bar(row=0, col=2, categories=["Yes", "No"], values=equal_step_count_values, 
                 title="", x_label="Equal reaction step\ncount?", y_label="", 
                 alpha=0.7, colour='blue', )
    
    plot.add_hist(row=0, col=3, data=results['reaction_smiles_similarity'], bins=20, 
                  title="", x_label="Reaction SMILES\nScore", y_label="", 
                  alpha=0.7, colour='blue', )
    
    plot.add_hist(row=1, col=0, data=results['solvent_similarity'], bins=20, 
                  title="", x_label="Solvents Score", y_label="Frequency", 
                  alpha=0.7, colour='blue', )
    
    plot.add_hist(row=1, col=1, data=results['time_similarity'], bins=20, 
                  title="", x_label="Time Score", y_label="", 
                  alpha=0.7, colour='blue', )
    
    plot.add_bar(row=1, col=2, categories=["Yes", "No"], values=equal_temp_count_values, 
                 title="", x_label="Equal (min, max)\ntemperatures?", y_label="", 
                 alpha=0.7, colour='blue', )
    
    plot.add_hist(row=1, col=3, data=results['yield_similarity'], bins=20, 
                  title="", x_label="Yield Data Score", y_label="", 
                  alpha=0.7, colour='blue', )
    
    plot.plot(savefig=True, filepath=os.path.join(output_dir, filename))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run scoring analysis on submission data')
    parser.add_argument('--submission-file', required=True, help='Path to submission JSON file')
    parser.add_argument('--output-dir', required=True, help='Directory to save results')
    parser.add_argument(
        '--similarity-method', default='inchikey', choices=['inchikey', 'tanimoto'],
        help="Reaction-SMILES similarity: 'inchikey' (default, exact InChIKey-set Jaccard) "
             "or 'tanimoto' (graded Morgan-Tanimoto — partial credit for near-identical molecules)",
    )
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    run_scoring(args.submission_file, args.output_dir, similarity_method=args.similarity_method)
