# OPRD-100: Ground-Truth Reaction Dataset and Validation Similarity Metrics

OPRD-100 is a ground-truth dataset of reaction data curated from 100 papers in Organic Process Research & Development (OPRD). This repository also provides code and example notebooks to compute and visualize similarity metrics for validating automated data extraction on the same literature.

## Key features
- Curated ground-truth reaction data from 100 OPRD papers.
- Re-extraction validation subset and scoring pipeline.
- Multiple similarity metrics for comparing structured reaction data.
- Example notebooks with code and plots illustrating the metrics.

# Compute Scoring metrics for your dataset!

There are two ways to evaluate your extraction method:

### Option 1: Use the Jupyter Notebook (Manual)
- Use the `example_scoring.ipynb` notebook in the notebooks folder
- Change the filepath to your own dataset
- Run all cells to generate scores and visualizations

### Option 2: Submit to the Leaderboard (Automated)
Submit your extracted data via pull request for automated scoring and public leaderboard inclusion:

1. **Extract reactions** from OPRD-100 papers using your method
2. **Format your data** following the [submission template](data/submissions/README.md)
3. **Create a pull request** with your submission file in `data/submissions/`
4. **Automated scoring** runs and adds your results to the leaderboard below

See detailed instructions in [data/submissions/README.md](data/submissions/README.md).

## 🏆 Leaderboard

Submissions are automatically scored against the OPRD-100 validation dataset. Scores represent similarity metrics between extracted and ground-truth data. The submission by GHodg1 is a test submission of the full OPRD-100 dataset scored against itself. Reaction location classifications are simplified in the scoring workflow and so are different from those in the manuscript. Only reactions from Scheme/Table/Experimental are considered, so some reactions are missing from the full 3.8K reactions in OPRD-100.

| Rank | Submitter | Combined | Experimental | Table | Scheme | Reactions | Date | Details |
|------|-----------|----------|--------------|-------|--------|-----------|------|----------|
| 1 | GHodg1 | 1.000 | 1.000 | 1.000 | 1.000 | 3537 | 2026-01-10 | [PR #6](../../pull/6) |
| 2 | test_workflow_submission | 0.933 | N/A | N/A | 0.933 | 5 | 2026-01-03 | [PR #5](../../pull/5) |

**Metrics explanation:**
- **Combined**: Average similarity across all reaction types (higher is better, max 1.0)
- **Experimental/Table/Scheme**: Average similarity for each data location type
- **Reactions**: Total number of reactions scored
## Reproducibility tips
- Keep raw, ground-truth, and validation JSONs immutable; write derived artifacts to a results/ folder.

## Similarity metrics
Unless stated otherwise, scores are binary (1/0) or Jaccard-style similarities computed on the exact reported text strings.

1) Reaction entry count
  - Equality of the number of reactions per location (for diagnostics; not in total_similarity).

2) SMILES similarity
  - String-based comparison of reported SMILES (exact match or set overlap where applicable).

3) Reaction step number equality
  - Binary equality of step counts.

4) Yield data similarity
  - Comparison of reported yields and their associations with products.

5) Reagent similarity
  - Product of two sub-scores:
    - Name match score
    - Amount-of-substance score

6) Solvent similarity
  - Product of two sub-scores:
    - Name match score
    - Amount-of-substance score

7) Reaction time similarity
  - 1 minus the normalized difference between reported times.

8) Temperature similarity
  - Binary comparison of reported min/max temperatures.

total_similarity = average of metrics 2–8 for a candidate reaction pair.

## Pairing strategy
- Optimal pairing via Hungarian algorithm using total_similarity as the cost/score.
- When counts differ, lowest-scoring unmatched reactions are excluded (affected 6 reactions in total during human validation).




## 🤖 Lenient Leaderboard (AI / OCSR Extraction)

Content-based scoring for automated extraction (AI / OCSR). Reactions are matched by chemical/field similarity (Hungarian assignment), ignoring `Location` labels, and the reaction-SMILES metric uses graded Tanimoto by default. Ranked by **Score** (per-match quality x coverage).

> The human re-extraction row (if present) marks the **maximum achievable** per-match quality under this method — a cross-referenceable target for AI systems.

| Rank | Submitter | Score | Quality | Coverage | Experimental | Table | Scheme | Matched | Date | Details |
|------|-----------|-------|---------|----------|--------------|-------|--------|---------|------|----------|
| 1 | Human re-extraction (reference) | 0.154 | 0.893 | 17.3% | 0.822 | 0.901 | 0.872 | 190/1100 | 2026-07-10 | [PR #0](../../pull/0) |

**Metrics explanation:**
- **Score**: Combined per-match similarity × coverage (the headline ranking metric)
- **Quality**: Mean similarity over matched reactions (the human-achievable ceiling)
- **Coverage**: Fraction of in-scope gold reactions that were matched
- **Experimental/Table/Scheme**: Per-type combined similarity over matched reactions
- **Matched**: Matched reactions / gold reactions in the attempted papers


## Contributing
- Open an issue to propose features, metrics, or bug fixes.

## Citation
If you use OPRD-100 or the validation metrics in a publication, please cite this repository. A BibTeX entry can be added here once the persistent identifier (e.g., DOI) is available.

## Contact
- For questions or feedback, please open an issue in this repository.