# OPRD-100: Ground-Truth Reaction Dataset and Validation Similarity Metrics

OPRD-100 is a ground-truth dataset of reaction data curated from 100 papers in Organic Process Research & Development (OPRD). This repository also provides code and example notebooks to compute and visualize similarity metrics for validating automated data extraction on the same literature.

## Key features
- Curated ground-truth reaction data from 100 OPRD papers.
- Re-extraction validation subset and scoring pipeline.
- Multiple similarity metrics for comparing structured reaction data.
- Example notebooks with code and plots illustrating the metrics.

# Scoring OPRD-100 extractions

OPRD-100 has two scoring workflows with deliberately different uses:

| Workflow | Pairing | Use |
|---|---|---|
| **Strict validation** | Exact `(Reference, Location.Type, Location.Num)` | Reproduce the paper's human re-extraction validation only |
| **Lenient automated scoring** | Content-based Hungarian matching within each paper and location type | Score LLM, OCSR, and rule-based extraction systems |

Do not use strict scoring to rank automated extraction systems. Automated systems can
recover the correct reaction while assigning a different Scheme, Table, Experimental,
or location number; strict matching would discard those reactions before their content
is assessed. The automated leaderboard therefore reports the lenient score.

There are two ways to run the supported scoring workflows:

### Option 1: Use the Jupyter Notebook (Manual)
- Use [`notebooks/example_scoring.ipynb`](notebooks/example_scoring.ipynb), which contains
  both workflows and explains when each one is appropriate.
- The strict section is fixed to the bundled human validation set from the paper.
- Set `automated_extraction_path` to your extracted reaction JSON to run the lenient section.

### Option 2: Submit to the Leaderboard (Automated)
Submit your extracted data via pull request for automated scoring and public leaderboard inclusion:

1. **Extract reactions** from OPRD-100 papers using your method
2. **Format your data** following the [submission template](data/submissions/README.md)
3. **Create a pull request** with your submission file in `data/submissions/`
4. **Automated scoring** runs with the **lenient (content-based)** method and adds your
   results to the leaderboard below

See detailed instructions in [data/submissions/README.md](data/submissions/README.md).

## 🏆 Leaderboard — Automated Extraction (AI / OCSR)

This is the primary OPRD-100 leaderboard. Automated extraction pipelines (LLM, OCSR, rule-based) are scored **only with the lenient, content-based method**: each predicted reaction is paired to the ground truth by chemical and field similarity (Hungarian assignment) rather than by exact location labels, because automated tools rarely reproduce a paper's `(Scheme/Table/Experimental, number)` bookkeeping. The reaction-SMILES metric uses exact InChIKey-set overlap (the same metric as the strict validation scorer), and stereochemistry is kept.

Entries are ranked by **Score = Quality × Coverage**, rewarding both accurate chemistry and completeness. The `Human re-extraction (reference)` row is the **ceiling** — the quality a careful human achieves under this method, and the target AI systems should aim for.

| Rank | Submitter | Score | Quality | Coverage | Experimental | Table | Scheme | Matched | Date | Details |
|------|-----------|-------|---------|----------|--------------|-------|--------|---------|------|----------|
| 1 | Human re-extraction (reference) | 0.857 | 0.884 | 96.9% | 0.817 | 0.892 | 0.858 | 190/196 | 2026-07-13 | [PR #0](../../pull/0) |

**How to read the scores:**
- **Score**: headline ranking metric, Quality × Coverage (max 1.0)
- **Quality**: mean per-reaction similarity over matched reactions (the human-achievable ceiling)
- **Coverage**: fraction of in-scope ground-truth reactions that were matched
- **Experimental/Table/Scheme**: Quality broken down by ground-truth source location
- **Matched**: matched reactions / ground-truth reactions in the attempted papers
## 🔬 Validation — Human Re-extraction (Strict Scoring)

This table is **not** a competitive leaderboard — it documents the strict validation that establishes OPRD-100 as a faithful ground-truth source. A human independently re-extracted reactions, scored with the **strict** method (reactions matched by exact `(Reference, Type, Num)` location key, stereochemistry required, binary InChIKey SMILES matching). High scores here show OPRD-100 is internally consistent and reproducible.

| Rank | Submitter | Combined | Experimental | Table | Scheme | Reactions | Date | Details |
|------|-----------|----------|--------------|-------|--------|-----------|------|----------|
| 1 | Human re-extraction (validation) | 0.882 | 0.817 | 0.892 | 0.847 | 189 | 2026-07-10 | [PR #0](../../pull/0) |

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

2) Reaction-SMILES similarity
  - By default, molecules on each side of the reaction are converted to InChIKeys and
    compared by set overlap. The optional Tanimoto mode uses Morgan fingerprints for
    graded molecular similarity.

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

## Contributing
- Open an issue to propose features, metrics, or bug fixes.

## Citation
If you use OPRD-100 or the validation metrics in a publication, please cite this repository. A BibTeX entry can be added here once the persistent identifier (e.g., DOI) is available.

## Contact
- For questions or feedback, please open an issue in this repository.