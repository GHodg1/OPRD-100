# Automated Scoring Scripts

This directory contains scripts for automated scoring of submissions via GitHub Actions.

## Scripts

### `run_scoring.py`
Runs the scoring analysis from `example_scoring.ipynb` programmatically on submission data.

**Usage:**
```bash
python scripts/run_scoring.py \
  --submission-file data/submissions/your_submission.json \
  --output-dir results/your_name_timestamp \
  [--similarity-method inchikey|tanimoto]
```

**What it does:**
- Loads submission data and initializes DataComparer
- Computes similarity scores for experimental, scheme, and table reactions
- Generates all plots from the example notebook
- Saves results as CSV files

**Reaction-SMILES similarity method** (`--similarity-method`, default `inchikey`):
- `inchikey` — binary InChIKey-set Jaccard of reactants and products (a molecule matches
  exactly or scores 0). This is the canonical leaderboard metric.
- `tanimoto` — graded Morgan-Tanimoto (radius 2, 2048 bits): near-identical molecules
  (e.g. ethyl vs butyl esters, stereochemistry, salt forms) receive partial credit instead
  of 0. Also drives reaction matching, so the Hungarian assignment picks the closest gold
  reaction. Only the reaction-SMILES sub-score changes; all other field metrics are identical.

### `run_validation_benchmark.py`
Scores the human re-extraction study (`data/validation_reactions.json`) against the
ground truth under **both** similarity methods and prints a per-metric comparison — use it
to confirm the Tanimoto option preserves strong scores on trusted data.

**Usage:**
```bash
python scripts/run_validation_benchmark.py
```

### `extract_scores.py`
Extracts key metrics from scoring results for leaderboard display.

**Usage:**
```bash
python scripts/extract_scores.py \
  --results-dir results/your_name_timestamp \
  --output results/your_name_timestamp/scores.json
```

**Output:** JSON file with mean, median, and std scores for all categories.

### `update_leaderboard.py`
Updates the README leaderboard table and leaderboard.json with new scores.

**Usage:**
```bash
python scripts/update_leaderboard.py \
  --submitter "Your Name" \
  --scores-file results/your_name_timestamp/scores.json \
  --pr-number 123
```

**What it does:**
- Adds or updates entry in `leaderboard.json`
- Sorts entries by combined score
- Updates the leaderboard table in README.md

## Manual Testing

You can test the entire pipeline locally:

```bash
# 1. Create a test submission
cp data/submissions/example_submission.json data/submissions/test_submission.json

# 2. Run scoring
python scripts/run_scoring.py \
  --submission-file data/submissions/test_submission.json \
  --output-dir results/test_$(date +%Y%m%d_%H%M%S)

# 3. Extract scores
python scripts/extract_scores.py \
  --results-dir results/test_* \
  --output results/test_*/scores.json

# 4. Update leaderboard (use a fake PR number for testing)
python scripts/update_leaderboard.py \
  --submitter "Test User" \
  --scores-file results/test_*/scores.json \
  --pr-number 0
```

## CI/CD Integration

These scripts are automatically run by the GitHub Actions workflow (`.github/workflows/score_submission.yml`) when:
- A pull request is opened or updated
- The PR modifies files in `data/submissions/`

The workflow will:
1. Extract the submitter name from the JSON file
2. Run scoring and generate plots
3. Extract metrics and update the leaderboard
4. Commit results back to the PR branch
5. Comment scores on the PR
