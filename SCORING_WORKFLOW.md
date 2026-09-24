# OPRD-100 Scoring Workflow

A detailed, step-by-step description of how a submission is scored against the
OPRD-100 ground truth.

There are **two scoring methodologies**:

| Mode | Matching strategy | Use for |
|------|-------------------|---------|
| **Strict** (default) | Reactions paired by exact `(Reference, Location.Type, Location.Num)` key | Human re-extraction, where the annotator revisits the same locations as the gold |
| **Lenient** (`--lenient`) | Reactions paired by chemical/field similarity (Hungarian assignment), ignoring `Location` labels | Automated extractors (LLM / OCR) whose location labels cannot match the gold verbatim |

Both modes reuse the **same per-field comparison logic**
(`CompareReactionEntries`); they differ only in *how reactions are paired* before
that logic runs.

---

## 0. Files involved

| File | Role |
|------|------|
| `scripts/run_scoring.py` | CLI entry point; orchestrates strict or lenient scoring and writes outputs |
| `src/validation.py` | `DataComparer` (strict pairing) + `CompareReactionEntries` (per-field scoring) |
| `src/lenient_scoring.py` | `score_lenient()` (content-based pairing) + `LenientScore` container |
| `src/security.py` | `validate_submission()` — security gate on the submission file |
| `src/plotting.py` | Plot helpers (strict mode only) |
| `data/OPRD-100.json` | The ground-truth reaction set |

---

## 1. Command-line entry point

```bash
# Strict (default)
python scripts/run_scoring.py \
  --submission-file data/submissions/my_submission.json \
  --output-dir results/my_run

# Lenient, graded Tanimoto reaction-SMILES similarity
python scripts/run_scoring.py \
  --submission-file data/submissions/my_submission.json \
  --output-dir results/my_run \
  --lenient --similarity-method tanimoto
```

Flags:

| Flag | Applies to | Meaning |
|------|-----------|---------|
| `--submission-file` | both | Path to the submission JSON |
| `--output-dir` | both | Where results/CSVs/plots are written |
| `--similarity-method {inchikey,tanimoto}` | both | Reaction-SMILES scoring method (default `inchikey`) |
| `--lenient` | — | Switch from strict to content-based matching |
| `--global-match` | lenient | Match globally per paper instead of within the same Location type |
| `--strip-stereo` | lenient | Remove stereochemistry from Reaction SMILES before comparison |

---

## 2. Submission format

A submission is a JSON object:

```jsonc
{
  "submitter_name": "…",
  "method_description": "…",
  "repository_url": "…",
  "paper_url": "…",
  "contact_email": "…",
  "reactions": [ /* list of OPRD-schema reaction dicts */ ]
}
```

Each reaction dict must contain at least `Reference`, `Location` (`{Type, Num}`),
`Reaction` (a `reactants>>products` SMILES or `null`), and `Steps` (a list of step
dicts, each with `Yield`, `Reagents`, `Solvents`, `Time`, `Temperature`).

---

## 3. Shared preamble (both modes)

`run_scoring.py`:

1. **Security validation** — `validate_submission()` scans the file for threats.
   On failure the run aborts with a non-zero exit code.
2. **Load submission** — reads `submitter_name`, `method_description`, and the
   `reactions` array.
3. **Write a temp copy** — the `reactions` list is written to
   `<output-dir>/submission_data.json`, which the scorer reads as the "prediction".

---

## 4. STRICT scoring pipeline

Runs when `--lenient` is **not** given. Implemented by `DataComparer` in
`src/validation.py`.

### 4.1 Build `DataComparer`

```python
dc = DataComparer(submission_data_path, similarity_method=...)
```

Inside `__init__`:

1. **Load** the gold (`data/OPRD-100.json`) and the prediction.
2. **Generate location pairs** — for each side, build the set of
   `(Reference, Type, Num)` tuples (`generate_location_pairs`). `Type`/`Num` are
   read after `remove_entry_from_location()` strips any per-entry sub-key.
3. **Record extras** — `incorrect_validation_pairs = validation_pairs − oprd_pairs`
   is the set of predicted locations absent from the gold. *(As of the current
   branch these are reported, not fatal; with `strict=True` passed to
   `DataComparer` they raise `ValueError`.)*
4. **Subset the gold** — `subset_json()` keeps only gold entries whose
   `(Reference, Type, Num)` also appears in the prediction.
5. **Split by location type** — both gold and prediction are split into
   `Experimental`, `Scheme`, and `Table` buckets (`split_by_location`), applying a
   precedence rule so an entry tagged with multiple sources is counted once
   (`Experimental` > `Table` > `Scheme`). Each bucket is sorted for determinism.

### 4.2 Match reactions within each bucket

For each location type, `compute_comparison_scores()` calls
`find_reaction_matches()`:

1. Group both sides by `(Reference, Type, Num)`.
2. For each key present on **both** sides, build a similarity matrix over the
   reactions in that group (`compare_data` → `total_similarity`).
3. Solve the **Hungarian assignment** on `cost = 1 − similarity`
   (`scipy.optimize.linear_sum_assignment`) to pair reactions optimally.
4. Emit `(pred_index, gold_index, similarity)` for each matched pair.

➡️ Only reactions whose location key exists on **both** sides can be matched.
Predicted reactions at locations not in the gold are never scored.

### 4.3 Score each matched pair

For every matched pair, `CompareReactionEntries(gold, pred, similarity_method)`
computes the per-field metrics (see §6). Results are collected into a DataFrame
per bucket, then concatenated into `combined_results`.

### 4.4 Outputs

- `experimental_results.csv`, `scheme_results.csv`, `table_results.csv`,
  `combined_results.csv`
- Score-distribution plots + entry-count comparison plots (PNG)
- `metadata.json` with submitter info and reaction counts

---

## 5. LENIENT scoring pipeline

Runs with `--lenient`. Implemented by `score_lenient()` in
`src/lenient_scoring.py`.

### 5.1 Load and normalise

1. Load gold and prediction.
2. If `--strip-stereo`: remove stereochemistry from every `Reaction` SMILES
   (in memory only — the gold file on disk is never modified).
3. **Normalise yields** — collapse per-product yield lists so `[70.0, null]`
   compares equal to `70.0` (`_normalize_entry_yields`). Gold stores some yields
   as per-product lists (e.g. a salt product records `[70.0, null]`); this removes
   that packaging noise.

### 5.2 Group reactions for matching

Each reaction is assigned a **primary type** via `primary_type()` using the
hierarchy `Experimental > Table > Scheme > Figure` (Figure kept distinct from
Scheme).

- **Type-aware (default):** group key is `(Reference, primary_type)` — a predicted
  Table row can only pair with a gold Table row, Scheme with Scheme, etc. This
  prevents spurious cross-type pairings.
- **Global (`--global-match`):** group key is `(Reference,)` — any reaction in a
  paper can pair with any other, ignoring type.

### 5.3 Hungarian matching within each group

For each `(Reference, type)` group present on both sides:

1. Build a pred × gold similarity matrix. Each cell is computed by
   `_compare_pair(gold, pred, similarity_method)`:
   - **Location is neutralised** — the prediction copy is given the gold's
     `Location`/`Reference` before comparison. This is essential: without it,
     `CompareReactionEntries` would zero every metric whenever the two `Location`
     dicts differ (see §7), defeating location-agnostic matching.
   - All per-field metrics (§6) are then computed via `CompareReactionEntries`.
2. Solve `linear_sum_assignment(1 − similarity)` to pick the optimal pairing.
3. Optionally drop pairs below `smiles_match_threshold` (default 0 = keep all).

Every accepted pairing becomes a row in the per-match DataFrame, tagged with its
`reference` and `primary_type`.

### 5.4 Aggregate + coverage

1. **Per-field means** — each metric is macro-averaged over all matched pairs.
2. **Coverage** — `num_matched / num_gold`, where `num_gold` counts gold reactions
   only in the papers the submission actually attempted (`Reference` present in the
   prediction). This rewards finding more of the paper's reactions.
3. **Coverage-weighted score** — `combined_mean_covered = combined_mean × coverage`.
4. **Per-type breakdown** — the same metrics recomputed per primary type
   (Scheme / Table / Experimental / Figure), each with its own matched/gold counts.

### 5.5 Outputs

- `lenient_scores.json` — all aggregate metrics, coverage, and the per-type breakdown
- `lenient_matches.csv` — one row per matched pair with every sub-metric
- `metadata.json` — submitter info, scoring mode, coverage

---

## 6. Per-field scoring (`CompareReactionEntries`)

Both modes score a matched pair identically. Seven sub-metrics, each in `[0, 1]`:

| Metric | How it is computed |
|--------|--------------------|
| `reaction_smiles_similarity` | See §6.1 — `inchikey` (set Jaccard) or `tanimoto` (graded) |
| `reaction_steps_similarity` | `1.0` if the two `Steps` lists have equal length, else `0.0` |
| `yield_similarity` | `dict_similarity_score` of the **final** step's `Yield` dict |
| `reagent_similarity` | `reagent_name_similarity × reagent_amount_similarity` |
| `solvent_similarity` | `solvent_name_similarity × solvent_amount_similarity` |
| `time_similarity` | Times summed to minutes across steps; `1 − |Δ| / max`. Both zero → `1.0`; one zero → `0.0` |
| `temperature_similarity` | Reduced to `[min, max]` (with `rt`→25, `reflux`→3000, etc.); `1.0` if equal, else `0.0` |

**Reagents / Solvents** (`compare_reagents_solvents`): names are standardised and
lower-cased; `name_score` is the set Jaccard of the name sets; `amount_score` is
the mean `dict_similarity_score` over the **shared** names' amount dicts.

**`dict_similarity_score`**: `1.0` for exact match; otherwise
`(# keys with matching non-null values) / (# union of non-null keys)`; `1.0` if
both dicts are entirely null.

**`total_similarity`**: the arithmetic mean of the seven metrics above.

### 6.1 Reaction-SMILES similarity methods

- **`inchikey`** (default): each side's molecules are converted to InChIKeys; the
  reactant sets and product sets are compared by Jaccard (`|∩| / |∪|`), then
  averaged. A molecule matches exactly or scores 0 — no partial credit.
- **`tanimoto`**: each molecule becomes a Morgan fingerprint (radius 2, 2048 bits);
  each side is scored by a symmetric best-match Tanimoto (`soft_set_similarity`),
  then reactant and product scores are averaged. Near-identical molecules (ethyl
  vs butyl ester, differing stereo, salt forms) receive **partial** credit.

---

## 7. The Location-matching rule (important)

`CompareReactionEntries` calls `compare_location_data()`. If the two reactions'
`Location` dicts (plus `Reference`) are **not identical**, it sets **every**
metric — including `total_similarity` — to `0`.

- In **strict** mode this is harmless: reactions are only ever paired within an
  identical `(Reference, Type, Num)` group, so the locations already match.
- In **lenient** mode this would be fatal (LLM location labels rarely match the
  gold), so `_compare_pair` **copies the gold's `Location`/`Reference` onto the
  prediction** before scoring. Matching has already been done by content, so this
  simply lets the genuine field comparisons through.

---

## 8. Choosing a mode

- Use **strict** for human re-extraction submissions that follow the OPRD-100
  annotation protocol (same reactions, same location labels).
- Use **lenient** for automated pipelines (LLM/OCR extraction) that discover
  reactions but cannot reproduce the exact location keys. Prefer
  `--similarity-method tanimoto` there, since binary InChIKey scoring gives 0 for
  chemically near-correct structures and understates real performance.

---

## 9. End-to-end example (lenient, automated extraction)

```bash
# 1. Wrap an extractor's reaction list as a submission
python - <<'PY'
import json
rxns = json.load(open("extracted.json"))          # list of OPRD reaction dicts
json.dump({"submitter_name": "MyExtractor",
           "method_description": "…",
           "reactions": rxns},
          open("submission.json", "w"), indent=2)
PY

# 2. Score it
python scripts/run_scoring.py \
  --submission-file submission.json \
  --output-dir results/myextractor \
  --lenient --similarity-method tanimoto

# 3. Inspect results
cat results/myextractor/lenient_scores.json
```
