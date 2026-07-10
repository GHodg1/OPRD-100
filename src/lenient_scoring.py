"""Lenient, content-based scorer for OPRD-100 submissions.

The strict :class:`~validation.DataComparer` matches reactions by their exact
``(Reference, Location.Type, Location.Num)`` key. That contract suits *human
re-extraction*, where the annotator revisits the same locations the gold used.
Automated (LLM/OCR) extractors cannot reproduce those keys reliably — they may
label a reaction ``"Scheme"`` where the gold says ``"Figure"``, assign different
``Num`` values, or collapse/split entries — so strict location matching leaves
most reactions unscored even when the chemistry is correct.

This module provides a **content-based** alternative that preserves the strict
scorer's per-field comparison logic but changes *how reactions are paired*:

* Reactions are matched by **similarity** (via
  :class:`~validation.CompareReactionEntries`) using the Hungarian algorithm,
  **ignoring** ``Location`` labels. Matching happens within a
  ``(Reference, primary_type)`` group by default (Scheme↔Scheme, Table↔Table,
  Experimental↔Experimental) so cross-type mispairings are avoided.
* The per-field sub-scores (reaction SMILES, steps, yield, reagents, solvents,
  time, temperature) are computed by the *same* ``CompareReactionEntries`` used
  by the strict scorer, honouring the chosen ``similarity_method``
  (``"inchikey"`` or ``"tanimoto"``).
* A **coverage** term is reported: ``matched / gold`` for the papers the
  submission attempted, plus ``combined_mean * coverage`` so under-extraction is
  penalised without collapsing the per-match quality signal.

The strict scorer is untouched; this is an additional, opt-in methodology.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.optimize import linear_sum_assignment

from validation import CompareReactionEntries, remove_entry_from_location

# Repo-root default for the OPRD-100 ground-truth file. This module lives at
# src/lenient_scoring.py, so parent.parent is the repo root.
DEFAULT_GOLD_PATH = Path(__file__).resolve().parent.parent / "data" / "OPRD-100.json"

_METRICS = (
    "reaction_smiles_similarity",
    "reaction_steps_similarity",
    "yield_similarity",
    "reagent_name_similarity",
    "reagent_amount_similarity",
    "reagent_similarity",
    "solvent_similarity",
    "time_similarity",
    "temperature_similarity",
    "total_similarity",
)


# ---------------------------------------------------------------------------
# Location-type helpers

def primary_type(location_type: str) -> str:
    """Resolve the primary bucket for a ``Location.Type`` string.

    Uses the hierarchy ``Experimental > Table > Scheme > Figure``. ``Figure`` is
    kept distinct from ``Scheme`` so figure-sourced gold reactions can be
    identified separately. Handles comma-separated, ``AND``-joined and mixed
    strings, case-insensitively.

        "Scheme AND Table"       -> "Table"
        "Table AND Experimental" -> "Experimental"
        "Figure"                 -> "Figure"
        "Scheme"                 -> "Scheme"
    """
    t = (location_type or "").lower()
    if "experimental" in t:
        return "Experimental"
    if "table" in t:
        return "Table"
    if "scheme" in t:
        return "Scheme"
    if "figure" in t:
        return "Figure"
    return (location_type or "").strip()


# ---------------------------------------------------------------------------
# Stereo-stripping helpers

def _strip_smiles_stereo(smi: str) -> str:
    """Canonical SMILES with all stereo removed; original string on parse failure."""
    if not smi:
        return smi
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return smi
        Chem.RemoveStereochemistry(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=False)
    except Exception:
        return smi


def _strip_reaction_stereo(reaction_smiles: str) -> str:
    """Strip stereo from every component of a ``reactants>>products`` SMILES."""
    if not reaction_smiles or ">>" not in reaction_smiles:
        return reaction_smiles
    lhs, rhs = reaction_smiles.split(">>", 1)

    def _strip_side(side: str) -> str:
        return ".".join(_strip_smiles_stereo(s) for s in side.split("."))

    return f"{_strip_side(lhs)}>>{_strip_side(rhs)}"


def _strip_entry_stereo(entry: dict) -> dict:
    """Shallow copy of *entry* with the ``"Reaction"`` field de-stereo'd."""
    rxn = entry.get("Reaction")
    if not rxn:
        return entry
    return {**entry, "Reaction": _strip_reaction_stereo(rxn)}


# ---------------------------------------------------------------------------
# Yield normalisation helpers

def _normalize_yield_value(v: object) -> object:
    """Collapse a per-product yield list to a comparable scalar/list.

    The gold stores many yields as per-product lists — e.g. a salt product
    records ``[70.0, null]`` (70 % for the product, null for the counter-ion)
    while other rows record the same datum as a bare scalar ``70.0``. Both are
    normalised to the same form so they compare as equal:

        * all-null list          -> None
        * single non-null list   -> that value ([70.0, null] -> 70.0)
        * multi non-null list    -> list of the non-null values
        * scalar                 -> unchanged
    """
    if isinstance(v, list):
        non_null = [x for x in v if x is not None]
        if not non_null:
            return None
        if len(non_null) == 1:
            return non_null[0]
        return non_null
    return v


def _normalize_entry_yields(entry: dict) -> dict:
    """Copy of *entry* with every step's Yield dict list-flattened."""
    steps = entry.get("Steps")
    if not steps:
        return entry
    new_steps = []
    for step in steps:
        y = step.get("Yield")
        if isinstance(y, dict):
            step = {**step, "Yield": {k: _normalize_yield_value(v) for k, v in y.items()}}
        new_steps.append(step)
    return {**entry, "Steps": new_steps}


# ---------------------------------------------------------------------------
# Per-pair comparison

def _compare_pair(gold: dict, pred: dict, similarity_method: str) -> dict[str, float]:
    """Compute all per-pair metrics using the strict ``CompareReactionEntries``.

    The lenient scorer matches reactions by content (Hungarian assignment), so
    the gold/pred ``Location`` labels are irrelevant here. ``CompareReactionEntries``
    zeroes *every* metric when the two ``Location`` dicts differ (see
    ``compare_location_data``), which would defeat location-agnostic scoring and
    make results depend on whether the extractor happened to reproduce the gold's
    exact Type/Num. Neutralise it by giving the prediction copy the gold's
    Location/Reference before comparison.

    ``similarity_method`` (``"inchikey"`` | ``"tanimoto"``) is forwarded to
    ``CompareReactionEntries`` so reaction-SMILES scoring matches the strict path.
    """
    pred = {**pred, "Location": gold.get("Location"), "Reference": gold.get("Reference")}
    cmp = CompareReactionEntries(gold, pred, similarity_method=similarity_method)
    return {m: float(getattr(cmp, m, 0.0)) for m in _METRICS}


def _safe_mean(s: pd.Series) -> float:
    return float(s.mean()) if len(s) > 0 else 0.0


# ---------------------------------------------------------------------------
# Result container

@dataclass
class LenientScore:
    """Scores from content-based (location-agnostic) matching.

    ``coverage`` is the fraction of in-scope gold reactions that were matched.
    ``combined_mean_covered = combined_mean * coverage`` down-weights the
    per-match quality by how much of the paper was actually found.
    """

    combined_mean: float
    combined_mean_covered: float
    coverage: float

    reaction_smiles_mean: float
    reaction_steps_mean: float
    yield_mean: float
    reagent_name_mean: float
    reagent_amount_mean: float
    solvent_mean: float
    time_mean: float
    temperature_mean: float

    num_extracted: int
    num_gold: int
    num_matched: int

    similarity_method: str = "inchikey"
    per_match_df: pd.DataFrame = field(repr=False, compare=False, default_factory=pd.DataFrame)
    gold_type_counts: dict[str, int] = field(default_factory=dict)
    pred_type_counts: dict[str, int] = field(default_factory=dict)

    def per_type_breakdown(self) -> dict[str, dict[str, float]]:
        """Per-primary-type metric breakdown (Scheme / Table / Experimental / …)."""
        breakdown: dict[str, dict[str, float]] = {}
        df = self.per_match_df
        all_types = set(self.gold_type_counts) | set(self.pred_type_counts)
        if "primary_type" in df.columns:
            all_types |= set(df["primary_type"].unique())

        for ptype in sorted(all_types):
            sub = df[df["primary_type"] == ptype] if "primary_type" in df.columns else df.iloc[0:0]
            n_gold = self.gold_type_counts.get(ptype, 0)
            n_pred = self.pred_type_counts.get(ptype, 0)
            n_matched = len(sub)
            breakdown[ptype] = {
                "combined_mean": _safe_mean(sub["total_similarity"]) if n_matched else 0.0,
                "reaction_smiles_mean": _safe_mean(sub["reaction_smiles_similarity"]) if n_matched else 0.0,
                "yield_mean": _safe_mean(sub["yield_similarity"]) if n_matched else 0.0,
                "reagent_name_mean": _safe_mean(sub["reagent_name_similarity"]) if n_matched else 0.0,
                "solvent_mean": _safe_mean(sub["solvent_similarity"]) if n_matched else 0.0,
                "time_mean": _safe_mean(sub["time_similarity"]) if n_matched else 0.0,
                "temperature_mean": _safe_mean(sub["temperature_similarity"]) if n_matched else 0.0,
                "num_matched": float(n_matched),
                "num_gold": float(n_gold),
                "num_predicted": float(n_pred),
                "coverage": (n_matched / n_gold) if n_gold else 0.0,
            }
        return breakdown

    def to_dict(self) -> dict:
        return {
            "combined_mean": self.combined_mean,
            "combined_mean_covered": self.combined_mean_covered,
            "coverage": self.coverage,
            "reaction_smiles_mean": self.reaction_smiles_mean,
            "reaction_steps_mean": self.reaction_steps_mean,
            "yield_mean": self.yield_mean,
            "reagent_name_mean": self.reagent_name_mean,
            "reagent_amount_mean": self.reagent_amount_mean,
            "solvent_mean": self.solvent_mean,
            "time_mean": self.time_mean,
            "temperature_mean": self.temperature_mean,
            "num_extracted": self.num_extracted,
            "num_gold": self.num_gold,
            "num_matched": self.num_matched,
            "similarity_method": self.similarity_method,
            "per_type_breakdown": self.per_type_breakdown(),
        }

    def write(self, output_dir: str | Path) -> Path:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        self.per_match_df.to_csv(out / "lenient_matches.csv", index=False)
        scores_path = out / "lenient_scores.json"
        scores_path.write_text(json.dumps(self.to_dict(), indent=2))
        return scores_path


# ---------------------------------------------------------------------------
# Main entry point

def score_lenient(
    prediction_path: str | Path,
    gold_path: str | Path | None = None,
    output_dir: str | Path | None = None,
    similarity_method: str = "inchikey",
    smiles_match_threshold: float = 0.0,
    strip_stereo: bool = False,
    normalize_yields: bool = True,
    type_aware_matching: bool = True,
    exclude_types: list[str] | None = None,
) -> LenientScore:
    """Score a submission against the OPRD-100 gold using content-based matching.

    Args:
        prediction_path: JSON file with a list of OPRD-schema reaction dicts.
        gold_path: Gold JSON. Defaults to repo-root ``data/OPRD-100.json``.
        output_dir: If given, write ``lenient_scores.json`` + ``lenient_matches.csv``.
        similarity_method: ``"inchikey"`` (binary InChIKey-set Jaccard) or
            ``"tanimoto"`` (graded Morgan-Tanimoto). Forwarded to
            ``CompareReactionEntries`` — identical semantics to the strict scorer.
        smiles_match_threshold: Minimum ``reaction_smiles_similarity`` for a match
            to count (0 = accept every Hungarian pairing).
        strip_stereo: Strip stereochemistry from Reaction SMILES before comparison
            (in-memory only; the gold file is never modified).
        normalize_yields: Flatten per-product yield lists (``[70.0, null]`` == ``70.0``).
        type_aware_matching: If ``True`` (default), match only within the same
            primary Location type (Scheme↔Scheme, Table↔Table, …). If ``False``,
            match globally per paper regardless of type.
        exclude_types: Primary Location types to remove from the gold before scoring
            (e.g. ``["Figure"]``). These reactions are excluded from both the matched
            set and the coverage denominator — useful when a type was not in scope
            for the extraction (e.g. Figures were not in the human validation set).

    Returns:
        A :class:`LenientScore` with aggregate metrics, coverage, per-type
        breakdown, and the raw per-match DataFrame.
    """
    if similarity_method not in ("inchikey", "tanimoto"):
        raise ValueError(
            f"similarity_method must be 'inchikey' or 'tanimoto', got {similarity_method!r}"
        )

    gold_file = str(gold_path) if gold_path is not None else str(DEFAULT_GOLD_PATH)
    pred_data: list[dict] = json.loads(Path(prediction_path).read_text())
    gold_data: list[dict] = json.loads(Path(gold_file).read_text())

    if strip_stereo:
        pred_data = [_strip_entry_stereo(e) for e in pred_data]
        gold_data = [_strip_entry_stereo(e) for e in gold_data]

    if normalize_yields:
        pred_data = [_normalize_entry_yields(e) for e in pred_data]
        gold_data = [_normalize_entry_yields(e) for e in gold_data]

    def _ptype(entry: dict) -> str:
        loc = remove_entry_from_location(entry.get("Location", {}) or {})
        return primary_type(loc.get("Type", "") or "")

    # Drop excluded types from the gold entirely — they won't count toward
    # coverage or scores (e.g. Figure reactions were not in scope for the
    # human validation and should be excluded from AI benchmark comparisons).
    if exclude_types:
        excluded = {t.lower() for t in exclude_types}
        gold_data = [e for e in gold_data if _ptype(e).lower() not in excluded]

    def _group_key(entry: dict) -> tuple:
        ref = entry.get("Reference") or ""
        return (ref, _ptype(entry)) if type_aware_matching else (ref,)

    pred_groups: dict[tuple, list[tuple[int, dict]]] = {}
    for idx, e in enumerate(pred_data):
        pred_groups.setdefault(_group_key(e), []).append((idx, e))
    gold_groups: dict[tuple, list[tuple[int, dict]]] = {}
    for idx, e in enumerate(gold_data):
        gold_groups.setdefault(_group_key(e), []).append((idx, e))

    all_rows: list[dict] = []

    for key, pred_items in pred_groups.items():
        gold_items = gold_groups.get(key, [])
        if not gold_items:
            # Predicted reactions of a (ref, type) with no gold counterpart — unmatched.
            continue

        sim: np.ndarray = np.zeros((len(pred_items), len(gold_items)))
        pair_cache: dict[tuple[int, int], dict[str, float]] = {}
        for i, (_pi, pr) in enumerate(pred_items):
            for j, (_gj, gr) in enumerate(gold_items):
                try:
                    m = _compare_pair(gr, pr, similarity_method)
                except Exception:
                    m = {k: 0.0 for k in _METRICS}
                pair_cache[(i, j)] = m
                sim[i, j] = m["total_similarity"]

        row_ind, col_ind = linear_sum_assignment(1.0 - sim)

        for r_local, c_local in zip(row_ind, col_ind):
            pi, _pr = pred_items[r_local]
            gj, _gr = gold_items[c_local]
            m = pair_cache[(r_local, c_local)]
            if smiles_match_threshold > 0 and m["reaction_smiles_similarity"] < smiles_match_threshold:
                continue
            row = {
                "reference": key[0],
                "primary_type": key[1] if type_aware_matching else "*",
                "pred_index": pi,
                "gold_index": gj,
            }
            row.update(m)
            all_rows.append(row)

    _df_cols = ["reference", "primary_type", "pred_index", "gold_index", *list(_METRICS)]
    df = pd.DataFrame(all_rows) if all_rows else pd.DataFrame(columns=_df_cols)

    n_matched = len(df)
    n_extracted = len(pred_data)
    # Coverage denominator: all gold reactions in papers the submission attempted.
    pred_refs = {e.get("Reference") or "" for e in pred_data}
    gold_in_scope = [e for e in gold_data if (e.get("Reference") or "") in pred_refs]
    n_gold = len(gold_in_scope)
    coverage = n_matched / n_gold if n_gold > 0 else 0.0

    gold_type_counts = dict(Counter(_ptype(e) for e in gold_in_scope))
    pred_type_counts = dict(Counter(_ptype(e) for e in pred_data))

    def mean(col: str) -> float:
        return _safe_mean(df[col]) if col in df.columns else 0.0

    combined_mean = mean("total_similarity")
    score = LenientScore(
        combined_mean=combined_mean,
        combined_mean_covered=combined_mean * coverage,
        coverage=coverage,
        reaction_smiles_mean=mean("reaction_smiles_similarity"),
        reaction_steps_mean=mean("reaction_steps_similarity"),
        yield_mean=mean("yield_similarity"),
        reagent_name_mean=mean("reagent_name_similarity"),
        reagent_amount_mean=mean("reagent_amount_similarity"),
        solvent_mean=mean("solvent_similarity"),
        time_mean=mean("time_similarity"),
        temperature_mean=mean("temperature_similarity"),
        num_extracted=n_extracted,
        num_gold=n_gold,
        num_matched=n_matched,
        similarity_method=similarity_method,
        per_match_df=df,
        gold_type_counts=gold_type_counts,
        pred_type_counts=pred_type_counts,
    )

    if output_dir is not None:
        score.write(output_dir)

    return score
