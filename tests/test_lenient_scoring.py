"""Regressions for global-match reporting and preserved benchmark policies."""

import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from extract_scores import extract_lenient_scores
from lenient_scoring import score_lenient


def _entry(reaction="CC>>CO", ptype="Scheme", reference="paper-a"):
    return {
        "Reference": reference,
        "Location": {"Type": ptype, "Num": "1"},
        "Reaction": reaction,
        "Steps": [{
            "Step": 1, "Yield": {}, "Reagents": [], "Solvents": [],
            "Time": None, "Temperature": None,
        }],
    }


def _score(tmp_path, gold, predictions, **kwargs):
    gold_path = tmp_path / "gold.json"
    prediction_path = tmp_path / "predictions.json"
    gold_path.write_text(json.dumps(gold))
    prediction_path.write_text(json.dumps(predictions))
    return score_lenient(prediction_path, gold_path=gold_path, **kwargs)


@pytest.mark.parametrize("type_aware", [False, True])
def test_per_type_results_and_exports_match_perfect_predictions(tmp_path, type_aware):
    entries = [
        _entry(),
        _entry("CCC>>CCO", ptype="Table"),
        _entry("CCCC>>CCCO", ptype="Experimental"),
    ]
    output = tmp_path / "results"
    score = _score(
        tmp_path, entries, entries,
        type_aware_matching=type_aware, output_dir=output,
    )
    assert score.combined_mean_covered == 1.0
    breakdown = score.per_type_breakdown()
    assert set(breakdown) == {"Scheme", "Table", "Experimental"}
    for values in breakdown.values():
        assert values["combined_mean"] == 1.0
        assert values["coverage"] == 1.0
        assert values["num_matched"] == values["num_gold"] == values["num_predicted"] == 1
    assert set(score.per_match_df["primary_type"]) == set(breakdown)
    saved = json.loads((output / "lenient_scores.json").read_text())
    assert saved["per_type_breakdown"] == breakdown
    flattened = extract_lenient_scores(output, output / "scores.json")
    assert flattened["scheme_mean"] == flattened["table_mean"] == flattened["exp_mean"] == 1


def test_global_cross_type_matches_are_attributed_to_gold_type(tmp_path):
    score = _score(
        tmp_path, [_entry(ptype="Table")], [_entry(ptype="Scheme")],
        type_aware_matching=False,
    )
    assert score.combined_mean_covered == 1.0
    breakdown = score.per_type_breakdown()
    assert set(breakdown) == {"Table"}
    assert breakdown["Table"]["num_matched"] == 1
    assert breakdown["Table"]["num_predicted"] == 1
    assert breakdown["Table"]["coverage"] == 1.0


def test_global_unmatched_predictions_keep_their_submitted_type(tmp_path):
    score = _score(
        tmp_path, [_entry(ptype="Table")],
        [_entry(), _entry("c1ccccc1>>c1ccccc1O")],
        type_aware_matching=False,
    )
    breakdown = score.per_type_breakdown()
    assert "*" not in breakdown
    assert breakdown["Table"]["num_matched"] == breakdown["Table"]["num_predicted"] == 1
    assert breakdown["Scheme"]["num_matched"] == 0
    assert breakdown["Scheme"]["num_predicted"] == 1
    assert sum(v["num_matched"] for v in breakdown.values()) == score.num_matched
    assert sum(v["num_predicted"] for v in breakdown.values()) == score.num_extracted


def test_global_rejected_matches_do_not_change_prediction_type_counts(tmp_path):
    score = _score(
        tmp_path, [_entry(ptype="Table")],
        [_entry("c1ccccc1>>c1ccccc1O")],
        type_aware_matching=False, smiles_match_threshold=1.0,
    )
    breakdown = score.per_type_breakdown()
    assert score.num_matched == 0
    assert "*" not in breakdown
    assert breakdown["Scheme"]["num_predicted"] == 1
    assert breakdown["Table"]["num_gold"] == 1
    assert breakdown["Table"]["num_predicted"] == 0


@pytest.mark.parametrize("type_aware", [False, True])
def test_empty_submission_has_zero_scores(tmp_path, type_aware):
    score = _score(
        tmp_path, [_entry()], [], type_aware_matching=type_aware
    )
    assert score.num_matched == 0
    assert score.combined_mean_covered == score.coverage == 0.0
    assert score.per_type_breakdown() == {}


def test_omitted_papers_still_do_not_count_toward_coverage(tmp_path):
    score = _score(
        tmp_path, [_entry(), _entry(reference="paper-b")], [_entry()]
    )
    assert score.num_gold == 1
    assert score.coverage == score.combined_mean_covered == 1.0


def test_extra_predictions_still_do_not_penalize_headline_score(tmp_path):
    score = _score(
        tmp_path, [_entry()], [_entry()] + [_entry("CCC>>CCO") for _ in range(9)]
    )
    assert score.num_extracted == 10
    assert score.num_matched == 1
    assert score.combined_mean_covered == 1.0


def test_invalid_tanimoto_smiles_warn_and_score_zero_in_export(tmp_path):
    with pytest.warns(RuntimeWarning, match="Invalid Tanimoto reaction SMILES"):
        score = _score(
            tmp_path, [_entry()], [_entry("CC.invalid>>CO")],
            similarity_method="tanimoto", output_dir=tmp_path / "results",
        )
    assert score.reaction_smiles_mean == 0.0
    assert score.combined_mean == pytest.approx(6 / 7)
    saved = json.loads((tmp_path / "results" / "lenient_scores.json").read_text())
    assert saved["reaction_smiles_mean"] == 0.0
