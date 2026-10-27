"""Tests for the optional graded Morgan-Tanimoto reaction-SMILES similarity method."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import pytest

from validation import CompareReactionEntries, DataComparer, SIMILARITY_METHODS


def _entry(reaction: str) -> dict:
    """Minimal reaction entry with matching Location so location gating passes."""
    return {
        "Reaction": reaction,
        "Reference": "test",
        "Location": {"Type": "Scheme", "Num": "1"},
        "Steps": [{"Step": 1, "Yield": {}, "Reagents": [], "Solvents": [], "Time": None, "Temperature": None}],
    }


def test_identical_reaction_scores_one_both_methods():
    rxn = "CC(=O)O.OCC>>CCOC(C)=O"
    a = _entry(rxn)
    for method in SIMILARITY_METHODS:
        cmp = CompareReactionEntries(a, _entry(rxn), similarity_method=method)
        assert cmp.reaction_smiles_similarity == pytest.approx(1.0)


def test_ethyl_vs_butyl_ester_partial_credit_only_for_tanimoto():
    # Same transformation, different alkyl group (ethyl vs butyl).
    ethyl = _entry("CCOC(=O)C(=O)NC(C)C(=O)OCC>>CCOC(=O)c1nc(C)c(OCC)o1")
    butyl = _entry("CCCCOC(=O)C(=O)NC(C)C(=O)OCCCC>>CCCCOC(=O)c1nc(C)c(OCCCC)o1")

    inchi = CompareReactionEntries(ethyl, butyl, similarity_method="inchikey")
    tani = CompareReactionEntries(ethyl, butyl, similarity_method="tanimoto")

    # Binary InChIKey comparison sees two entirely different molecule sets → 0.
    assert inchi.reaction_smiles_similarity == pytest.approx(0.0)
    # Graded Tanimoto gives substantial partial credit (same scaffold, different chain).
    assert 0.4 < tani.reaction_smiles_similarity < 1.0


def test_completely_different_reactions_low_for_both():
    a = _entry("CC(=O)O.OCC>>CCOC(C)=O")
    b = _entry("c1ccccc1N>>c1ccccc1[N+](=O)[O-]")
    for method in SIMILARITY_METHODS:
        cmp = CompareReactionEntries(a, b, similarity_method=method)
        assert cmp.reaction_smiles_similarity < 0.5


def test_datacomparer_rejects_unknown_method(tmp_path):
    sub = tmp_path / "sub.json"
    sub.write_text("[]")
    with pytest.raises(ValueError):
        DataComparer(str(sub), similarity_method="bogus")


@pytest.mark.parametrize("method", SIMILARITY_METHODS)
@pytest.mark.parametrize("key_name,name", [("Reagents", "NaOH"), ("Solvents", "THF")])
@pytest.mark.parametrize("prediction_mass,expected", [("1 g", 1.0), ("2 g", 0.0)])
@pytest.mark.parametrize("mixed_case", [False, True])
def test_entity_amounts_use_normalized_names(
    method, key_name, name, prediction_mass, expected, mixed_case
):
    gold = _entry("CC>>CO")
    pred = _entry("CC>>CO")
    gold["Steps"][0][key_name] = [
        {key_name[:-1]: name, "Amounts": {"Mass": "1 g"}}
    ]
    pred["Steps"][0][key_name] = [
        {key_name[:-1]: f" {name.lower()} " if mixed_case else name,
         "Amounts": {"Mass": prediction_mass}}
    ]

    comparison = CompareReactionEntries(gold, pred, similarity_method=method)
    names, amounts = comparison.compare_reagents_solvents(gold, pred, key_name=key_name)
    assert names == pytest.approx(1.0)
    assert amounts == pytest.approx(expected)
    metric = "reagent_similarity" if key_name == "Reagents" else "solvent_similarity"
    assert getattr(comparison, metric) == pytest.approx(expected)


@pytest.mark.parametrize("key_name,name", [("Reagents", "NaOH"), ("Solvents", "THF")])
def test_case_variants_share_one_amount_bucket(key_name, name):
    gold = _entry("CC>>CO")
    pred = _entry("CC>>CO")
    gold["Steps"][0][key_name] = [
        {key_name[:-1]: name, "Amounts": {"Mass": "2 g"}}
    ]
    pred["Steps"][0][key_name] = [
        {key_name[:-1]: name, "Amounts": {"Mass": "1 g"}},
        {key_name[:-1]: f" {name.lower()} ", "Amounts": {"Mass": "1000 mg"}},
    ]
    comparison = CompareReactionEntries(gold, pred)
    assert comparison.compare_reagents_solvents(gold, pred, key_name=key_name) == (1.0, 1.0)


@pytest.mark.parametrize("invalid", [
    "CC.invalid>>CO",
    "CC>>CO.invalid",
    "invalid>>invalid",
    "CC.>>CO",
    "CC>>.CO",
    ">>",
])
@pytest.mark.parametrize("invalid_gold", [False, True])
def test_invalid_tanimoto_components_have_no_molecular_credit(invalid, invalid_gold):
    gold = _entry(invalid if invalid_gold else "CC>>CO")
    pred = _entry(invalid)
    with pytest.warns(RuntimeWarning, match="Invalid Tanimoto reaction SMILES"):
        comparison = CompareReactionEntries(gold, pred, similarity_method="tanimoto")
    assert comparison.reaction_smiles_similarity == 0.0
    assert comparison.total_similarity == pytest.approx(6 / 7)


def test_tanimoto_three_part_smiles_remain_unsupported():
    comparison = CompareReactionEntries(
        _entry("CC>O>CO"), _entry("CC>O>CO"), similarity_method="tanimoto"
    )
    assert comparison.reaction_smiles_similarity == 0.0


@pytest.mark.xfail(strict=True, reason="TODO: make Tanimoto fingerprints stereochemistry-aware")
def test_tanimoto_distinguishes_opposite_enantiomers():
    comparison = CompareReactionEntries(
        _entry("CC>>N[C@H](C)C(=O)O"),
        _entry("CC>>N[C@@H](C)C(=O)O"),
        similarity_method="tanimoto",
    )
    assert comparison.reaction_smiles_similarity < 1.0
