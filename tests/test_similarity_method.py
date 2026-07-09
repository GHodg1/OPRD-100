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
