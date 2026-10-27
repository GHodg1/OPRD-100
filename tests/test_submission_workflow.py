"""Data-only score publication and isolated CLI smoke tests."""

import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from submission_workflow import (
    COUNT_FIELDS, SCORE_FIELDS, load_publication, publish, summarize, validate_publication,
)


def _publication():
    return {
        "request": {
            "pr_number": 12, "head_sha": "a" * 40, "base_sha": "b" * 40,
            "repository": "GHodg1/OPRD-100",
        },
        "submitter": "Smoke Test",
        "repository_url": "https://github.com/example/extractor",
        "scores": {
            "mode": "lenient", "similarity_method": "inchikey",
            **{key: 1.0 for key in SCORE_FIELDS},
            **{key: 1 for key in COUNT_FIELDS},
        },
    }


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -0.1, 1.1, "1", True])
def test_invalid_scores_cannot_be_published(value):
    publication = _publication()
    publication["scores"]["combined_mean"] = value
    with pytest.raises(ValueError, match="Invalid score"):
        validate_publication(publication)


@pytest.mark.parametrize("field,value", [
    ("num_matched", -1), ("num_gold", True), ("num_extracted", 0.5),
    ("num_matched", 2),
])
def test_invalid_counts_cannot_be_published(field, value):
    publication = _publication()
    publication["scores"][field] = value
    with pytest.raises(ValueError):
        validate_publication(publication)


def test_publication_allowlists_data_and_drops_contact_details():
    publication = _publication()
    publication["contact_email"] = "private@example.org"
    publication["scores"]["unexpected"] = "private@example.org"
    cleaned = validate_publication(publication)
    assert "private@example.org" not in json.dumps(cleaned)
    assert set(cleaned["scores"]) == {
        "mode", "similarity_method", *SCORE_FIELDS, *COUNT_FIELDS,
    }


@pytest.mark.parametrize("submitter", ["", "  ", "$(echo bad)", "Team | fake row"])
def test_invalid_public_submitter_names_are_rejected(submitter):
    publication = _publication()
    publication["submitter"] = submitter
    with pytest.raises(ValueError):
        validate_publication(publication)


def test_oversized_artifact_is_rejected(tmp_path):
    (tmp_path / "publication.json").write_text(" " * (64 * 1024 + 1))
    with pytest.raises(ValueError, match="64 KB"):
        load_publication(tmp_path)


def test_publish_changes_only_leaderboard_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "README.md").write_text("# Benchmark\n")
    (tmp_path / "publication.json").write_text(json.dumps(_publication()))
    publish(tmp_path)
    leaderboard = json.loads((tmp_path / "leaderboard_lenient.json").read_text())
    entry = leaderboard["entries"][0]
    assert entry["submitter"] == "Smoke Test"
    assert entry["repository_url"] == "https://github.com/example/extractor"
    assert "Smoke Test" in (tmp_path / "README.md").read_text()


def test_real_cli_smoke_scores_and_publishes_only_in_temporary_directory(tmp_path, monkeypatch):
    source = json.loads((ROOT / "data" / "OPRD-100.json").read_text())
    reference = source[0]["Reference"]
    gold = [
        entry for entry in source
        if entry["Reference"] == reference and entry["Location"]["Type"] == "Scheme"
    ]
    assert gold
    submission = {
        "submitter_name": "Workflow smoke test",
        "repository_url": "https://github.com/example/extractor",
        "reactions": copy.deepcopy(gold),
    }
    (tmp_path / "submission.json").write_text(json.dumps(submission))
    (tmp_path / "request.json").write_text(json.dumps(_publication()["request"]))
    results = tmp_path / "results"
    commands = [
        [ROOT / "scripts/run_scoring.py", "--lenient", "--submission-file",
         tmp_path / "submission.json", "--output-dir", results],
        [ROOT / "scripts/extract_scores.py", "--mode", "lenient",
         "--results-dir", results, "--output", results / "scores.json"],
        [ROOT / "scripts/submission_workflow.py", "summarize", "--directory", tmp_path],
    ]
    for command in commands:
        subprocess.run(
            [sys.executable, *map(str, command)], check=True,
            capture_output=True, text=True, timeout=120,
        )
    publication = load_publication(tmp_path)
    assert publication["scores"]["num_matched"] == len(gold)
    assert publication["scores"]["combined_mean"] == pytest.approx(1.0)
    assert publication["scores"]["reaction_smiles_mean"] == pytest.approx(1.0)
    in_scope = [
        entry for entry in source
        if entry["Reference"] == reference
        and any(kind in entry["Location"]["Type"] for kind in ("Scheme", "Table", "Experimental"))
    ]
    assert publication["scores"]["num_gold"] == len(in_scope)
    assert publication["scores"]["coverage"] == pytest.approx(len(gold) / len(in_scope))
    assert publication["scores"]["combined_mean_covered"] == pytest.approx(len(gold) / len(in_scope))
    assert not (tmp_path / "leaderboard_lenient.json").exists()
    assert not (tmp_path / "README.md").exists()
    monkeypatch.chdir(tmp_path)
    (tmp_path / "README.md").write_text("# Disposable smoke-test leaderboard\n")
    publish(tmp_path)
    saved = json.loads((tmp_path / "leaderboard_lenient.json").read_text())
    assert saved["entries"][0]["scores"] == publication["scores"]
    assert saved["entries"][0]["repository_url"] == submission["repository_url"]
