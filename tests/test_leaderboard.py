"""Optional code-repository links in leaderboard publication."""

import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from update_leaderboard import (
    MODE_CONFIG,
    _render_lenient_table,
    _render_strict_table,
    update_leaderboard,
)


def _scores():
    return {
        "combined_mean": 0.9,
        "combined_mean_covered": 0.8,
        "coverage": 0.9,
    }


@pytest.mark.parametrize("mode", ["strict", "lenient"])
@pytest.mark.parametrize("repository_url", [
    "https://github.com/team/extractor", "", "  ", None,
])
def test_optional_repository_link_flows_from_metadata(tmp_path, monkeypatch, mode, repository_url):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "README.md").write_text("# Benchmark\n\n## Contributing\n")
    scores_path = tmp_path / "scores.json"
    scores_path.write_text(json.dumps(_scores()))
    (tmp_path / "metadata.json").write_text(json.dumps({
        "repository_url": repository_url,
        "contact_email": "private@example.org",
        "method_description": "Not part of the leaderboard record",
    }))

    update_leaderboard("Example team", scores_path, 42, mode=mode)

    saved = json.loads((tmp_path / MODE_CONFIG[mode]["leaderboard_path"]).read_text())
    entry = saved["entries"][0]
    assert entry["repository_url"] == (repository_url.strip() or None if repository_url else None)
    assert entry["scores"] == _scores()
    assert entry["rank"] == 1
    assert "private@example.org" not in json.dumps(saved)
    readme = (tmp_path / "README.md").read_text()
    assert "Code repository" in readme
    assert "[PR #42]" in readme
    if repository_url and repository_url.strip():
        assert f"[Code](<{repository_url.strip()}>)" in readme
    else:
        assert "Not shared" in readme
        assert "[Code]" not in readme


@pytest.mark.parametrize("mode", ["strict", "lenient"])
def test_legacy_runs_without_metadata_or_urls_still_render(tmp_path, monkeypatch, mode):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "README.md").write_text("# Benchmark\n")
    scores_path = tmp_path / "scores.json"
    scores_path.write_text(json.dumps(_scores()))
    update_leaderboard("Example team", scores_path, 42, mode=mode)
    saved = json.loads((tmp_path / MODE_CONFIG[mode]["leaderboard_path"]).read_text())
    del saved["entries"][0]["repository_url"]
    render = _render_lenient_table if mode == "lenient" else _render_strict_table
    assert "Not shared" in render(saved)


def test_resubmission_can_withdraw_link_without_changing_score(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "README.md").write_text("# Benchmark\n")
    scores_path = tmp_path / "scores.json"
    scores_path.write_text(json.dumps(_scores()))
    metadata_path = tmp_path / "metadata.json"
    metadata_path.write_text(json.dumps({"repository_url": "https://gitlab.com/team/code"}))
    update_leaderboard("Team", scores_path, 1, mode="lenient")
    metadata_path.write_text(json.dumps({"repository_url": ""}))
    update_leaderboard("Team", scores_path, 2, mode="lenient")
    saved = json.loads((tmp_path / "leaderboard_lenient.json").read_text())
    assert len(saved["entries"]) == 1
    assert saved["entries"][0]["repository_url"] is None
    assert saved["entries"][0]["scores"] == _scores()
    readme = (tmp_path / "README.md").read_text()
    assert "Not shared" in readme
    assert "gitlab.com" not in readme


@pytest.mark.parametrize("url", [
    "javascript:alert(1)", "file:///tmp/code", "https://user:token@github.com/team/code",
    "https://", "https://github.com/team/\ncode", 42,
])
def test_invalid_repository_url_is_rejected_before_publication(tmp_path, monkeypatch, url):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "README.md").write_text("# Benchmark\n")
    scores_path = tmp_path / "scores.json"
    scores_path.write_text(json.dumps(_scores()))
    (tmp_path / "metadata.json").write_text(json.dumps({"repository_url": url}))
    with pytest.raises(ValueError, match="repository_url"):
        update_leaderboard("Team", scores_path, 1, mode="lenient")
    assert not (tmp_path / "leaderboard_lenient.json").exists()
    assert (tmp_path / "README.md").read_text() == "# Benchmark\n"


def test_repository_url_cannot_break_markdown_table():
    leaderboard = {"entries": [{
        "rank": 1, "submitter": "Team", "date": "2026-10-27", "pr_number": 42,
        "scores": _scores(),
        "repository_url": "https://example.org/team/code?branch=(main)|<tag>",
    }]}
    table = _render_lenient_table(leaderboard)
    assert "[Code](<https://example.org/team/code?branch=%28main%29%7C%3Ctag%3E>)" in table
