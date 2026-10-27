"""
Script to update the README.md leaderboard(s) and leaderboard JSON with new
submission scores.

Two independent leaderboards are maintained:

* **strict** (``leaderboard.json`` / "🏆 Leaderboard") — the location-keyed
  validation scoring used for human re-extraction. Ranked by ``combined_mean``.
* **lenient** (``leaderboard_lenient.json`` / "🤖 Lenient Leaderboard") — the
  content-based scoring used for automated (AI / OCSR) extraction, where location
  labels cannot be reproduced verbatim. Ranked by ``combined_mean_covered``
  (per-match quality x coverage).
"""
import argparse
import json
import re
from datetime import datetime
from pathlib import Path
from urllib.parse import quote, urlsplit

# Per-mode configuration: which JSON file, which README section, how to rank.
MODE_CONFIG = {
    "strict": {
        "leaderboard_path": "leaderboard.json",
        "sort_key": lambda s: s.get("combined_mean", 0.0),
        "section_pattern": r"## 🔬 Validation.*?(?=\n## |\Z)",
    },
    "lenient": {
        "leaderboard_path": "leaderboard_lenient.json",
        "sort_key": lambda s: s.get("combined_mean_covered", 0.0),
        "section_pattern": r"## 🏆 Leaderboard.*?(?=\n## |\Z)",
    },
}


def _normalize_repository_url(value):
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError("repository_url must be an HTTP(S) URL or left blank")
    url = value.strip()
    if not url:
        return None
    parsed = urlsplit(url)
    if (
        parsed.scheme not in ("http", "https")
        or not parsed.hostname
        or parsed.username is not None
        or parsed.password is not None
        or any(char.isspace() or ord(char) < 32 or ord(char) == 127 for char in url)
    ):
        raise ValueError("repository_url must be an HTTP(S) URL without credentials or whitespace")
    return url


def _repository_link(entry):
    url = _normalize_repository_url(entry.get("repository_url"))
    if url is None:
        return "Not shared"
    return f"[Code](<{quote(url, safe=':/?#[]@!$&*+,;=%~')}>)"


def update_leaderboard(submitter, scores_file, pr_number, mode="strict"):
    """Update the README leaderboard and the leaderboard JSON for the given mode."""
    cfg = MODE_CONFIG[mode]
    print(f"Updating {mode} leaderboard for {submitter}")

    with open(scores_file, "r") as f:
        scores = json.load(f)

    metadata_path = Path(scores_file).with_name("metadata.json")
    try:
        with metadata_path.open() as f:
            metadata = json.load(f)
    except FileNotFoundError:
        metadata = {}
    repository_url = _normalize_repository_url(metadata.get("repository_url"))

    leaderboard_path = cfg["leaderboard_path"]
    try:
        with open(leaderboard_path, "r") as f:
            leaderboard = json.load(f)
    except FileNotFoundError:
        leaderboard = {"entries": []}

    # Update or add this submitter's entry.
    existing_entry = None
    for i, entry in enumerate(leaderboard["entries"]):
        if entry["submitter"] == submitter:
            existing_entry = i
            break

    entry = {
        "submitter": submitter,
        "date": datetime.now().strftime("%Y-%m-%d"),
        "pr_number": pr_number,
        "repository_url": repository_url,
        "scores": scores,
    }

    if existing_entry is not None:
        leaderboard["entries"][existing_entry] = entry
        print(f"Updated existing entry for {submitter}")
    else:
        leaderboard["entries"].append(entry)
        print(f"Added new entry for {submitter}")

    # Rank (descending by the mode's ranking metric).
    leaderboard["entries"].sort(key=lambda x: cfg["sort_key"](x["scores"]), reverse=True)
    for i, e in enumerate(leaderboard["entries"], 1):
        e["rank"] = i

    with open(leaderboard_path, "w") as f:
        json.dump(leaderboard, f, indent=2)
    print(f"Leaderboard saved to {leaderboard_path}")

    if mode == "lenient":
        table = _render_lenient_table(leaderboard)
    else:
        table = _render_strict_table(leaderboard)
    _write_readme_section(table, cfg["section_pattern"])
    print("README.md updated")


def _render_strict_table(leaderboard):
    table = "## 🔬 Validation — Human Re-extraction (Strict Scoring)\n\n"
    table += (
        "This table is **not** a competitive leaderboard — it documents the strict "
        "validation that establishes OPRD-100 as a faithful ground-truth source. A human "
        "independently re-extracted reactions, scored with the **strict** method (reactions "
        "matched by exact `(Reference, Type, Num)` location key, stereochemistry required, "
        "binary InChIKey SMILES matching). High scores here show OPRD-100 is internally "
        "consistent and reproducible.\n\n"
    )
    table += "| Rank | Submitter | Combined | Experimental | Table | Scheme | Reactions | Date | Code repository | Details |\n"
    table += "|------|-----------|----------|--------------|-------|--------|-----------|------|-----------------|----------|\n"
    for entry in leaderboard["entries"]:
        s = entry["scores"]
        exp = f"{s['exp_mean']:.3f}" if s.get("num_experimental", 0) > 0 else "N/A"
        table_score = f"{s['table_mean']:.3f}" if s.get("num_table", 0) > 0 else "N/A"
        scheme = f"{s['scheme_mean']:.3f}" if s.get("num_scheme", 0) > 0 else "N/A"
        table += (
            f"| {entry['rank']} | {entry['submitter']} | {s['combined_mean']:.3f} | "
            f"{exp} | {table_score} | {scheme} | {s.get('total_reactions', 0)} | "
            f"{entry['date']} | {_repository_link(entry)} | "
            f"[PR #{entry['pr_number']}](../../pull/{entry['pr_number']}) |\n"
        )
    table += "\n**Metrics explanation:**\n"
    table += "- **Combined**: Average similarity across all reaction types (higher is better, max 1.0)\n"
    table += "- **Experimental/Table/Scheme**: Average similarity for each data location type\n"
    table += "- **Reactions**: Total number of reactions scored\n\n"
    return table


def _render_lenient_table(leaderboard):
    table = "## 🏆 Leaderboard — Automated Extraction (AI / OCSR)\n\n"
    table += (
        "This is the primary OPRD-100 leaderboard. Automated extraction pipelines (LLM, "
        "OCSR, rule-based) are scored with the **lenient, content-based** method: each "
        "predicted reaction is paired to the ground truth by chemical and field similarity "
        "(Hungarian assignment) rather than by exact location labels, because automated "
        "tools rarely reproduce a paper's `(Scheme/Table/Experimental, number)` bookkeeping. "
        "The reaction-SMILES metric uses exact InChIKey-set overlap (the same metric as the "
        "strict validation scorer), and stereochemistry is kept.\n\n"
        "Entries are ranked by **Score = Quality × Coverage**, rewarding both accurate "
        "chemistry and completeness. The `Human re-extraction (reference)` row is the "
        "**ceiling** — the quality a careful human achieves under this method, and the target "
        "AI systems should aim for.\n\n"
        "Code sharing is optional. Include `repository_url` in your submission to "
        "publish a code link, or omit it / leave it blank to keep your code private. "
        "Any URL included in a public submission is public. Not sharing code does "
        "not affect scoring or ranking.\n\n"
    )
    table += "| Rank | Submitter | Score | Quality | Coverage | Experimental | Table | Scheme | Matched | Date | Code repository | Details |\n"
    table += "|------|-----------|-------|---------|----------|--------------|-------|--------|---------|------|-----------------|----------|\n"
    for entry in leaderboard["entries"]:
        s = entry["scores"]
        cov = s.get("coverage", 0.0)
        matched = f"{s.get('num_matched', 0)}/{s.get('num_gold', 0)}"
        exp = f"{s.get('exp_mean', 0.0):.3f}"
        table_score = f"{s.get('table_mean', 0.0):.3f}"
        scheme = f"{s.get('scheme_mean', 0.0):.3f}"
        table += (
            f"| {entry['rank']} | {entry['submitter']} | "
            f"{s.get('combined_mean_covered', 0.0):.3f} | {s.get('combined_mean', 0.0):.3f} | "
            f"{cov * 100:.1f}% | {exp} | {table_score} | {scheme} | {matched} | "
            f"{entry['date']} | {_repository_link(entry)} | "
            f"[PR #{entry['pr_number']}](../../pull/{entry['pr_number']}) |\n"
        )
    table += "\n**How to read the scores:**\n"
    table += "- **Score**: headline ranking metric, Quality × Coverage (max 1.0)\n"
    table += "- **Quality**: mean per-reaction similarity over matched reactions (the human-achievable ceiling)\n"
    table += "- **Coverage**: fraction of in-scope ground-truth reactions that were matched\n"
    table += "- **Experimental/Table/Scheme**: Quality broken down by ground-truth source location\n"
    table += "- **Matched**: matched reactions / ground-truth reactions in the attempted papers\n\n"
    return table


def _write_readme_section(table, section_pattern):
    readme_path = "README.md"
    with open(readme_path, "r") as f:
        readme = f.read()

    if re.search(section_pattern, readme, re.DOTALL):
        readme = re.sub(section_pattern, table.rstrip(), readme, flags=re.DOTALL)
    else:
        # Insert before "Contributing" if present, otherwise append.
        contributing_pattern = r"\n## Contributing"
        if re.search(contributing_pattern, readme):
            readme = re.sub(contributing_pattern, f"\n{table}\n## Contributing", readme)
        else:
            readme += f"\n\n{table}"

    with open(readme_path, "w") as f:
        f.write(readme)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update leaderboard with new scores")
    parser.add_argument("--submitter", required=True, help="Submitter name")
    parser.add_argument("--scores-file", required=True, help="Path to scores JSON file")
    parser.add_argument("--pr-number", required=True, type=int, help="Pull request number")
    parser.add_argument(
        "--mode", default="strict", choices=["strict", "lenient"],
        help="Which leaderboard to update: 'strict' (leaderboard.json) or "
             "'lenient' (leaderboard_lenient.json).",
    )
    args = parser.parse_args()

    update_leaderboard(args.submitter, args.scores_file, args.pr_number, mode=args.mode)
