"""Prepare and validate data-only score artifacts before privileged publication."""

import argparse
import json
import math
from pathlib import Path
import re

from update_leaderboard import _normalize_repository_url, update_leaderboard


SCORE_FIELDS = (
    "combined_mean", "combined_mean_covered", "coverage", "exp_mean", "table_mean",
    "scheme_mean", "reagent_name_mean", "reagent_amount_mean", "reaction_smiles_mean",
    "solvent_mean", "time_mean", "temperature_mean", "yield_mean", "reaction_steps_mean",
)
COUNT_FIELDS = ("num_matched", "num_gold", "num_extracted")


def validate_publication(publication: dict) -> dict:
    request = publication["request"]
    if type(request["pr_number"]) is not int or request["pr_number"] < 1:
        raise ValueError("Invalid publication PR number")
    for key in ("head_sha", "base_sha"):
        if not isinstance(request[key], str) or not re.fullmatch(r"[0-9a-f]{40}", request[key]):
            raise ValueError(f"Invalid publication {key}")
    if not isinstance(request["repository"], str) or not re.fullmatch(
        r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", request["repository"]
    ):
        raise ValueError("Invalid publication repository")
    submitter = publication["submitter"]
    if not isinstance(submitter, str) or not re.fullmatch(r"[A-Za-z0-9 _-]{1,50}", submitter):
        raise ValueError("Invalid public submitter name")
    if not submitter.strip():
        raise ValueError("Public submitter name cannot be blank")
    scores = publication["scores"]
    for field in SCORE_FIELDS:
        value = scores[field]
        if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"Invalid score: {field}")
    for field in COUNT_FIELDS:
        if type(scores[field]) is not int or scores[field] < 0:
            raise ValueError(f"Invalid count: {field}")
    if scores["num_matched"] > min(scores["num_gold"], scores["num_extracted"]):
        raise ValueError("Matched count exceeds gold or extracted count")
    return {
        "request": {key: request[key] for key in ("pr_number", "head_sha", "base_sha", "repository")},
        "submitter": submitter,
        "repository_url": _normalize_repository_url(publication.get("repository_url")),
        "scores": {
            "mode": "lenient", "similarity_method": "inchikey",
            **{key: scores[key] for key in (*SCORE_FIELDS, *COUNT_FIELDS)},
        },
    }


def summarize(directory: Path) -> dict:
    request = json.loads((directory / "request.json").read_text())
    results = directory / "results"
    metadata = json.loads((results / "metadata.json").read_text())
    scores = json.loads((results / "scores.json").read_text())
    if scores.get("similarity_method") != "inchikey" or scores.get("mode") != "lenient":
        raise ValueError("Official submissions require lenient InChIKey scoring")
    submitter = re.sub(r"[^A-Za-z0-9 _-]", "", metadata["submitter_name"])[:50].strip()
    publication = validate_publication({
        "request": request,
        "submitter": submitter,
        "repository_url": metadata.get("repository_url"),
        "scores": scores,
    })
    (directory / "publication.json").write_text(json.dumps(publication, indent=2))
    return publication


def load_publication(directory: Path) -> dict:
    source = directory / "publication.json"
    if source.stat().st_size > 64 * 1024:
        raise ValueError("Publication summary exceeds 64 KB")
    return validate_publication(json.loads(source.read_text()))


def publish(directory: Path) -> None:
    publication = load_publication(directory)
    scores_path = directory / "scores.json"
    scores_path.write_text(json.dumps(publication["scores"]))
    (directory / "metadata.json").write_text(json.dumps({
        "repository_url": publication["repository_url"],
    }))
    update_leaderboard(
        publication["submitter"], scores_path,
        publication["request"]["pr_number"], mode="lenient",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("summarize", "validate", "publish"))
    parser.add_argument("--directory", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "summarize":
        summarize(args.directory)
    elif args.command == "publish":
        publish(args.directory)
    else:
        load_publication(args.directory)
