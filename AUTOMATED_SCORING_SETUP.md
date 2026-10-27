# Automated scoring setup

OPRD-100 accepts public extraction submissions from both forks and branches in
this repository. GitHub-hosted runners evaluate data; no public service needs to
run on the maintainer's server.

## Submission lifecycle

1. A participant opens a PR adding or modifying **one** JSON file directly in
   `data/submissions/`. Names may contain letters, digits, underscores, and dashes.
   The PR must not change any other files.
2. `Validate PR File Changes` validates participant PRs. Code/documentation PRs
   without a participant JSON are maintenance PRs and pass this check without
   being scored. Bundled example/test JSON files are not participant submissions.
3. `Score Submission` uses `pull_request_target` to run the trusted **base
   revision** of the scorer. It downloads only the selected JSON at the event's
   exact head revision; it never checks out or executes the submitter's code.
4. The evaluation job has only read permissions. It validates the JSON, scores
   with `--lenient --similarity-method inchikey`, and uploads score summaries and
   per-match diagnostics as a 30-day artifact.
5. A separate publication job checks the artifact's PR/revision provenance and
   validates numeric scores. It posts or updates one bot comment containing
   provisional results. It does **not** commit to the participant's branch.
6. On merge, the closed-PR event re-evaluates the accepted submission and updates
   only `README.md` and `leaderboard_lenient.json` on `main`. Closing a PR without
   merging does not change the leaderboard.

The strict scoring track remains available locally for reproducing human
re-extraction validation. It is not the automated-extraction leaderboard.
The scorer's attempted-paper coverage and treatment of extra predictions are
unchanged.

## Files and permissions

| File | Purpose |
|---|---|
| `.github/workflows/validate_pr_changes.yml` | Read-only maintenance/submission classification |
| `.github/workflows/score_submission.yml` | Trusted evaluation, comments, and post-merge publication |
| `.github/scripts/submissions.cjs` | Changed-file selection, revision-pinned download, stale-result checks, comments |
| `scripts/run_scoring.py` | Validates and scores participant JSON |
| `scripts/extract_scores.py` | Produces the flat leaderboard metrics |
| `scripts/submission_workflow.py` | Allowlists and validates publication data |
| `scripts/update_leaderboard.py` | Renders leaderboard JSON and Markdown |

Only the publication job receives `contents: write` and `pull-requests: write`.
It checks out trusted `main`, never a fork. User names, paths, and JSON values are
not interpolated into shell commands. Raw submission files, metadata containing
contact details, and generated result directories are not committed by CI.
An optional repository URL is included in the public leaderboard.

This is a **public** evaluation route: submitted JSON, workflow logs, comments,
and artifacts are public. Removing contact details from the publication summary
does not make the original submission private.

## Deployment prerequisites

- Merge these workflows, helpers, and tests into `main` first. The target-event
  scoring workflow deliberately cannot execute a proposed workflow change from
  an unmerged PR.
- GitHub Actions must be enabled, and repository/organisation policy must allow
  the publication job's requested permissions.
- Required checks should distinguish `validate-changes` from the optional
  submission jobs. Maintenance PRs have no scoring/publication job to require.
- `main` rules must permit the intended bot publication. If protection blocks the
  push, the job fails explicitly: do not bypass protections. Arrange an approved
  publication route with the repository administrator before launch.
- Concurrent merged submissions may race when pushing leaderboard updates.
  A rejected non-fast-forward push fails rather than overwriting newer results;
  rerun the failed workflow after the earlier update completes.

## Disposable hosted smoke test

After deployment:

1. From a new branch, add a small, clearly labelled submission such as
   `data/submissions/workflow_smoke.json`. Use a fresh filename, not a bundled
   example.
2. Open a PR to `main`. Verify file validation, scoring, the artifact, and the
   provisional comment. Re-run evaluation to check that the bot updates its
   comment rather than adding duplicates.
3. Verify that `main`'s leaderboard and the PR branch receive no bot commits.
4. Push an updated submission and verify that the comment corresponds to the new
   head SHA, not an earlier result.
5. Close the PR **without merging**, then delete its branch. No rollback on
   `main` is needed.
6. Repeat from a fork to verify GitHub's actual cross-repository permissions.

The post-merge update should be exercised with a genuine accepted submission, or
in a separate disposable repository if testing the final push itself. Do not
merge a dummy participant just to exercise publication on the real leaderboard.

## Local checks

```bash
python -m pytest -q tests/test_submission_workflow.py tests/test_leaderboard.py \
  tests/test_lenient_scoring.py tests/test_similarity_method.py \
  tests/test_security.py tests/test_metal_catalysts.py
node --test tests/submissions.test.cjs
actionlint .github/workflows/score_submission.yml .github/workflows/validate_pr_changes.yml
```

The Python smoke test runs the real scoring and extraction CLIs, validates the
publication artifact, and updates a disposable leaderboard in a temporary
directory. JavaScript tests mock GitHub's API to cover fork downloads, maintenance
PRs, changed-file selection, stale revisions, comment updates, and the merge gate.
These checks do **not** substitute for running the deployed workflows on GitHub.

## Failure handling

- Invalid or mixed submission PRs fail with an explicit file-selection error.
- Scoring failures/timeouts stop publication; they do not produce success-shaped
  score records. Inspect the evaluation job's logs for validation/scoring errors.
- Stale or closed-unmerged PRs cannot publish old results.
- Summaries with invalid scores, counts, URLs, names, or revision identifiers are
  rejected before leaderboard changes.
- Artifacts are attached to the workflow run, not stored in Git. Download them
  before their retention period expires if a longer audit trail is needed.
