const fs = require('node:fs');
const path = require('node:path');

const MAX_BYTES = 10 * 1024 * 1024;
const FIXTURES = new Set(['example_submission.json', 'time_test_submission.json']);
const COMMENT_MARKER = '<!-- oprd-submission-score -->';

function selectSubmission(files) {
  const candidates = files.filter(file =>
    file.filename.startsWith('data/submissions/') &&
    file.filename.endsWith('.json') &&
    !FIXTURES.has(path.posix.basename(file.filename)));
  if (candidates.length === 0) return null;
  if (candidates.length !== 1 || files.length !== 1) {
    throw new Error('A submission PR must add or modify exactly one JSON file and no other files.');
  }
  const file = candidates[0];
  if (!/^data\/submissions\/[A-Za-z0-9][A-Za-z0-9_-]*\.json$/.test(file.filename) ||
      !['added', 'modified'].includes(file.status)) {
    throw new Error('Use an added or modified JSON directly in data/submissions/, with an alphanumeric, dash or underscore filename.');
  }
  return file;
}

function directory() {
  if (!process.env.RUNNER_TEMP) throw new Error('RUNNER_TEMP is required.');
  return path.join(process.env.RUNNER_TEMP, 'oprd-scoring');
}

async function prepare({github, context, core}) {
  core.setOutput('is_submission', 'false');
  const pr = context.payload.pull_request;
  const files = await github.paginate(github.rest.pulls.listFiles, {
    ...context.repo, pull_number: pr.number, per_page: 100,
  });
  const file = selectSubmission(files);
  if (!file) {
    core.info('Maintenance PR: no participant submission to score.');
    return;
  }
  const source = pr.merged ? pr.base.repo : pr.head.repo;
  if (!source) throw new Error('The submission repository is no longer available.');
  const response = await github.rest.repos.getContent({
    owner: source.owner.login,
    repo: source.name,
    path: file.filename,
    ref: pr.merged ? pr.merge_commit_sha : pr.head.sha,
  });
  const fileData = response.data;
  if (fileData.type !== 'file' || fileData.size > MAX_BYTES ||
      !/^[0-9a-f]{40}$/.test(fileData.sha)) {
    throw new Error('Submission must be a JSON text file no larger than 10 MB.');
  }
  const blob = await github.rest.git.getBlob({
    owner: source.owner.login, repo: source.name, file_sha: fileData.sha,
  });
  if (blob.data.encoding !== 'base64' || blob.data.size > MAX_BYTES) {
    throw new Error('Submission must be a base64-encoded Git blob no larger than 10 MB.');
  }
  const bytes = Buffer.from(blob.data.content, 'base64');
  if (bytes.length > MAX_BYTES) throw new Error('Submission exceeds 10 MB.');
  const contents = bytes.toString('utf8');
  JSON.parse(contents);
  fs.mkdirSync(directory(), {recursive: true});
  fs.writeFileSync(path.join(directory(), 'submission.json'), contents);
  fs.writeFileSync(path.join(directory(), 'request.json'), JSON.stringify({
    pr_number: pr.number,
    head_sha: pr.head.sha,
    base_sha: context.sha,
    repository: `${context.repo.owner}/${context.repo.repo}`,
  }));
  core.setOutput('is_submission', 'true');
}

function publication() {
  return JSON.parse(fs.readFileSync(path.join(directory(), 'publication.json'), 'utf8'));
}

async function currentPR({github, context}) {
  const result = publication();
  const expected = context.payload.pull_request;
  if (result.request.pr_number !== expected.number ||
      result.request.head_sha !== expected.head.sha ||
      result.request.base_sha !== context.sha ||
      result.request.repository !== `${context.repo.owner}/${context.repo.repo}`) {
    throw new Error('Score artifact does not belong to this PR revision.');
  }
  const {data: pr} = await github.rest.pulls.get({
    ...context.repo, pull_number: expected.number,
  });
  return {
    result,
    pr,
    current: pr.head.sha === result.request.head_sha &&
      (pr.state === 'open' || pr.merged),
  };
}

async function gate({github, context, core}) {
  const {pr, current} = await currentPR({github, context});
  core.setOutput('current', String(current));
  core.setOutput('publish_leaderboard', String(
    current && context.payload.action === 'closed' && pr.merged && pr.base.ref === 'main'));
  if (!current) core.info('PR changed or was closed without merging; skipping stale results.');
}

async function comment({github, context, core}) {
  const {result, pr, current} = await currentPR({github, context});
  if (!current) {
    core.info('PR changed or was closed without merging; skipping stale results.');
    return;
  }
  const s = result.scores;
  const published = context.payload.action === 'closed' && pr.merged;
  const body = `${COMMENT_MARKER}
## Scoring results

**${result.submitter}** - lenient content-based scoring, InChIKey molecular similarity.

| Score (Quality x Coverage) | Quality | Coverage | Experimental | Table | Scheme | Matched |
|---|---|---|---|---|---|---|
| ${s.combined_mean_covered.toFixed(3)} | ${s.combined_mean.toFixed(3)} | ${(s.coverage * 100).toFixed(1)}% | ${s.exp_mean.toFixed(3)} | ${s.table_mean.toFixed(3)} | ${s.scheme_mean.toFixed(3)} | ${s.num_matched}/${s.num_gold} |

${published ? 'This merged submission has been published to the leaderboard.' :
    'Provisional result: the public leaderboard changes only after this PR is merged.'}

[Score artifact and diagnostics](https://github.com/${context.repo.owner}/${context.repo.repo}/actions/runs/${context.runId})

Scorer/ground-truth revision: \`${result.request.base_sha}\`.
Submission revision: \`${result.request.head_sha}\`.
Submissions, workflow logs, and artifacts are public.`;
  const comments = await github.paginate(github.rest.issues.listComments, {
    ...context.repo, issue_number: pr.number, per_page: 100,
  });
  const previous = comments.find(item =>
    item.user?.login === 'github-actions[bot]' && item.body?.startsWith(COMMENT_MARKER));
  if (previous) {
    await github.rest.issues.updateComment({...context.repo, comment_id: previous.id, body});
  } else {
    await github.rest.issues.createComment({...context.repo, issue_number: pr.number, body});
  }
}

module.exports = {selectSubmission, prepare, gate, comment};
