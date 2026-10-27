const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const {test} = require('node:test');
const workflow = require('../.github/scripts/submissions.cjs');

const changed = (filename = 'data/submissions/team.json', status = 'added') => ({filename, status});
const pr = () => ({
  number: 12, state: 'open', merged: false, merge_commit_sha: 'd'.repeat(40),
  head: {sha: 'a'.repeat(40), repo: {owner: {login: 'participant'}, name: 'fork'}},
  base: {sha: 'b'.repeat(40), ref: 'main', repo: {owner: {login: 'owner'}, name: 'benchmark'}},
});

function fixture(t) {
  const previous = process.env.RUNNER_TEMP;
  process.env.RUNNER_TEMP = fs.mkdtempSync(path.join(os.tmpdir(), 'oprd-workflow-'));
  const temp = process.env.RUNNER_TEMP;
  t.after(() => {
    fs.rmSync(temp, {recursive: true});
    if (previous === undefined) delete process.env.RUNNER_TEMP;
    else process.env.RUNNER_TEMP = previous;
  });
  const context = {
    repo: {owner: 'owner', repo: 'benchmark'}, runId: 99, sha: 'b'.repeat(40),
    payload: {action: 'opened', pull_request: pr()},
  };
  const outputs = {};
  const calls = [];
  const github = {
    paginate: async method => {
      if (method === github.rest.pulls.listFiles) return [changed()];
      return [];
    },
    rest: {
      pulls: {listFiles() {}, get: async () => ({data: context.payload.pull_request})},
      repos: {getContent: async params => {
        calls.push(params);
        return {data: {type: 'file', sha: 'c'.repeat(40), size: 2}};
      }},
      git: {getBlob: async () => ({data: {
        encoding: 'base64', content: Buffer.from('{}').toString('base64'), size: 2,
      }})},
      issues: {
        listComments() {},
        createComment: async params => calls.push(params),
        updateComment: async params => calls.push(params),
      },
    },
  };
  const core = {setOutput: (key, value) => { outputs[key] = value; }, info() {}};
  return {context, github, core, outputs, calls, temp};
}

function artifact(f) {
  const directory = path.join(f.temp, 'oprd-scoring');
  fs.mkdirSync(directory, {recursive: true});
  fs.writeFileSync(path.join(directory, 'publication.json'), JSON.stringify({
    request: {
      pr_number: 12, head_sha: 'a'.repeat(40), base_sha: f.context.sha,
      repository: 'owner/benchmark',
    },
    submitter: 'Team',
    scores: {
      combined_mean: 1, combined_mean_covered: 1, coverage: 1,
      exp_mean: 1, table_mean: 1, scheme_mean: 1, num_matched: 1, num_gold: 1,
    },
  }));
}

test('maintenance and bundled fixtures are not participant submissions', () => {
  for (const files of [
    [changed('src/validation.py')], [changed('data/submissions/README.md')],
    [changed('data/submissions/example_submission.json')], [],
  ]) assert.equal(workflow.selectSubmission(files), null);
});

test('select only the changed participant JSON, never an old fixture', () => {
  assert.equal(workflow.selectSubmission([changed()]).filename, 'data/submissions/team.json');
});

test('reject mixed, multiple, renamed, deleted and unsafe submissions', () => {
  for (const files of [
    [changed(), changed('README.md')],
    [changed(), changed('data/submissions/other.json')],
    [changed('data/submissions/team.json', 'removed')],
    [changed('data/submissions/team.json', 'renamed')],
    [changed("data/submissions/$(echo bad).json")],
    [changed('data/submissions/nested/team.json')],
    [changed('data/submissions/../team.json')],
  ]) assert.throws(() => workflow.selectSubmission(files));
});

test('fork preparation fetches only JSON at the event head SHA', async t => {
  const f = fixture(t);
  await workflow.prepare(f);
  assert.equal(f.outputs.is_submission, 'true');
  assert.equal(f.calls[0].owner, 'participant');
  assert.equal(f.calls[0].repo, 'fork');
  assert.equal(f.calls[0].ref, 'a'.repeat(40));
  assert.equal(fs.readFileSync(path.join(f.temp, 'oprd-scoring/submission.json'), 'utf8'), '{}');
});

test('merged submissions are fetched from the base merge commit', async t => {
  const f = fixture(t);
  f.context.payload.pull_request.merged = true;
  f.context.sha = 'd'.repeat(40);
  await workflow.prepare(f);
  assert.equal(f.calls[0].owner, 'owner');
  assert.equal(f.calls[0].ref, 'd'.repeat(40));
  const request = JSON.parse(fs.readFileSync(path.join(f.temp, 'oprd-scoring/request.json'), 'utf8'));
  assert.equal(request.base_sha, f.context.sha);
});

test('maintenance PR preparation does not download data', async t => {
  const f = fixture(t);
  f.github.paginate = async () => [changed('README.md')];
  await workflow.prepare(f);
  assert.equal(f.outputs.is_submission, 'false');
  assert.equal(f.calls.length, 0);
});

test('reject oversized file metadata before fetching the blob', async t => {
  const f = fixture(t);
  f.github.rest.repos.getContent = async () => ({data: {
    type: 'file', size: 11 * 1024 * 1024, sha: 'c'.repeat(40),
  }});
  await assert.rejects(() => workflow.prepare(f), /10 MB/);
});

test('unmerged smoke-test PR may comment but cannot publish leaderboard', async t => {
  const f = fixture(t);
  artifact(f);
  await workflow.gate(f);
  assert.equal(f.outputs.current, 'true');
  assert.equal(f.outputs.publish_leaderboard, 'false');
  await workflow.comment(f);
  assert.match(f.calls[0].body, /Provisional result/);
});

test('only the merged closed event publishes to main', async t => {
  const f = fixture(t);
  artifact(f);
  f.context.payload.action = 'closed';
  f.context.payload.pull_request.merged = true;
  f.context.payload.pull_request.state = 'closed';
  await workflow.gate(f);
  assert.equal(f.outputs.publish_leaderboard, 'true');
});

test('stale head revisions cannot comment or publish', async t => {
  const f = fixture(t);
  artifact(f);
  f.github.rest.pulls.get = async () => ({data: {...pr(), head: {sha: 'e'.repeat(40)}}});
  await workflow.gate(f);
  await workflow.comment(f);
  assert.equal(f.outputs.current, 'false');
  assert.equal(f.outputs.publish_leaderboard, 'false');
  assert.equal(f.calls.length, 0);
});

test('closed unmerged PRs cannot comment or publish', async t => {
  const f = fixture(t);
  artifact(f);
  f.context.payload.pull_request.state = 'closed';
  await workflow.gate(f);
  await workflow.comment(f);
  assert.equal(f.outputs.current, 'false');
  assert.equal(f.calls.length, 0);
});

test('artifact provenance must match the event', async t => {
  const f = fixture(t);
  artifact(f);
  f.context.payload.pull_request.number = 13;
  await assert.rejects(() => workflow.gate(f), /does not belong/);
});

test('repeat scoring updates the existing bot comment', async t => {
  const f = fixture(t);
  artifact(f);
  f.github.paginate = async () => [{
    id: 123, user: {login: 'github-actions[bot]'}, body: '<!-- oprd-submission-score --> old',
  }];
  await workflow.comment(f);
  assert.equal(f.calls[0].comment_id, 123);
});
