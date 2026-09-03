#!/usr/bin/env node

import assert from 'node:assert/strict';
import fs from 'node:fs';
import http from 'node:http';
import os from 'node:os';
import path from 'node:path';
import { spawn } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const testDir = path.dirname(fileURLToPath(import.meta.url));
const skillDir = path.resolve(testDir, '..');
const executor = path.join(skillDir, 'scripts', 'apply-project-operation.mjs');
const fixtureRoot = path.join(testDir, 'fixtures', 'bridge-single');
const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'project-admin-bridge-'));

const state = {
  membership: false,
  itemId: 'ITEM_1',
  fields: new Map(),
  mutations: 0,
};

function reply(res, data) {
  res.writeHead(200, { 'content-type': 'application/json' });
  res.end(JSON.stringify({ data }));
}

const server = http.createServer((req, res) => {
  let body = '';
  req.setEncoding('utf8');
  req.on('data', chunk => { body += chunk; });
  req.on('end', () => {
    const { query, variables = {} } = JSON.parse(body);

    if (query.includes('query($login: String!, $number: Int!)')) {
      return reply(res, {
        user: {
          projectV2: {
            id: 'PVT_1',
            title: 'example',
            fields: {
              nodes: [
                { id: 'FIELD_CLASS', name: 'Class', dataType: 'SINGLE_SELECT', options: [
                  { id: 'OPT_TASK', name: 'Task' },
                  { id: 'OPT_ENH', name: 'Enhancement' },
                ] },
                { id: 'FIELD_STATUS', name: 'Status', dataType: 'SINGLE_SELECT', options: [
                  { id: 'OPT_TODO', name: 'Todo' },
                  { id: 'OPT_DONE', name: 'Done' },
                ] },
                { id: 'FIELD_DATE', name: 'Target date', dataType: 'DATE' },
              ],
              pageInfo: { hasNextPage: false },
            },
          },
        },
        organization: null,
      });
    }

    if (query.includes('query($owner: String!, $repo: String!, $number: Int!)')) {
      return reply(res, {
        repository: {
          issue: {
            id: 'ISSUE_1',
            url: 'https://github.com/example/repo/issues/7',
            state: 'OPEN',
            updatedAt: '2026-09-03T08:00:00Z',
            projectItems: {
              nodes: state.membership ? [{ id: state.itemId, project: { id: 'PVT_1', number: 40 } }] : [],
              pageInfo: { hasNextPage: false },
            },
          },
        },
      });
    }

    if (query.includes('query($item: ID!, $field: String!)')) {
      const value = state.fields.get(variables.field) ?? null;
      let fieldValueByName = null;
      if (value !== null) {
        if (variables.field === 'Target date') {
          fieldValueByName = { __typename: 'ProjectV2ItemFieldDateValue', date: value };
        } else {
          fieldValueByName = { __typename: 'ProjectV2ItemFieldSingleSelectValue', name: value, optionId: 'OPTION' };
        }
      }
      return reply(res, { node: { fieldValueByName } });
    }

    if (query.includes('addProjectV2ItemById')) {
      state.membership = true;
      state.mutations += 1;
      return reply(res, { addProjectV2ItemById: { item: { id: state.itemId } } });
    }

    if (query.includes('updateProjectV2ItemFieldValue') && query.includes('$option: String!')) {
      const optionById = { OPT_TASK: 'Task', OPT_ENH: 'Enhancement', OPT_TODO: 'Todo', OPT_DONE: 'Done' };
      const fieldById = { FIELD_CLASS: 'Class', FIELD_STATUS: 'Status' };
      state.fields.set(fieldById[variables.field], optionById[variables.option]);
      state.mutations += 1;
      return reply(res, { updateProjectV2ItemFieldValue: { projectV2Item: { id: state.itemId } } });
    }

    if (query.includes('updateProjectV2ItemFieldValue') && query.includes('$date: Date!')) {
      state.fields.set('Target date', variables.date);
      state.mutations += 1;
      return reply(res, { updateProjectV2ItemFieldValue: { projectV2Item: { id: state.itemId } } });
    }

    if (query.includes('clearProjectV2ItemFieldValue')) {
      const fieldById = { FIELD_CLASS: 'Class', FIELD_STATUS: 'Status', FIELD_DATE: 'Target date' };
      state.fields.delete(fieldById[variables.field]);
      state.mutations += 1;
      return reply(res, { clearProjectV2ItemFieldValue: { projectV2Item: { id: state.itemId } } });
    }

    res.writeHead(500, { 'content-type': 'application/json' });
    res.end(JSON.stringify({ errors: [{ message: `unexpected query: ${query}` }] }));
  });
});

await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
const address = server.address();
const endpoint = `http://127.0.0.1:${address.port}`;

async function run(manifest) {
  const file = path.join(tmp, `manifest-${Math.random().toString(16).slice(2)}.json`);
  fs.writeFileSync(file, `${JSON.stringify(manifest, null, 2)}\n`);
  return await new Promise((resolve, reject) => {
    const child = spawn(process.execPath, [executor, file, fixtureRoot], {
      env: { ...process.env, PROJECTS_TOKEN: 'test-token', GITHUB_GRAPHQL_URL: endpoint },
    });
    let stdout = '';
    let stderr = '';
    child.stdout.setEncoding('utf8');
    child.stderr.setEncoding('utf8');
    child.stdout.on('data', chunk => { stdout += chunk; });
    child.stderr.on('data', chunk => { stderr += chunk; });
    child.on('error', reject);
    child.on('close', code => resolve({ status: code, stdout, stderr }));
  });
}

try {
  const manifest = {
    version: 1,
    project: { owner: 'example', ownerType: 'user', number: 40 },
    issue: { repository: 'example/repo', number: 7, expectedUpdatedAt: '2026-09-03T08:00:00Z' },
    operations: [
      { kind: 'ensure_project_membership' },
      { kind: 'set_single_select', field: 'Class', value: 'Enhancement', expected: null },
      { kind: 'set_date', field: 'Target date', value: '2026-09-08', expected: null },
    ],
  };

  let result = await run(manifest);
  assert.equal(result.status, 0, result.stderr);
  const report = JSON.parse(result.stdout);
  assert.equal(report.status, 'verified');
  assert.deepEqual(report.results.map(r => r.status), ['applied_verified', 'applied_verified', 'applied_verified']);
  assert.equal(state.fields.get('Class'), 'Enhancement');
  assert.equal(state.fields.get('Target date'), '2026-09-08');
  const mutationCount = state.mutations;

  result = await run({
    ...manifest,
    operations: [
      { kind: 'ensure_project_membership' },
      { kind: 'set_single_select', field: 'Class', value: 'Enhancement', expected: 'Enhancement' },
      { kind: 'set_date', field: 'Target date', value: '2026-09-08', expected: '2026-09-08' },
    ],
  });
  assert.equal(result.status, 0, result.stderr);
  assert.equal(state.mutations, mutationCount, 'idempotent rerun unexpectedly mutated state');

  result = await run({
    ...manifest,
    operations: [{ kind: 'set_single_select', field: 'Class', value: 'Task', expected: 'Task' }],
  });
  assert.equal(result.status, 4);
  assert.match(result.stderr, /BLOCKED: stale field Class/);
  assert.equal(state.fields.get('Class'), 'Enhancement');

  result = await run({
    ...manifest,
    issue: { repository: 'other/repo', number: 7 },
    operations: [{ kind: 'ensure_project_membership' }],
  });
  assert.equal(result.status, 4);
  assert.match(result.stderr, /does not match contract issue repository/);

  console.log('project-admin bridge tests passed');
} finally {
  await new Promise(resolve => server.close(resolve));
  fs.rmSync(tmp, { recursive: true, force: true });
}
