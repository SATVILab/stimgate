#!/usr/bin/env node

import fs from 'node:fs';
import path from 'node:path';

const [manifestPath, repositoryRoot = process.cwd()] = process.argv.slice(2);
if (!manifestPath) {
  console.error('usage: apply-project-operation.mjs MANIFEST_JSON [REPOSITORY_ROOT]');
  process.exit(2);
}

const token = process.env.PROJECTS_TOKEN;
if (!token) {
  console.error('PROJECTS_TOKEN is required');
  process.exit(3);
}

const graphqlUrl = process.env.GITHUB_GRAPHQL_URL || 'https://api.github.com/graphql';

class Blocked extends Error {
  constructor(message) {
    super(message);
    this.name = 'Blocked';
  }
}

function fail(message) {
  throw new Blocked(message);
}

function readJson(file) {
  try {
    return JSON.parse(fs.readFileSync(file, 'utf8'));
  } catch (error) {
    fail(`manifest is not valid JSON: ${error.message}`);
  }
}

function assertExactKeys(object, allowed, label) {
  if (!object || typeof object !== 'object' || Array.isArray(object)) {
    fail(`${label} must be an object`);
  }
  for (const key of Object.keys(object)) {
    if (!allowed.includes(key)) fail(`${label} contains unsupported key ${key}`);
  }
}

function validateManifest(m) {
  assertExactKeys(m, ['version', 'project', 'issue', 'operations'], 'manifest');
  if (m.version !== 1) fail('manifest.version must be 1');

  assertExactKeys(m.project, ['owner', 'ownerType', 'number'], 'manifest.project');
  if (typeof m.project.owner !== 'string' || !/^[A-Za-z0-9-]+$/.test(m.project.owner)) {
    fail('manifest.project.owner is invalid');
  }
  if (!['user', 'organization'].includes(m.project.ownerType)) {
    fail('manifest.project.ownerType must be user or organization');
  }
  if (!Number.isInteger(m.project.number) || m.project.number < 1) {
    fail('manifest.project.number must be a positive integer');
  }

  assertExactKeys(m.issue, ['repository', 'number', 'expectedUpdatedAt'], 'manifest.issue');
  if (typeof m.issue.repository !== 'string' || !/^[A-Za-z0-9_.-]+\/[A-Za-z0-9_.-]+$/.test(m.issue.repository)) {
    fail('manifest.issue.repository must be owner/name');
  }
  if (!Number.isInteger(m.issue.number) || m.issue.number < 1) {
    fail('manifest.issue.number must be a positive integer');
  }
  if ('expectedUpdatedAt' in m.issue && m.issue.expectedUpdatedAt !== null && typeof m.issue.expectedUpdatedAt !== 'string') {
    fail('manifest.issue.expectedUpdatedAt must be a string or null');
  }

  if (!Array.isArray(m.operations) || m.operations.length < 1 || m.operations.length > 8) {
    fail('manifest.operations must contain between 1 and 8 operations');
  }
  for (const [index, op] of m.operations.entries()) {
    if (!op || typeof op !== 'object' || Array.isArray(op)) fail(`operation ${index + 1} must be an object`);
    if (op.kind === 'ensure_project_membership') {
      assertExactKeys(op, ['kind'], `operation ${index + 1}`);
      continue;
    }
    if (op.kind === 'set_single_select') {
      assertExactKeys(op, ['kind', 'field', 'value', 'expected'], `operation ${index + 1}`);
      if (typeof op.field !== 'string' || !op.field) fail(`operation ${index + 1} field is required`);
      if (typeof op.value !== 'string' || !op.value) fail(`operation ${index + 1} value is required`);
      if ('expected' in op && op.expected !== null && typeof op.expected !== 'string') {
        fail(`operation ${index + 1} expected must be a string or null`);
      }
      continue;
    }
    if (op.kind === 'set_date') {
      assertExactKeys(op, ['kind', 'field', 'value', 'expected'], `operation ${index + 1}`);
      if (typeof op.field !== 'string' || !op.field) fail(`operation ${index + 1} field is required`);
      if (typeof op.value !== 'string' || !/^\d{4}-\d{2}-\d{2}$/.test(op.value)) {
        fail(`operation ${index + 1} value must be YYYY-MM-DD`);
      }
      if ('expected' in op && op.expected !== null && (typeof op.expected !== 'string' || !/^\d{4}-\d{2}-\d{2}$/.test(op.expected))) {
        fail(`operation ${index + 1} expected must be YYYY-MM-DD or null`);
      }
      continue;
    }
    if (op.kind === 'clear_field') {
      assertExactKeys(op, ['kind', 'field', 'expected'], `operation ${index + 1}`);
      if (typeof op.field !== 'string' || !op.field) fail(`operation ${index + 1} field is required`);
      if ('expected' in op && op.expected !== null && typeof op.expected !== 'string') {
        fail(`operation ${index + 1} expected must be a string or null`);
      }
      continue;
    }
    fail(`operation ${index + 1} has unsupported kind ${op.kind}`);
  }
  return m;
}

function parseTable(text) {
  const values = new Map();
  for (const line of text.split(/\r?\n/)) {
    const match = line.match(/^\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|\s*$/);
    if (!match) continue;
    const key = match[1].trim();
    const value = match[2].trim();
    if (!key || key === 'Key' || /^-+$/.test(key)) continue;
    values.set(key, value);
  }
  return values;
}

function loadContracts(root) {
  const files = [];
  const rootContract = path.join(root, '.projects', 'project.md');
  if (fs.existsSync(rootContract)) files.push(rootContract);
  const leafDir = path.join(root, '.projects', 'projects');
  if (fs.existsSync(leafDir)) {
    for (const name of fs.readdirSync(leafDir).sort()) {
      if (name.endsWith('.md')) files.push(path.join(leafDir, name));
    }
  }
  return files.map(file => ({ file, values: parseTable(fs.readFileSync(file, 'utf8')) }));
}

function resolveContract(root, manifest) {
  const matches = loadContracts(root).filter(({ values }) =>
    values.get('Project owner') === manifest.project.owner &&
    Number(values.get('Project number')) === manifest.project.number
  );
  if (matches.length !== 1) {
    fail(`expected exactly one local contract for ${manifest.project.owner} Project ${manifest.project.number}; found ${matches.length}`);
  }
  const issueRepository = matches[0].values.get('Issue repository');
  if (issueRepository && issueRepository !== manifest.issue.repository) {
    fail(`manifest issue repository ${manifest.issue.repository} does not match contract issue repository ${issueRepository}`);
  }
  return matches[0];
}

async function graphql(query, variables) {
  const response = await fetch(graphqlUrl, {
    method: 'POST',
    headers: {
      authorization: `Bearer ${token}`,
      'content-type': 'application/json',
      'user-agent': 'github-project-admin-bridge/1',
    },
    body: JSON.stringify({ query, variables }),
  });
  if (!response.ok) fail(`GitHub GraphQL request failed with HTTP ${response.status}`);
  const payload = await response.json();
  if (payload.errors?.length) fail(`GitHub GraphQL error: ${payload.errors.map(e => e.message).join('; ')}`);
  return payload.data;
}

const PROJECT_QUERY = `
query($login: String!, $number: Int!) {
  user(login: $login) {
    projectV2(number: $number) {
      id
      title
      fields(first: 100) {
        nodes {
          ... on ProjectV2Field { id name dataType }
          ... on ProjectV2SingleSelectField { id name dataType options { id name } }
        }
        pageInfo { hasNextPage }
      }
    }
  }
  organization(login: $login) {
    projectV2(number: $number) {
      id
      title
      fields(first: 100) {
        nodes {
          ... on ProjectV2Field { id name dataType }
          ... on ProjectV2SingleSelectField { id name dataType options { id name } }
        }
        pageInfo { hasNextPage }
      }
    }
  }
}`;

const ISSUE_QUERY = `
query($owner: String!, $repo: String!, $number: Int!) {
  repository(owner: $owner, name: $repo) {
    issue(number: $number) {
      id
      url
      state
      updatedAt
      projectItems(first: 100, includeArchived: false) {
        nodes { id project { id number } }
        pageInfo { hasNextPage }
      }
    }
  }
}`;

const FIELD_VALUE_QUERY = `
query($item: ID!, $field: String!) {
  node(id: $item) {
    ... on ProjectV2Item {
      fieldValueByName(name: $field) {
        __typename
        ... on ProjectV2ItemFieldSingleSelectValue { name optionId }
        ... on ProjectV2ItemFieldDateValue { date }
        ... on ProjectV2ItemFieldTextValue { text }
        ... on ProjectV2ItemFieldNumberValue { number }
      }
    }
  }
}`;

async function getProject(manifest) {
  const data = await graphql(PROJECT_QUERY, { login: manifest.project.owner, number: manifest.project.number });
  const project = manifest.project.ownerType === 'user' ? data.user?.projectV2 : data.organization?.projectV2;
  if (!project) fail(`Project ${manifest.project.owner}/${manifest.project.number} was not found or is inaccessible`);
  if (project.fields.pageInfo.hasNextPage) fail('Project has more than 100 fields; bridge refuses incomplete discovery');
  return project;
}

async function getIssue(manifest) {
  const [owner, repo] = manifest.issue.repository.split('/');
  const data = await graphql(ISSUE_QUERY, { owner, repo, number: manifest.issue.number });
  const issue = data.repository?.issue;
  if (!issue) fail(`issue ${manifest.issue.repository}#${manifest.issue.number} was not found or is inaccessible`);
  if (issue.projectItems.pageInfo.hasNextPage) fail('issue belongs to more than 100 Projects; bridge refuses incomplete discovery');
  return issue;
}

function projectItem(issue, projectId) {
  const matches = issue.projectItems.nodes.filter(item => item.project?.id === projectId);
  if (matches.length > 1) fail('issue unexpectedly has duplicate items in the same Project');
  return matches[0] || null;
}

function findField(project, name) {
  const matches = project.fields.nodes.filter(field => field?.name === name);
  if (matches.length !== 1) fail(`expected exactly one Project field named ${name}; found ${matches.length}`);
  return matches[0];
}

async function getFieldValue(itemId, fieldName) {
  const data = await graphql(FIELD_VALUE_QUERY, { item: itemId, field: fieldName });
  const value = data.node?.fieldValueByName ?? null;
  if (!value) return null;
  if (value.__typename === 'ProjectV2ItemFieldSingleSelectValue') return value.name ?? null;
  if (value.__typename === 'ProjectV2ItemFieldDateValue') return value.date ?? null;
  if (value.__typename === 'ProjectV2ItemFieldTextValue') return value.text ?? null;
  if (value.__typename === 'ProjectV2ItemFieldNumberValue') return value.number == null ? null : String(value.number);
  fail(`field ${fieldName} has unsupported value type ${value.__typename}`);
}

function checkExpected(op, current) {
  if ('expected' in op && op.expected !== current) {
    fail(`stale field ${op.field}: expected ${JSON.stringify(op.expected)}, observed ${JSON.stringify(current)}`);
  }
}

async function addMembership(projectId, issueId) {
  const data = await graphql(`
    mutation($project: ID!, $content: ID!) {
      addProjectV2ItemById(input: { projectId: $project, contentId: $content }) { item { id } }
    }`, { project: projectId, content: issueId });
  return data.addProjectV2ItemById?.item?.id || fail('membership mutation returned no item id');
}

async function setSingleSelect(projectId, itemId, field, optionId) {
  await graphql(`
    mutation($project: ID!, $item: ID!, $field: ID!, $option: String!) {
      updateProjectV2ItemFieldValue(input: {
        projectId: $project, itemId: $item, fieldId: $field,
        value: { singleSelectOptionId: $option }
      }) { projectV2Item { id } }
    }`, { project: projectId, item: itemId, field: field.id, option: optionId });
}

async function setDate(projectId, itemId, field, date) {
  await graphql(`
    mutation($project: ID!, $item: ID!, $field: ID!, $date: Date!) {
      updateProjectV2ItemFieldValue(input: {
        projectId: $project, itemId: $item, fieldId: $field,
        value: { date: $date }
      }) { projectV2Item { id } }
    }`, { project: projectId, item: itemId, field: field.id, date });
}

async function clearField(projectId, itemId, field) {
  await graphql(`
    mutation($project: ID!, $item: ID!, $field: ID!) {
      clearProjectV2ItemFieldValue(input: { projectId: $project, itemId: $item, fieldId: $field }) {
        projectV2Item { id }
      }
    }`, { project: projectId, item: itemId, field: field.id });
}

async function main() {
  const manifest = validateManifest(readJson(manifestPath));
  const contract = resolveContract(repositoryRoot, manifest);
  const project = await getProject(manifest);
  let issue = await getIssue(manifest);

  if (manifest.issue.expectedUpdatedAt && issue.updatedAt !== manifest.issue.expectedUpdatedAt) {
    fail(`stale issue: expected updatedAt ${manifest.issue.expectedUpdatedAt}, observed ${issue.updatedAt}`);
  }

  let item = projectItem(issue, project.id);
  const results = [];

  for (const op of manifest.operations) {
    if (op.kind === 'ensure_project_membership') {
      if (item) {
        results.push({ kind: op.kind, status: 'already_conformant', itemId: item.id });
      } else {
        const itemId = await addMembership(project.id, issue.id);
        issue = await getIssue(manifest);
        item = projectItem(issue, project.id);
        if (!item || item.id !== itemId) fail('membership readback failed');
        results.push({ kind: op.kind, status: 'applied_verified', itemId });
      }
      continue;
    }

    if (!item) fail(`operation ${op.kind} requires Project membership; add ensure_project_membership explicitly`);
    const field = findField(project, op.field);
    const before = await getFieldValue(item.id, op.field);
    checkExpected(op, before);

    if (op.kind === 'set_single_select') {
      if (field.dataType !== 'SINGLE_SELECT' || !Array.isArray(field.options)) {
        fail(`field ${op.field} is not a single-select field`);
      }
      const options = field.options.filter(option => option.name === op.value);
      if (options.length !== 1) fail(`expected exactly one option ${op.value} in field ${op.field}; found ${options.length}`);
      if (before === op.value) {
        results.push({ kind: op.kind, field: op.field, before, after: before, status: 'already_conformant' });
      } else {
        await setSingleSelect(project.id, item.id, field, options[0].id);
        const after = await getFieldValue(item.id, op.field);
        if (after !== op.value) fail(`readback failed for ${op.field}: expected ${op.value}, observed ${after}`);
        results.push({ kind: op.kind, field: op.field, before, after, status: 'applied_verified' });
      }
      continue;
    }

    if (op.kind === 'set_date') {
      if (field.dataType !== 'DATE') fail(`field ${op.field} is not a date field`);
      if (before === op.value) {
        results.push({ kind: op.kind, field: op.field, before, after: before, status: 'already_conformant' });
      } else {
        await setDate(project.id, item.id, field, op.value);
        const after = await getFieldValue(item.id, op.field);
        if (after !== op.value) fail(`readback failed for ${op.field}: expected ${op.value}, observed ${after}`);
        results.push({ kind: op.kind, field: op.field, before, after, status: 'applied_verified' });
      }
      continue;
    }

    if (op.kind === 'clear_field') {
      if (before === null) {
        results.push({ kind: op.kind, field: op.field, before, after: null, status: 'already_conformant' });
      } else {
        await clearField(project.id, item.id, field);
        const after = await getFieldValue(item.id, op.field);
        if (after !== null) fail(`readback failed for ${op.field}: expected null, observed ${after}`);
        results.push({ kind: op.kind, field: op.field, before, after, status: 'applied_verified' });
      }
      continue;
    }
  }

  const report = {
    status: 'verified',
    contract: path.relative(repositoryRoot, contract.file),
    project: manifest.project,
    issue: { repository: manifest.issue.repository, number: manifest.issue.number, url: issue.url },
    results,
  };
  process.stdout.write(`${JSON.stringify(report, null, 2)}\n`);
}

main().catch(error => {
  if (error instanceof Blocked) {
    console.error(`BLOCKED: ${error.message}`);
    process.exit(4);
  }
  console.error(error.stack || String(error));
  process.exit(1);
});
