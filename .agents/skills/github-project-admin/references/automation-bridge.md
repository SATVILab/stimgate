# Issue-triggered Project administration bridge

Use this bridge only when an authorised GitHub Project mutation cannot be performed by the current chat/provider surface. It turns a bounded operation issue into a constrained Codex manifest, applies that manifest in a deterministic executor with a separate Project credential, and closes the operation issue only after independent readback.

## Security model

The operation issue is an execution queue item, not a mirror of the real task. It must contain no credentials or private material that is inappropriate for its repository.

The workflow is triggered only by the exact `automation:project-admin` label. It ignores issue comments as instructions. The operation issue body must contain the marker:

```text
<!-- project-admin-operation:v1 -->
```

Codex runs read-only and without `PROJECTS_TOKEN`. It may inspect the checked-out repository contract and the canonical skill, then emit only the v1 JSON manifest. A deterministic Node executor receives `PROJECTS_TOKEN` afterwards, re-discovers live Project and issue state, applies only the supported field/membership operations, and performs targeted readback.

Supported v1 operations are:

- ensure Project membership;
- set one single-select Project field to an exact option name;
- set one date Project field to an ISO `YYYY-MM-DD` value;
- clear one supported Project field.

The executor rejects a Project that is not declared by the checked-out `.projects` contract. It also rejects an issue repository that disagrees with the resolved contract. Cross-repository contributor items therefore remain outside this first bridge until the repository contract explicitly supports that topology.

## Operation issue format

Create the issue in the repository that owns the resolved Project contract. Use a concise title such as `Project admin: set example/repo#7 to P1` and a body like:

```markdown
<!-- project-admin-operation:v1 -->
## Requested mutation

Set `example/repo#7` Priority to `P1` in Project 40.

## Exact target

- Issue: `example/repo#7`
- Project owner: `example`
- Project owner type: `user`
- Project number: `40`

## Observed stale-sensitive state

- Issue updatedAt: `2026-09-03T08:00:00Z`
- Priority: `P2`

## Authority

This operation issue contains only the mutation explicitly authorised by the user. Preserve unrelated issue and Project state.
```

Include only stale-sensitive values that were actually inspected. Do not invent an `updatedAt` value or an expected field value merely to fill the template. The Codex manifest may omit an expected value when the originating surface could not inspect that Project field; the deterministic executor then performs the first live read immediately before the narrow write. When a field mutation requires Project membership and membership is authorised by the request or resolved contract, put `ensure_project_membership` before the field operation; the executor never adds membership implicitly.

Apply `automation:project-admin` only after the issue body is complete. The workflow removes the trigger label when it claims the operation. A blocked retry requires removing `automation:blocked` and applying `automation:project-admin` again.

## Repository setup

The repository must have Actions enabled and a label named `automation:project-admin`. Create it once:

```text
gh label create automation:project-admin --color 5319e7 --description "Run the bounded Project administration bridge"
```

Add these Actions secrets through **Settings → Secrets and variables → Actions**, or use `gh secret set` so the secret value is entered only into the secure terminal prompt:

```text
gh secret set OPENAI_API_KEY
gh secret set PROJECTS_TOKEN
```

`OPENAI_API_KEY` is used only by the pinned OpenAI Codex Action. `PROJECTS_TOKEN` is not passed to Codex. It is used only by the deterministic executor after Codex has produced a valid manifest.

For a user-owned Project, a classic personal access token with `project` and the repository access needed for the target issues is the least surprising `PROJECTS_TOKEN`. Add `read:org` when organisation discovery requires it. Prefer a shorter expiry and the minimum repository access that still covers the declared Project work. Never place either token in an issue, prompt, comment or log.

The reusable workflow automatically creates the reserved result labels `automation:running`, `automation:done` and `automation:blocked` when needed.

## Reusing the workflow from another repository

A repository can keep a tiny caller workflow and reuse the implementation from `MiguelRodo/projects`. Pin the reusable workflow to a reviewed commit SHA:

```yaml
name: Project administration bridge

on:
  issues:
    types: [labeled]

permissions:
  contents: read
  issues: write

jobs:
  bridge:
    if: github.event.label.name == 'automation:project-admin'
    uses: MiguelRodo/projects/.github/workflows/project-admin-bridge.yml@<reviewed-commit-sha>
    with:
      issue_number: ${{ github.event.issue.number }}
    secrets:
      OPENAI_API_KEY: ${{ secrets.OPENAI_API_KEY }}
      PROJECTS_TOKEN: ${{ secrets.PROJECTS_TOKEN }}
```

The reusable workflow uses GitHub's `job.workflow_repository` and `job.workflow_sha` context to check out the exact bridge implementation that defined the running job. It separately checks out the caller repository as the Project-contract target.

## Chat/provider handoff

When the current surface cannot perform an authorised Project mutation:

1. inspect the target and stale-sensitive state as far as the current surface allows;
2. create one operation issue using the format above in the repository that owns the resolved Project contract;
3. apply `automation:project-admin` only when the workflow is known to be installed and the issue is complete;
4. report the operation issue as queued, not completed;
5. treat the Project mutation as complete only after the operation issue closes with a verified executor report.

If the bridge is not installed, the trigger label is unavailable, or the workflow is blocked, fall back to the smallest executable command/readback handoff. Do not claim that creating an operation issue itself changed the Project.
