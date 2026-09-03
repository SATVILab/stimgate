#!/usr/bin/env bash
set -Eeuo pipefail

test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
skill_dir="$(cd "$test_dir/.." && pwd)"
repo_root="$(cd "$skill_dir/../.." && pwd)"

node --check "$skill_dir/scripts/apply-project-operation.mjs"
python -m json.tool "$skill_dir/schemas/project-admin-operation-v1.schema.json" >/dev/null
node "$test_dir/project-admin-bridge.mjs"

grep -Fq "github.event.label.name == 'automation:project-admin'" \
  "$repo_root/.github/workflows/project-admin-bridge.yml"
grep -Fq 'openai/codex-action@86365089eb2b84e0a8fb0717b304f8bdcb13b20e' \
  "$repo_root/.github/workflows/project-admin-bridge.yml"
grep -Fq 'actions/checkout@d23441a48e516b6c34aea4fa41551a30e30af803' \
  "$repo_root/.github/workflows/project-admin-bridge.yml"
grep -Fq "permission-profile: ':read-only'" "$repo_root/.github/workflows/project-admin-bridge.yml"
grep -Fq 'PROJECTS_TOKEN: ${{ secrets.PROJECTS_TOKEN }}' \
  "$repo_root/.github/workflows/project-admin-bridge.yml"
grep -Fq '<!-- project-admin-operation:v1 -->' \
  "$skill_dir/references/automation-bridge.md"

echo 'project-admin bridge static tests passed'
