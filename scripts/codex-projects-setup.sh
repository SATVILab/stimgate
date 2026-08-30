#!/usr/bin/env bash
# Codex Cloud may invoke setup scripts with shell tracing enabled. Disable it
# before expanding any environment variable so credentials cannot reach logs.
set +x
set -Eeuo pipefail

readonly GH_CLI_VERSION="2.98.0"
readonly GH_CLI_SHA256_AMD64="3b8ac6b30336802fc1a858d7c084e11cdf24ac1a761ca90b68022d7d729208de"
readonly GH_CLI_SHA256_ARM64="cf689084f3a3618f7eae4a2420d335d74626d65f5e594b9828d125d69f800d86"
readonly API_VERSION="2026-03-10"
readonly PROJECT_OWNER="SATVILab"
readonly PROJECT_NUMBER="30"
readonly PROJECT_TITLE="stimgate"
readonly REPOSITORY="SATVILab/stimgate"

setup_tmp_dir=""

die() {
  echo "ERROR: $*" >&2
  exit 1
}

ensure_tmp_dir() {
  if [[ -z "$setup_tmp_dir" ]]; then
    setup_tmp_dir="$(mktemp -d)"
  fi
}

cleanup() {
  if [[ -n "$setup_tmp_dir" && "$setup_tmp_dir" == /tmp/* ]]; then
    rm -rf -- "$setup_tmp_dir"
  fi
}
trap cleanup EXIT

install_gh() {
  local architecture archive_name archive_path expected_sha extracted_binary
  local install_dir path_line shell_rc user_directory=""

  case "$(uname -m)" in
    x86_64)
      architecture="amd64"
      expected_sha="$GH_CLI_SHA256_AMD64"
      ;;
    aarch64 | arm64)
      architecture="arm64"
      expected_sha="$GH_CLI_SHA256_ARM64"
      ;;
    *)
      die "unsupported architecture for automatic gh installation: $(uname -m)"
      ;;
  esac

  command -v curl >/dev/null 2>&1 || die "curl is required to install gh"
  command -v getent >/dev/null 2>&1 || die "getent is required to locate the user tool directory"
  command -v install >/dev/null 2>&1 || die "install is required to install gh"
  command -v tar >/dev/null 2>&1 || die "tar is required to install gh"
  command -v sha256sum >/dev/null 2>&1 || die "sha256sum is required to verify gh"

  ensure_tmp_dir
  archive_name="gh_${GH_CLI_VERSION}_linux_${architecture}.tar.gz"
  archive_path="$setup_tmp_dir/$archive_name"

  curl --fail --location --retry 3 --silent --show-error --output "$archive_path" \
    "https://github.com/cli/cli/releases/download/v${GH_CLI_VERSION}/${archive_name}"
  printf '%s  %s\n' "$expected_sha" "$archive_path" | sha256sum --check --status ||
    die "downloaded gh archive failed SHA-256 verification"

  tar --no-same-owner -xzf "$archive_path" -C "$setup_tmp_dir"
  extracted_binary="$setup_tmp_dir/gh_${GH_CLI_VERSION}_linux_${architecture}/bin/gh"
  [[ -x "$extracted_binary" ]] || die "downloaded gh archive did not contain the expected binary"

  if [[ -d /usr/local/bin && -w /usr/local/bin ]]; then
    install_dir="/usr/local/bin"
  elif [[ -d /usr/bin && -w /usr/bin ]]; then
    install_dir="/usr/bin"
  else
    user_directory="$(getent passwd "$(id -u)" | cut -d: -f6)"
    [[ -n "$user_directory" && -d "$user_directory" ]] ||
      die "could not resolve a writable user directory for gh"
    install_dir="$user_directory/.local/bin"
  fi

  mkdir -p "$install_dir"
  install -m 0755 "$extracted_binary" "$install_dir/gh"

  if [[ -n "$user_directory" && "$install_dir" == "$user_directory/.local/bin" ]]; then
    path_line="$(printf 'export PATH=%q:$PATH' "$install_dir")"
    for shell_rc in "$user_directory/.bashrc" "$user_directory/.bash_profile"; do
      touch "$shell_rc"
      grep -Fqx "$path_line" "$shell_rc" || printf '%s\n' "$path_line" >>"$shell_rc"
    done
  fi
  export PATH="$install_dir:$PATH"
}

command -v gh >/dev/null 2>&1 || install_gh

command -v python >/dev/null 2>&1 || die "python is required for JSON validation"
command -v gh >/dev/null 2>&1 || die "gh installation failed"

[[ -n "${GH_TOKEN:-}" ]] ||
  die "Configure GH_TOKEN as a Codex Cloud environment variable, not a setup-only secret"

login="$(gh api user --jq .login)"
[[ -n "$login" ]] || die "GitHub authentication returned no login"

repository_name="$(gh repo view "$REPOSITORY" --json nameWithOwner --jq .nameWithOwner)"
[[ "$repository_name" == "$REPOSITORY" ]] ||
  die "expected repository $REPOSITORY but GitHub returned $repository_name"
gh issue list --repo "$REPOSITORY" --limit 1 --json number >/dev/null

ensure_tmp_dir

# Query the known organisation owner directly. This avoids the gh Project owner
# resolution path that has returned the misleading error "unknown owner type".
gh api graphql \
  -f query='query($owner: String!, $number: Int!) {
    organization(login: $owner) {
      projectV2(number: $number) {
        id
        number
        title
        repositories(first: 100) {
          nodes { nameWithOwner }
          pageInfo { hasNextPage }
        }
        fields(first: 100) {
          nodes {
            __typename
            ... on ProjectV2Field { id name dataType }
            ... on ProjectV2SingleSelectField { id name options { id name } }
            ... on ProjectV2IterationField { id name }
          }
          pageInfo { hasNextPage }
        }
      }
    }
  }' \
  -F owner="$PROJECT_OWNER" \
  -F number="$PROJECT_NUMBER" \
  >"$setup_tmp_dir/project.json"

# Issue Type and Priority are organisation-native issue metadata, not custom
# Project fields. Reading these endpoints also proves that GH_TOKEN has the
# required read:org scope.
gh api \
  -H "Accept: application/vnd.github+json" \
  -H "X-GitHub-Api-Version: $API_VERSION" \
  "/orgs/$PROJECT_OWNER/issue-types" \
  >"$setup_tmp_dir/issue-types.json"

gh api \
  -H "Accept: application/vnd.github+json" \
  -H "X-GitHub-Api-Version: $API_VERSION" \
  "/orgs/$PROJECT_OWNER/issue-fields" \
  >"$setup_tmp_dir/issue-fields.json"

python - \
  "$setup_tmp_dir/project.json" \
  "$setup_tmp_dir/issue-types.json" \
  "$setup_tmp_dir/issue-fields.json" \
  "$PROJECT_OWNER" \
  "$PROJECT_NUMBER" \
  "$PROJECT_TITLE" \
  "$REPOSITORY" <<'PY'
import json
import pathlib
import sys


def fail(message):
    raise SystemExit(f"ERROR: {message}")


(
    project_path,
    issue_types_path,
    issue_fields_path,
    project_owner,
    project_number,
    project_title,
    repository,
) = sys.argv[1:]

project_response = json.loads(pathlib.Path(project_path).read_text())
organisation = project_response.get("data", {}).get("organization")
if not organisation:
    fail(f"organisation {project_owner} is not readable")

project = organisation.get("projectV2")
if not project:
    fail(f"Project {project_owner}/{project_number} was not found or is not readable")
if project.get("number") != int(project_number):
    fail(f"expected Project number {project_number}, got {project.get('number')!r}")
if project.get("title") != project_title:
    fail(f"expected Project title {project_title!r}, got {project.get('title')!r}")

repositories = project.get("repositories", {})
if repositories.get("pageInfo", {}).get("hasNextPage"):
    fail("Project has more than 100 linked repositories; update the preflight before continuing")
repository_names = {
    node.get("nameWithOwner") for node in repositories.get("nodes", [])
}
if repository not in repository_names:
    fail(f"Project {project_number} is not linked to {repository}")

fields_connection = project.get("fields", {})
if fields_connection.get("pageInfo", {}).get("hasNextPage"):
    fail("Project has more than 100 fields; update the preflight before continuing")
fields = {
    field.get("name"): field
    for field in fields_connection.get("nodes", [])
    if field.get("name")
}

for field_name in ("Status", "Workstream"):
    field = fields.get(field_name)
    if not field:
        fail(f"Project field {field_name!r} is missing")
    if field.get("__typename") != "ProjectV2SingleSelectField":
        fail(f"Project field {field_name!r} is not a single-select field")
    if not field.get("options"):
        fail(f"Project field {field_name!r} has no options")

target_date = fields.get("Target date")
if not target_date:
    fail("Project field 'Target date' is missing")
if target_date.get("__typename") != "ProjectV2Field" or target_date.get("dataType") != "DATE":
    fail("Project field 'Target date' is not a date field")

issue_types = json.loads(pathlib.Path(issue_types_path).read_text())
issue_type_names = sorted(
    issue_type.get("name") for issue_type in issue_types if issue_type.get("name")
)
if not issue_type_names:
    fail(f"organisation {project_owner} exposes no issue types")

issue_fields = json.loads(pathlib.Path(issue_fields_path).read_text())
priority_fields = [field for field in issue_fields if field.get("name") == "Priority"]
if len(priority_fields) != 1:
    fail(f"expected exactly one organisation issue field named 'Priority', got {len(priority_fields)}")
priority = priority_fields[0]
if priority.get("data_type") != "single_select":
    fail("organisation issue field 'Priority' is not single-select")

priority_options = {
    option.get("name") for option in priority.get("options", []) if option.get("name")
}
expected_priority_options = {"Urgent", "High", "Medium", "Low"}
if priority_options != expected_priority_options:
    fail(
        "organisation Priority options changed; expected "
        f"{sorted(expected_priority_options)}, got {sorted(priority_options)}"
    )

print(f"Verified {project_owner} Project {project_number}: {project_title}.")
print(f"Verified linked repository: {repository}.")
print("Verified Project fields: Status, Workstream and Target date.")
print(f"Verified organisation issue types: {', '.join(issue_type_names)}.")
print("Verified organisation Priority options: Urgent, High, Medium and Low.")
PY

echo "Codex Cloud Project-administration preflight passed for $login."
echo "This setup intentionally does not install R or package dependencies."
