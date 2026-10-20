# GitHub Project configuration

| Key | Value |
| --- | --- |
| Contract version | 1 |
| Sub-project storage | field |
| Sub-project field | Sub-project |
| Sub-project field type | SINGLE_SELECT |
| Mode | single |
| Issue repository | SATVILab/stimgate |
| Project owner | SATVILab |
| Owner type | organization |
| Project number | 30 |
| Project title | stimgate |
| Routing | membership |
| Privacy | public repository with a private organisation Project |

## Field locations

| Common dimension | Provider location | Provider field |
| --- | --- | --- |
| Class | organization issue type | Issue Type |
| Priority | organization issue field | Priority |
| Status | project field | Status |
| Due date | project field | Target date |
| Parent | native issue relationship | Parent issue |
| Sub-project | project field | Sub-project |

## Other value mappings

- Class uses the exact live organisation Issue Type name. No aliases are declared.
- Status uses exact live Project option names. No aliases are declared.
- Due date uses an ISO 8601 date in the `Target date` Project field.

## Routing and grouping

- The Project is linked directly to this one repository, so no `project:*` routing label is required or permitted merely for membership.
- Use the Project’s `Sub-project` single-select field for genuine local grouping; leave it blank when no grouping is declared.
- Do not invent a Sub-project option merely to route an issue.
- Labels must not duplicate Issue Type, Priority or Status.

## Governance

- Collaboration mode: collaborative administration in a public repository.
- This is a shared organisation repository. Other collaborators may edit issues and Project state.
- Assignment is explicit only. Preserve all collaborator, bot and automation changes outside the requested delta.
- Private operator sources, personal tasks and credentials are forbidden in repository content, issues, comments, commits and logs.
- The repository has no external scope-design source requirement for ordinary administration.
- R, Bioconductor, package and analysis dependencies are unrelated to Project administration.

## Storage authority

- The field includes a schema-only `Not assigned` option. Leave unclassified item values blank; clearing does not select that option.
- `Sub-project` is the sole writable grouping representation. Legacy grouping labels are read-only migration artifacts until verified cleanup.
- Resolve field and option IDs freshly; require the exact configured name, SINGLE_SELECT type and vocabulary mapping. Stop if the expected field is missing or inaccessible.
- Multiple memberships require explicit Project context unless exactly one relevant managed membership resolves it. An explicit add-membership request may name its target before placement.
