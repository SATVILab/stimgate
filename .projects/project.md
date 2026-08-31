# GitHub Project configuration

| Key | Value |
| --- | --- |
| Contract version | 1 |
| Mode | single |
| Issue repository | SATVILab/stimgate |
| Project owner | SATVILab |
| Owner type | organization |
| Project number | 30 |
| Project title | stimgate |
| Routing | linked repository; no project:* label |
| Privacy | public repository with a private organisation Project |

## Field locations

| Common dimension | Provider location | Provider field |
| --- | --- | --- |
| Class | organization issue type | Issue Type |
| Priority | organization issue field | Priority |
| Status | project field | Status |
| Workstream | project field | Workstream |
| Due date | project field | Target date |
| Parent | native issue relationship | Parent issue |
| Sub-project | repository label | subproject:* |

## Priority mapping

| Common value | Provider value |
| --- | --- |
| P0 | Urgent |
| P1 | High |
| P2 | Medium |
| P3 | Low |

The mapping is exact in both directions. In particular, Low is P3 and must never be read as P2.

## Other value mappings

- Class uses the exact live organisation Issue Type name. No aliases are declared.
- Status and Workstream use exact live Project option names. No aliases are declared.
- Due date uses an ISO 8601 date in the `Target date` Project field.

## Routing and grouping

- The Project is linked directly to this one repository, so no `project:*` routing label is required or permitted merely for membership.
- Use at most one existing `subproject:*` label when a genuine local grouping requires it.
- Do not invent a sub-project label merely to route an issue.
- Labels must not duplicate Issue Type, Priority, Status or Workstream.

## Governance

- This is a shared organisation repository. Other collaborators may edit issues and Project state.
- Assignment is explicit only. Preserve all collaborator, bot and automation changes outside the requested delta.
- Private operator sources, personal tasks and credentials are forbidden in repository content, issues, comments, commits and logs.
- The repository has no external scope-design source requirement for ordinary administration.
- R, Bioconductor, package and analysis dependencies are unrelated to Project administration.
