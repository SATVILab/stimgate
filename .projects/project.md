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
| Issue write-up style | tidy |
| Issue prose style | natural-direct |
| Chat implementation label | pj:implement-chat |

## Field locations

| Common dimension | Provider location | Provider field |
| --- | --- |
| Class | organization issue type | Issue Type |
| Priority | organization issue field | Priority |
| Status | project field | Status |
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

The mapping is exact in both directions. The SATVILab provider names remain unchanged because changing the organisation-wide Priority field would require organisation-admin rights and is not required for canonical P0-P3 semantics.

## Class vocabulary

Use the current shared organisation Issue Type vocabulary:

- Task
- Bug
- Enhancement
- Data
- Analysis
- Deliverable
- Documentation
- Epic

`Deliverable` replaces the retired `Report` value. Retired `Raw data` and `Processed data` issues have been migrated to `Data` and those values are not part of the active vocabulary.

## Other value mappings

- Class uses the exact live organisation Issue Type name. No aliases are declared.
- Status uses exact live Project option names. No aliases are declared.
- Due date uses an ISO 8601 date in the `Target date` Project field.

## Routing and grouping

- The Project is linked directly to this one repository, so no `project:*` routing label is required or permitted merely for membership.
- Use at most one existing `subproject:*` label when a genuine local grouping requires it.
- Do not invent a sub-project label merely to route an issue.
- Labels must not duplicate Issue Type, Priority or Status.

## Governance

- This is a shared organisation repository. Other collaborators may edit issues and Project state.
- Assignment is explicit only. Preserve all collaborator, bot and automation changes outside the requested delta.
- Private operator sources, personal tasks and credentials are forbidden in repository content, issues, comments, commits and logs.
- The repository has no external scope-design source requirement for ordinary administration.
- R, Bioconductor, package and analysis dependencies are unrelated to Project administration.
