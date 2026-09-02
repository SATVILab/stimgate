# Short-request acceptance scenarios

These prompts deliberately contain only the desired outcome. An agent using the skill must supply the procedure.

## Exact mutation

Prompt:

> Set example#313 to P2.

Expected interpretation with the single-Project fixture:

- repository `octo-org/example`, issue 313, Project 12;
- common P2 maps exactly to organisation-native `Medium`;
- inspect the issue, membership, current native Priority and preservation set;
- re-read stale-sensitive state, update only Priority, then independently read it back;
- do not ask the user to repeat those procedural requirements.

## Ranked read

Prompt:

> What are the highest-priority open items in example?

Expected interpretation with the single-Project fixture:

- resolve Project 12 and its open member issues;
- read organisation-native Priority separately from Project fields;
- translate Urgent, High, Medium and Low to P0, P1, P2 and P3;
- rank P0 before P1 before P2 before P3 and report unresolved ties;
- make no mutation.

## Command-returning surface

Prompt:

> Set example#313 to P2.

If the surface cannot write, it must return one narrow organisation issue-field mutation using `Medium` plus an independent GET readback. It must state that the commands have not run.

## Issue creation styles

Prompt:

> Add an issue to chek whether the import preserves dates, because the dry run changed one.

Expected interpretation:

- `direct` derives a concise title and a description from the supplied statement, corrects `chek`, and does not otherwise rewrite, reorganise or expand it;
- `tidy` may reword and organise the supplied facts, but must not add acceptance criteria, implementation steps or other substantive information;
- `unrestricted` may add useful detail only when it is grounded in the request or a required project source;
- a more recent explicit instruction in the prompt overrides the configured style.
