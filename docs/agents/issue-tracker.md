# Issue tracker: GitHub

Issues and specs for this repo live as GitHub issues. Use the `gh` CLI for all operations.

## Conventions

- **Create an issue**: `gh issue create --title "..." --body "..."`. Use a heredoc for multi-line bodies.
- **Read an issue**: `gh issue view <number> --json number,title,body,labels,comments`.
- **List issues**: `gh issue list --state open --json number,title,body,labels,comments --jq '[.[] | {number, title, body, labels: [.labels[].name], comments: [.comments[].body]}]'` with appropriate `--label` and `--state` filters.
- **Blocking**: set native dependencies with `gh issue create --blocked-by <n>,<n> ...` or `gh issue edit <child> --add-blocked-by <blocker>` (`gh` 2.94+). Where dependencies aren't available, put a `Blocked by: #<n>, #<n>` line at the top of the child body. A ticket is unblocked when every blocker is closed.
- **Comment on an issue**: `gh issue comment <number> --body "..."`
- **Apply / remove labels**: `gh issue edit <number> --add-label "..."` / `--remove-label "..."`
- **Close**: post the explanation with `gh issue comment <number> --body "..."`, then `gh issue close <number>`.

Infer the repo from `git remote -v`; `gh` does this automatically when run inside a clone.

## Branch-scope convention

Work on the `astra-hackathon` branch is **branch-only** and is marked with three signals so
it cannot be mistaken for shared-codebase work:

- label **`ci/cd`**
- title prefix **`[astra-hackathon]`**
- assigned to **`MathisCertenais`**

## Pull requests as a triage surface

**PRs as a request surface: no.**
