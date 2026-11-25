# AGENTS.md — ddf-pipeline

ddf-pipeline carries out several iterations of direction-dependent self-calibration on
LOFAR imaging data using DDFacet and KillMS. This file tells the agents working on this
repository who does what, how work flows, and which conventions to respect.

## Roles and routing

The **orchestrator** talks to the user, plans, routes, integrates and verifies. It is not
the best at any one discipline; it decides who does the work and checks the result.

| Role | Agent | Use it for |
| --- | --- | --- |
| Scientist | `scientist` | radio-astronomy and DDF-calibration correctness in `scripts/`, `utils/` |
| Architect | `architect` | module boundaries, data model, technical trade-offs, ADRs |
| Developer | `developer` | implementing features, fixes, refactors |
| Quality | `quality` | running the checks and acceptance criteria before anything is called done |
| Security | `security` | secrets, credentials, licence and dependency risk |
| Writer | `writer` | the documentation that ships (`AGENTS.md`, `CONTEXT.md`, `GLOSSARY.md`, mkdocs) |
| HPC | `hpc` | SLURM launchers, singularity image, Jean Zay submission, data install scripts and movement, the smoke job |
| Explore | `explore` | finding files and answering "how does this work" questions |
| General | `general` | multi-step research and execution |

Routing rules:

- Code change → `developer`, then `quality`. Always both, in that order.
- Structural decision / new subsystem / schema change → `architect` first, `developer` second.
- Numbers, formulas, "is this real science" → `scientist`.
- Credentials, personal data, network exposure → `security` before `developer` touches the design.
- Documentation that ships → `writer`.
- SLURM, singularity, cluster execution → `hpc`.
- Two agents disagree → report both positions with a recommendation; never silently pick one.

## Branches and releases

Code is split across three canonical repositories: `mhardcastle/ddf-pipeline`,
`cyriltasse/DDFacet`, `cyriltasse/killMS`. Work flows `astra-hackathon` →
`HackathonRennes_no_futures` → `HackathonRennes` (a reviewed merge); `master` is not the
integration point. CI is per-branch; the smoke test runs on a validated tag. Record
validated commits in `VERSIONS.md`. See `docs/adr/0004-branching-and-release-model.md`.

## Branch-scope convention

Work on this branch (`astra-hackathon`) is **branch-only**. Every ticket is labelled
`ci/cd`, prefixed `[astra-hackathon]`, and assigned to `MathisCertenais`. Do not touch
files outside the current deliverables without asking.

## Pipeline

Work flows through a fixed pipeline. Each stage's output is the next stage's input.

0. **Bootstrap** — `setup-matt-pocock-skills`: configure the tracker, labels, domain layout (run once).
1. **grill-with-docs** — `grilling` + `domain-modeling`: shared understanding + the domain model → `GLOSSARY.md`, `docs/adr/`.
2. **to-spec** — `to-spec`: one spec (problem / solution / user stories / implementation / testing / out of scope) → a GitHub issue.
3. **spec-review** — fresh-context subagents (`quality` + `scientist` + `architect`) review the spec with only the spec + `GLOSSARY.md`; findings fold back into the spec.
4. **to-tickets** — `to-tickets`: tracer-bullet slices with blocking edges → one GitHub issue per ticket.
5. **implement** — `implement-spec` / `implement` + `tdd`: build the tickets on `astra-hackathon`.
6. **code-review** — `code-review`: two-axis (Standards + Spec) review of the diff.

## Domain vocabulary

Use the terms in `GLOSSARY.md` exactly; do not drift to synonyms. If a term is missing,
flag it rather than inventing one.

## Decisions

Architectural decisions live in `docs/adr/`. Read the ones that touch your area; a change
that contradicts an ADR must be surfaced, not silently made.

## Agent skills

### Issue tracker

Issues and specs live as GitHub issues via the `gh` CLI. See `docs/agents/issue-tracker.md`.

### Triage labels

The five canonical roles map to these labels: `needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`; branch-scoped work also carries `ci/cd`. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context layout: root `GLOSSARY.md` + `docs/adr/`. See `docs/agents/domain.md`.
