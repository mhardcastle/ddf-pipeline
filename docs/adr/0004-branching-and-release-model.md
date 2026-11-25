# Branching and release model: canonical repos, a merge chain, tags

The code is split across three **canonical** GitHub repositories — `mhardcastle/ddf-pipeline`
(the pipeline), `cyriltasse/DDFacet` (imaging) and `cyriltasse/killMS` (calibration) — which
are the reference everyone clones. Work happens on `astra-hackathon` branches, merges into
`HackathonRennes_no_futures` (fast-forward), then into `HackathonRennes` (a reviewed merge).
`master` is a parallel lineage that does not contain the Hackathon work and is **not** the
integration point.

CI is branch-agnostic: every push and pull request runs the fast checks (unit tests,
compile, docs). The heavy smoke test runs only on a validated tag. Branch protection is
proposed on `HackathonRennes` after team agreement. The `HackathonRennes_no_futures ->
HackathonRennes` merge is **not** clean — it conflicts on `install_uv.sh`,
`prepare_sim_exe.sh` and `sim_exe.slurm`, which this work renamed or moved, so it is a
deliberate, tested step.

A validated run is recorded in `VERSIONS.md` (the three commits) and pinned by a tag per
repository.
