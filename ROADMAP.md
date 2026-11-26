# Roadmap

Test backlog after the `astra-hackathon` debrief. Tiers and oracles are in
`docs/adr/0005-test-strategy-across-repositories.md`; metrics in
`docs/quality-control.md`.

## Done

- [x] `TestAtomicSave` (DDFacet) — concurrent write/read race on `NodesCat*.npy`.
- [x] `test_ms_quality` (ddf-pipeline) — pure-numpy metric functions.
- [x] `test_mslist_distribution` (ddf-pipeline) — nodes > MS fallback (`de84e2b`).
- [x] `test_repackage_dataset` (ddf-pipeline) — stale `table.lock` stripped, data kept.

## Next — unit / CI

- [ ] `ms_quality` casacore test: a synthetic MS written in `tmp_path`, in a separate CI job with `python-casacore`.
- [ ] Audit the other `DicoNode2mslist[Node]` sites in `pipeline.py` (~195, 688, 921, 1240, 1385) for the same KeyError; route them through `utils/pernode.py`.
- [ ] `mpi_manager` unit tests: node scatter, `DicoNode2mslist`, `filterHost` with mocked ranks/hosts.

## Next — acceptance / cluster

- [ ] DDFacet TOCTOU: test that `ClassCasaImage.ToFits` tolerates an already-deleted output (casacore, Jean Zay).
- [ ] Migrate DDFacet/killMS tests to pytest (ADR-0005), regenerating reference data from the cleaned sim.
- [ ] Smoke assertion: the run directory contains the `DDF_monitor*` / `DDF_register*` / `DDF_io_monitor*` CSVs.
- [ ] Submit-dir test for the SLURM launchers (fake `SLURM_SUBMIT_DIR`).

## Next — end-to-end / multi-node

- [ ] `scripts/run_tests_jean_zay.sh` (1 node smoke, 2–3 nodes MPI).
- [x] 3-node e2e on the cleaned dataset — job `803303`, smoke `PASS` (locks + atomic save).
- [ ] Real 12-node run once multi-node is stable.

## Known risks that are not unit-testable

- Multi-node file races: only the e2e run proves them; unit tests prove the arithmetic/atomicity.
- `/dev/shm` and `vm.max_map_count` limits on the cluster (environmental).
- CINES/Adastra account access (infrastructure).
