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
- [ ] Audit killMS/DDFacet for other shared `SaveFile` / read-after-write races (the `ClipCal.last` and `NodesCat` classes).

## Fundamentals harness (ADR-0006)

- [x] Layer A (`harness/fundamentals/`): Fourier (2), visibility/VZ (4), imaging/PSF/weighting (5), deconvolution (6), Jones/RIME/closure (7), calibration (8) — 15 tests, CI.
- [ ] Layer B (`harness/synthetic/`): tiny synthetic MS (python-casacore) + known gains; run KillMS/DDFacet; known-answer oracle on Jean Zay.
- [ ] Layer C (deferred): QC oracle injecting known sources/gains into the sim.

## Shared-file audit (concurrent multi-node paths)

Fixed: `NodesCatAll/Imaged.npy`, `DicoFacet`, `facetCoord.txt` (DDFacet, atomic save); `ClipCal.last` (killMS, no re-read); FITS image output (`ToFits`, atomic write).

Still to audit / decide:
- [ ] Monitors: `last_IOMonitor.obj` / `last_MemMonitor.obj` (per-node, fixed `SaveFile` write+read).
- [ ] `last_DDFacet.finalised.obj` in `ddf_image`'s MPI branch (written/read by each rank if `ddf_image` runs on all ranks).
- [ ] `DicoModel` (`ClassModelMachine.ToFile`): verify it is written once (rank 0) or make it atomic.
- [ ] killMS `SmoothSols` / `InterpSols` / `MergeSols` / `AQWeight` fixed `SaveFile`s, if run per node.

Benign in the concurrent path (written, not read): `last_killMS.obj`, `last_DDFacet.obj`, `ddf-pipeline.last`.
- [ ] `mpi_manager` unit tests: node scatter, `DicoNode2mslist`, `filterHost` with mocked ranks/hosts.

## Next — acceptance / cluster

- [ ] DDFacet TOCTOU: test that `ClassCasaImage.ToFits` tolerates an already-deleted output (casacore, Jean Zay).
- [ ] Migrate DDFacet/killMS tests to pytest (ADR-0005), regenerating reference data from the cleaned sim.
- [ ] Smoke assertion: the run directory contains the `DDF_monitor*` / `DDF_register*` / `DDF_io_monitor*` CSVs.
- [ ] Submit-dir test for the SLURM launchers (fake `SLURM_SUBMIT_DIR`).

## Next — end-to-end / multi-node

- [ ] `scripts/run_tests_jean_zay.sh` (1 node smoke, 2–3 nodes MPI).
- [x] 3-node e2e on the cleaned dataset — job `803303`, smoke `PASS` (locks + atomic save).
- [x] 6-node e2e — smoke `PASS` after the ClipCal / DicoFacet / FITS atomic fixes.
- [ ] 12-node e2e.
- [ ] Real 12-node run on `lotss-big` once multi-node is stable.

## Known risks that are not unit-testable

- Multi-node file races: only the e2e run proves them; unit tests prove the arithmetic/atomicity.
- `/dev/shm` and `vm.max_map_count` limits on the cluster (environmental).
- CINES/Adastra account access (infrastructure).
