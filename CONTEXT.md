# CONTEXT.md — ddf-pipeline

## What this is

ddf-pipeline reduces LOFAR (continuum, Stokes I) imaging data to a science-ready state
through several iterations of direction-dependent self-calibration, using DDFacet
(imaging) and KillMS (calibration). It is used by the LOFAR Surveys KSP for the LoTSS
survey. Written by Hardcastle, Shimwell, Tasse and Williams; GPL-licensed.

## Architecture

`scripts/pipeline.py` orchestrates a loop: it drives DDFacet (imaging) and KillMS
(calibration) over several self-calibration iterations. Inputs are a measurement set
(`.ms`), a configuration file (`.cfg`), and list files (`mslist.txt`, `big-mslist.txt`).
Outputs are FITS imaging products and modeled measurement sets.

## Layout

| Path | Purpose |
| --- | --- |
| `scripts/` | Orchestration (`pipeline.py`), stage scripts, `compare_runs.py` |
| `utils/` | Parameter parsing (`parset.py`), helpers (`auxcodes.py`), options, quality |
| `examples/` | `.cfg` configuration files |
| `slurm/` | Cluster launchers (`jean-zay.slurm`, `adastra.slurm`) and their shared `env.sh` |
| `tests/` | pytest unit tests |
| `.github/` | GitHub Actions CI (`.github/workflows/ci.yml`) |
| `others/` | Default parset, survey DB |
| `misc/` | `DDF.sh`, pointing lists |
| `docs/` | Manual, developer notes, datasets (`data.md`), ADRs, agent config |

## Build and run

Install (on a machine with Python ≥3.10 and uv):

```bash
./install_uv_jean_zay.sh   # modules + pinned DDFacet/killMS + uv sync --extra mpi-support
./install_data_jean_zay.sh # simulated "small" dataset + catalogs (SHA256-pinned)
# ./install_data_big.sh    # LoTSS "big" dataset for a realistic run (not CI)
source init_dev_venv.sh    # activates the venv, sets PYTHONPATH / LD_LIBRARY_PATH
```

Run (Jean Zay, simulated dataset; the full run needs ~40 min):

```bash
sbatch --time=00:40:00 slurm/jean-zay.slurm   # 1 MPI rank, 40 threads → scripts/pipeline.py examples/Simul-Hackaton.cfg
sbatch --time=00:40:00 slurm/adastra.slurm    # same on Adastra GENOA: 1 MPI rank, 192 threads
```

Entry point: `scripts/pipeline.py <config.cfg>` (also registered as the `pipeline.py`
console script).

## Branches and releases

Code lives in three canonical repositories (`mhardcastle/ddf-pipeline`, `cyriltasse/DDFacet`,
`cyriltasse/killMS`). Work flows `astra-hackathon` → `HackathonRennes_no_futures` →
`HackathonRennes`; `master` is not the integration point. **`master` is production and is
never pushed**; CI is on the exploratory `HackathonRennes` lineage. The smoke test runs on a
validated tag. See `docs/adr/0004-branching-and-release-model.md` and `VERSIONS.md`.

## Configuration

Parameters live in `utils/parset.py` (`option_list`), grouped by category (`machine`,
`data`, `solutions`, `image`, …). A `.cfg` file sets them; they surface in the `o`
options dict.

## Environment

- Requires Python ≥3.10, uv, and sibling checkouts of DDFacet and KillMS.
- DDFacet is pinned to `675ce535`, KillMS to `c8f17884` (`install_uv_jean_zay.sh`).
- Heavy dependencies: `python-casacore`, `everybeam`, `ddfacet`, `killms`.
- MCP servers (`ddf`, `hpc`) are **lab-network-only**; treat them as optional, never required.

## Docs

- `docs/manual.md` — user manual (legacy, being superseded by the mkdocs site).
- `docs/data.md` — datasets, downloads and the first run.
- `docs/developer_notes.md` — developer notes (legacy).
- `GLOSSARY.md` — domain vocabulary.
- `docs/adr/` — architectural decisions.
