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
| `scripts/` | Orchestration (`pipeline.py`) and the pipeline stage scripts |
| `utils/` | Parameter parsing (`parset.py`), helpers (`auxcodes.py`), options, quality |
| `examples/` | `.cfg` configuration files |
| `others/` | Default parset, survey DB |
| `misc/` | `DDF.sh`, pointing lists |
| `docs/` | Manual, developer notes, ADRs, agent config |

## Build and run

Install (on a machine with Python ≥3.10 and uv):

```bash
./install_uv.sh          # installs uv, clones DDFacet + killMS, uv sync --extra mpi-support
./prepare_sim_exe.sh     # downloads catalogs + the simulated dataset + Simul-Hackaton.cfg
source init_dev_venv.sh  # activates the venv, sets PYTHONPATH / LD_LIBRARY_PATH
```

Run (Jean Zay, simulated dataset):

```bash
sbatch sim_exe.slurm     # 40 cpus, 6 tasks, 20 min → scripts/pipeline.py examples/Simul-Hackaton.cfg
```

Entry point: `scripts/pipeline.py <config.cfg>` (also registered as the `pipeline.py`
console script).

## Configuration

Parameters live in `utils/parset.py` (`option_list`), grouped by category (`machine`,
`data`, `solutions`, `image`, …). A `.cfg` file sets them; they surface in the `o`
options dict.

## Environment

- Requires Python ≥3.10, uv, and sibling checkouts of DDFacet and KillMS.
- Heavy dependencies: `python-casacore`, `everybeam`, `ddfacet`, `killms`.
- MCP servers (`ddf`, `hpc`) are **lab-network-only**; treat them as optional, never required.

## Docs

- `docs/manual.md` — user manual (legacy, being superseded by the mkdocs site).
- `docs/developer_notes.md` — developer notes (legacy).
- `GLOSSARY.md` — domain vocabulary.
- `docs/adr/` — architectural decisions.
