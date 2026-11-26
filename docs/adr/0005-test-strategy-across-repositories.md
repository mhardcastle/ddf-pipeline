# Test strategy across the three repositories

Testing was ad hoc: DDFacet and killMS carry `nosetests` + `unittest.TestCase` golden-master
tests in Docker, while ddf-pipeline has only a statistical smoke check. We are starting from
scratch with **pytest** and a three-tier taxonomy, shared by all three repositories:

- **unit** — algorithm-level, no data, exact assertions. Runs on a laptop and in CI.
- **acceptance** — intermediate steps (imaging, weighting, `ClipCal`, solutions) compared to a
  reference; oracle is **pixel-exact** because it runs one process at a fixed thread count.
- **end-to-end** — the ddf-pipeline workflow and **multi-node MPI**; oracle is **statistical**
  (ADR-0003), because an iterative, threaded, MPI run is not pixel-reproducible.

We keep only the *concepts* of DDFacet's `ClassCompareFITSImage` (golden master, max squared
error, MSE, runtime, cycle count) — not its nose/unittest code, which we rewrite in pytest.

## Environments

Docker is forbidden on Jean Zay and Adastra, so the heavy tiers run there through **SLURM**
with a pinned **`uv`** install (`install_uv_jean_zay.sh`). Since we regenerate our own
reference images on that same stack, the pinned commits + lockfile give the reproducibility
the pixel-exact oracle needs; a **Singularity** image is documented as an optional
cross-site portability path, not a requirement.

| Tier | Oracle | Where | Reference data | Trigger |
| --- | --- | --- | --- | --- |
| unit | exact | laptop + GitHub Actions | none (synthetic arrays) | every push |
| acceptance (short) | pixel-exact | GitHub Actions | small MS + FITS, pinned SHA | every push |
| acceptance (long) | pixel-exact | Jean Zay / Adastra SLURM | small MS + FITS, pinned SHA | `run_tests_jean_zay.sh` |
| end-to-end | statistical | Jean Zay / Adastra SLURM | cleaned `Sim.tgz` | `run_tests_jean_zay.sh` |

Tests live in `cyriltasse/DDFacet`, `cyriltasse/killMS` and `mhardcastle/ddf-pipeline` on
their `astra-hackathon` branches, never on `master`.

## Considered options

- **Keep nose / `ClassCompareFITSImage`** — rejected: nose is deprecated, the existing DDFacet
  acceptance data lives on a private Jenkins server (`/var/lib/jenkins/test-data`), and little
  test effort exists to preserve. We rebuild in pytest and regenerate the reference data from
  the cleaned simulated dataset (we are its author, so we control the truth).
- **Singularity as the primary test runtime** — rejected for now: the pinned `uv` install
  already reproduces the stack, and we generate the references on it. Singularity is kept for
  portability between sites if that becomes necessary.

## Consequences

- CI cannot run the MPI or long tiers: those are manual SLURM jobs, so a green CI does not
  mean the pipeline is scientifically valid (see `docs/quality-control.md`).
- A pixel-exact acceptance test pins `NCPU` (DDFacet) or `OPENBLAS_NUM_THREADS=1` (killMS);
  changing either changes the reference.
- The reference MS/FITS are generated once on a validated stack and pinned by SHA256; a
  regenerated reference is a reviewed, deliberate change.
