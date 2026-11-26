# Testing

Three tiers, shared by the three repositories (decision: ADR-0005):

| Tier | What it covers | Oracle | Runs on |
| --- | --- | --- | --- |
| **unit** | algorithm-level, no data | exact assertions | laptop + GitHub Actions |
| **acceptance** | intermediate steps vs a stored reference | pixel-exact (fixed threads) | CI (short) / Jean Zay + Adastra SLURM (long) |
| **end-to-end** | full workflow + multi-node MPI | statistical (ADR-0003) | Jean Zay + Adastra SLURM |

Everything is **pytest**. Docker is used only on laptops and in GitHub Actions; Jean Zay and
Adastra run the heavy tiers through SLURM with the pinned `uv` install
(`install_uv_jean_zay.sh`). A Singularity image is an optional portable alternative.

## Layout — ddf-pipeline

```
tests/                         # current (flat); target is tests/unit/ + tests/e2e/
  test_parset.py               # existing: option defaults/types
  test_ms_quality.py           # new: pure-numpy metric functions (residual, flags, weights, sols)
  data/smoke_reference.json    # e2e golden statistics
  # planned:
  # unit/test_mpi_manager.py   # node scatter, DicoNode2mslist, filterHost (mocked, no real MPI)
  # e2e/                       # smoke_check.py + smoke_reference.json
```

The e2e run will be submitted with `scripts/run_tests_jean_zay.sh` (planned; SLURM, 1 node for
the smoke, 2–3 nodes for MPI).

## Layout — DDFacet and killMS

The existing nose/`unittest` tests are rewritten in pytest; only their concepts and metric
vocabulary are kept.

```
DDFacet/Tests/
  unit/                     # was FastUnitTests: Fitter, StokesConverter, libraries
  acceptance/               # was ShortAcceptanceTests + VeryLongAcceptanceTests
  <TestName>.parset.cfg
  <TestName>.<product>.fits # generated reference (pinned, not committed)

killMS/TestHarness/
  unit/                     # new: solver, gains, ClipCal, kMS.py
  acceptance/               # was LongAcceptanceTests
```

An acceptance test is one process at a fixed thread count (`NCPU` for DDFacet,
`OPENBLAS_NUM_THREADS=1` for killMS) compared against a reference generated on the same
pinned stack; regenerate the reference deliberately and re-pin its SHA256.

## Running

```bash
# unit (anywhere)
pytest tests/unit -q

# acceptance (CI or a SLURM job)
DDFACET_TEST_DATA_DIR=/path/ref DDFACET_TEST_NCPU=8 pytest DDFacet/Tests/acceptance -q

# end-to-end (SLURM only)
scripts/run_tests_jean_zay.sh            # smoke, 1 node   (planned)
scripts/run_tests_jean_zay.sh --nodes 3  # MPI smoke       (planned)
```

A green unit/acceptance run proves the code and its conventions, not the science: see
`docs/quality-control.md`.
