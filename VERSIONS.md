# Validated versions

The version set validated on Jean Zay. Update this file (and the per-repo tag) whenever the
validated set changes. `install_uv_jean_zay.sh` pins these.

| Repository | Reference | Notes |
| --- | --- | --- |
| `mhardcastle/ddf-pipeline` | `astra-hackathon` (`de84e2b`, per-node distribution extracted) | pipeline |
| `cyriltasse/DDFacet` | `675ce535` + `astra-hackathon` (`f495fd9`: TOCTOU unlink fix, atomic `NodesCat` save) | imaging |
| `cyriltasse/killMS` | `c8f17884` + `astra-hackathon` (local patch) | calibration |

Status:

- **1 node**: validated — `00:32:13`, smoke `PASS`.
- **3 nodes**: validated — smoke `PASS`, job `803303` on the cleaned dataset (no `table.lock`;
  atomic `NodesCat` save). Stats: mean 1.4 %, std 1.6 %, p1 5.9 %, p99 5.7 %, max 5.7 %.
  An earlier 3-node run (`750593`) also passed (`00:46:38`, MaxRSS ≈ 23 GB), but the same
  configuration failed intermittently on stale table locks and the non-atomic `NodesCat`
  write; both are fixed.
- **6 / 12 nodes**: not worth running on the simulation (overhead-bound). The real-dataset run
  is the meaningful multi-node test.
