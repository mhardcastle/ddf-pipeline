# Validated versions

The version set validated on Jean Zay. Update this file (and the per-repo tag) whenever the
validated set changes. `install_uv_jean_zay.sh` pins these.

| Repository | Reference | Notes |
| --- | --- | --- |
| `mhardcastle/ddf-pipeline` | `astra-hackathon` (KeyError fix `de84e2b`) | pipeline |
| `cyriltasse/DDFacet` | `675ce535` + `astra-hackathon` (TOCTOU unlink fix) | imaging |
| `cyriltasse/killMS` | `c8f17884` + `astra-hackathon` | calibration |

Status:

- **1 node**: validated — `00:32:13`, smoke `PASS`.
- **3 nodes**: validated — `00:46:38`, smoke `PASS`, MaxRSS ≈ 23 GB. **Negative scaling**: 3
  nodes are *slower* than 1, because the global DI/Predict steps are not distributed.
- **6 / 12 nodes**: not worth running on the simulation (overhead-bound). The multi-node
  path is proven at 3 nodes; the real-dataset run is the meaningful test.
