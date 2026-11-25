# Validated versions

The version set validated on Jean Zay. Update this file (and the per-repo tag) whenever the
validated set changes. `install_uv_jean_zay.sh` pins these.

| Repository | Reference | Notes |
| --- | --- | --- |
| `mhardcastle/ddf-pipeline` | `astra-hackathon` (KeyError fix `de84e2b`) | pipeline |
| `cyriltasse/DDFacet` | `675ce535` + `astra-hackathon` (TOCTOU unlink fix) | imaging |
| `cyriltasse/killMS` | `c8f17884` + `astra-hackathon` | calibration |

Status:

- **1 node**: validated — 32 min, smoke `PASS`.
- **3 / 6 nodes**: still failing (exit 143, to re-measure with the DDFacet fix).
- **12 nodes**: `KeyError` fixed; `ToFits` TOCTOU race fixed in DDFacet; the run currently
  exceeds the 40-minute wall-time, so it needs a larger `--time`.
