# Validated versions

The version set validated on Jean Zay. Update this file (and the per-repo tag) whenever the
validated set changes. `install_uv_jean_zay.sh` pins these.

| Repository | Reference | Notes |
| --- | --- | --- |
| `mhardcastle/ddf-pipeline` | `astra-hackathon` (`de53885`: per-node dict; `3671241`: pkill non-fatal; matrix copy/config hardening) | pipeline |
| `cyriltasse/DDFacet` | `675ce535` + `astra-hackathon` (`f495fd9`: NodesCat atomic; `3ce3799`: DicoFacet/facetCoord atomic; `acdc8b0`: FITS atomic) | imaging |
| `cyriltasse/killMS` | `c8f17884` + `astra-hackathon` (`554b0d9`: ClipCal shared `SaveFile` race) | calibration |

Status:

- **1 node**: validated — `00:32:13`, smoke `PASS`.
- **3 nodes**: validated — job `803303`, `00:34:37`, MaxRSS ≈ 28.5 GB, smoke `PASS` on the
  cleaned dataset (no `table.lock`; atomic `NodesCat` save). Stats: mean 1.4 %, std 1.6 %,
  p1 5.9 %, p99 5.7 %, max 5.7 %. An earlier 3-node run (`750593`) also passed
  (`00:46:38`, MaxRSS ≈ 23 GB), but the same configuration failed intermittently on stale
  table locks and the non-atomic `NodesCat` write; both are fixed. 3 nodes ≈ 1 node
  (`00:32:13`): no scaling gain, as expected for the overhead-bound simulation.
- **6 nodes**: validated — smoke `PASS` after the `ClipCal`/`DicoFacet`/FITS fixes (image stats
  ≤3.5 % of the reference). A first 6-node attempt (job `809116`) had failed on the `ClipCal`
  shared `SaveFile` race, fixed in killMS `554b0d9`.
- **12 nodes**: validated — 12a and 12b both smoke `PASS` (mean 4.6 % / 1.9 % of the reference;
  all statistics within tolerance). The run exercises the collective imaging across all 12 ranks
  and exits cleanly (the `pkill` no-match is no longer fatal).
