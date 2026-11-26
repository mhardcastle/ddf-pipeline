# Quality control

The primary output of ddf-pipeline is the **calibrated Measurement Set(s)** — the
direction-dependent gain solutions and the corrected visibilities. The restored FITS images
are a by-product. Quality control therefore targets the **visibilities and solutions first**,
and the images second.

The smoke test proves the plumbing runs end to end; it says nothing about flux scale,
astrometry, dynamic range or convergence (see `GLOSSARY.md`). This page defines the metrics
that turn "the run finished" into "the run produced defensible science", and where each one
is implemented.

## Principle: what is being measured

Two kinds of check, with different evidential value:

- **Internal / self-consistency** — no external truth needed (residual RMS, solution
  statistics, closure, flags, weights). Cheap, scriptable from one MS or `.npz`.
- **Absolute** — needs a truth (a reference catalogue, or the injected model of the
  simulated dataset). This is where flux scale and astrometry live.

The simulated dataset (`Sim.tgz`, driven by `examples/Simul-Hackaton.cfg`) injects a known
model: `scripts/simulation-pipeline.py` predicts a known image into visibilities and adds it
to real data, so the injected sources are the ground truth. Metrics that depend on truth can
be **pinned against the injected model**; the real sky underneath has unknown true gains, so
its solutions and residuals can only be checked for *stability*, not absolutes.

## Calibrated-MS metrics (primary)

The residual is the master quantity: the difference between the corrected visibilities and
the model. In this pipeline the corrected column is `o['colname']` (normally
`DATA_DI_CORRECTED`; `DATA_SUB` when `do_wide`, `SCALED_DATA` when `bootstrap` — see
`docs/colname.txt`) and the model column is `DD_PREDICT`. `ClipCal.py` already forms
`DATA_DI_CORRECTED − DD_PREDICT`.

| # | Metric | Definition | Initial threshold | Detects / cannot detect |
| --- | --- | --- | --- | --- |
| 1 | **Normalised residual RMS** | `std( (colname − DD_PREDICT) · √w )` over unflagged visibilities, per correlation. A correctly calibrated, correct model gives `≈ 1` (Re/Im `≈ 1/√2` each). | `≤ 2.0` (initial; pin against the sim) | Diverged calibration, unmodeled sources, RFI, bad flags/weights. **Cannot** detect a global flux-scale error (see caveat 1). |
| 2 | **Solution amplitude** | Median and max of `|g|` per antenna/direction from `Sols.G` in `killMS.*.sols.npz`, plus the per-antenna scale in `gu.npy`. | median `∈ [0.7, 1.4]`, max `∈ [0.5, 2.0]`; the pipeline already uses `CutGainsMinMax 0.1,1.5` | Dead/broken antennas, non-convergence, a diverged direction. **Cannot** detect a coherent global scale (selfcal ambiguity). |
| 3 | **Solution phase** | Unwrapped `arg(g)` continuity, `±2π` wrap count, per-facet slope. | smooth, wrap-free (empirical) | Ionospheric wraps, one bad facet, solver non-convergence. **Cannot** detect common-mode phase. |
| 4 | **Closure phase** | `C_ijk = arg(V_ij) + arg(V_jk) − arg(V_ik)`; invariant to antenna-based gains. | use as **stability** (unchanged across iterations) and **residual** (scatter → thermal), never "≈ 0" | Baseline/non-closing errors, bandpass, model error. **Cannot** detect *any* antenna-based error (by construction). |
| 5 | **Closure amplitude** | `A_ijkl = |V_ij|·|V_kl| / (|V_ik|·|V_jl|)` on a bright compact source; `≈ 1` for a point source. | `|A − 1|` scatter ≈ thermal (empirical) | Unmodeled flux, baseline amplitude errors. **Cannot** detect antenna-based amplitude errors. |
| 6 | **Flag fraction** | `FLAG` fraction per channel/baseline/time. | no channel > 90 %; consistent with reference run (**report, don't score**) | RFI, dead antennas/baselines, broken flag step. **Cannot** grade calibration quality. |
| 7 | **Weight distribution** | Statistics of `WEIGHT`/`WEIGHT_SPECTRUM` (`σ_i = 1/√w`). | no outliers, no unintended zero weights | Mis-set weights, a visibility that should be flagged. **Cannot** measure calibration accuracy (weight is a *proxy*). |
| 8 | **Model-vs-data amplitude** | Observed/model apparent-flux ratio vs uv distance (`utils/histmsamp.py`). | ratio `≈ 1` over the well-covered uv range; `~10 %` flux-scale goal | Gross flux-scale errors, missing short baselines. **Cannot** detect sub-10 % flux errors. |
| 9 | **Solution coverage** | Fraction of (antenna, direction, time, freq) slots with a finite solution; NaN/inf count. | 100 % coverage, 0 NaN (deterministic) | Non-convergence, dropped directions, a silently-crashed solver. **Cannot** tell if a converged solution is *right*. |

## Image metrics (by-product)

Most of these are already implemented in `scripts/quality_pipeline.py` and its helpers
(`utils/dr_checker.py`, `utils/fit_sourcecounts.py`, `utils/cent-rms.py`).

| # | Metric | Definition | Initial threshold | Source |
| --- | --- | --- | --- | --- |
| 1 | Central RMS noise | std in source-free regions with outlier rejection (`utils/auxcodes.get_rms`). | ≤ ~100 µJy/beam for a full pointing; sim vs `tests/data/smoke_reference.json` | Shimwell et al. 2022 |
| 2 | Flux scale | median `S_LOFAR / S_ref` over compact isolated sources; sim: recovered/injected ratio. | `∈ [0.8, 1.2]` (sim truth), `∈ [0.7, 1.3]` (catalogue) | Scaife & Heald 2012; Shimwell et al. 2022 |
| 3 | Astrometry | median RA/Dec offset vs reference or injected positions. | `< 0.2″` (survey), `< 1.5″` (smoke bound) | Shimwell et al. 2022 |
| 4 | Dynamic range | peak flux / near-source noise (`utils/dr_checker.py`). | median > 50 (source/field dependent) | Shimwell et al. 2022, 2025 |
| 5 | Completeness / purity | recovered source counts vs the 150 MHz counts; completeness of injected sources. | 90 % at 0.8 mJy/beam (survey) | Intema et al. 2011; Mohan & Rafferty 2015 |
| 6 | Artefact level | spurious sources (sidelobe ghosts) in the residual map. | empirical density threshold | Tasse et al. 2021 |

## Caveats (metrics that mislead in this context)

1. **Self-calibrated residual RMS is biased (over-fitting).** Selfcal derives the model from
   the data, so `DATA − MODEL` can be driven to thermal noise even with a wrong model or
   wrong flux scale — the error is absorbed into the gains (Smirnov & Tasse 2015). Residual
   RMS alone can false-pass; pair it with an absolute check (flux scale, or injected-source
   recovery).
2. **Closure is blind to antenna-based gains by construction** (Jennison 1958; Readhead &
   Wilkinson 1978). It cannot detect a global gain/flux error — that is precisely what it is
   invariant to. Use it for (i) stability across iterations and (ii) residual-closure scatter,
   never as "phase ≈ 0".
3. **"Solution amplitude ≈ 1" is necessary, not sufficient.** It does not rule out a wrong
   flux scale or sky model; the solver converges to an internally-consistent but possibly
   wrong answer. Convergence ≠ correctness.
4. **Flag fraction is not a quality score.** A heavily-flagged, correctly-calibrated MS is
   better than a lightly-flagged, wrongly-calibrated one. Report it; only fail on catastrophic
   per-channel/baseline flagging.
5. **Weight is a proxy, not a measurement.** `σ_i = 1/√w` encodes the *assumed* noise, not the
   achieved calibration quality.
6. **Dynamic range is source- and definition-dependent.** Use Shimwell et al. 2022's
   definition (noise increase near a source); `dr_checker.py`'s median is a within-run
   comparator, not an absolute grade.
7. **Central RMS in a field with bright sources is not thermal noise.** Use source-free
   regions and outlier rejection, or artefacts contaminate the number.
8. **Single-pointing source counts are biased** by the bright population and survey area;
   meaningful only with a defined area and flux scale (and, on the sim, the injected list).

## Recommended minimal set (smoke)

Extend `smoke_check.py`'s artefact/log/image-stat checks with, in value-per-effort order
(all numeric bounds are **calibration targets to fix on the sim**, not published constants):

1. Normalised residual RMS `≤ 2.0`.
2. Solution amplitude `∈ [0.7, 1.4]` median / `∈ [0.5, 2.0]` max, and 100 % coverage, 0 NaN.
3. Flag fraction: no channel > 90 % (report, don't score).
4. Flux scale vs injected truth `∈ [0.8, 1.2]`.
5. Astrometry vs injected truth `< 1.5″`.
6. Central RMS vs `smoke_reference.json` (already wired).
7. Dynamic range median > 50.
8. Closure-phase stability across iterations (optional, cheap).

## Where each metric lives

- **MS-level** (1, 4–9): `scripts/ms_quality.py` (new) — reads the MS and the
  `killMS.*.sols.npz`/`gu.npy`, emits a JSON summary with the same tolerance machinery as the
  smoke check.
- **Image RMS** (image #1): `smoke_check.py` (already) against `tests/data/smoke_reference.json`.
- **Flux scale, astrometry, completeness, artefacts** (image #2–#6): `scripts/quality_pipeline.py`
  + catalogues — release-level, not CI (needs `bdsf`, external catalogues, `SurveysDB`).

## References

1. Smirnov & Tasse 2015, MNRAS 449, 2668 — https://doi.org/10.1093/mnras/stv418
2. Tasse et al. 2018, A&A 611, A87 — https://doi.org/10.1051/0004-6361/201732474
3. Jennison 1958, MNRAS 118, 276 — https://doi.org/10.1093/mnras/118.3.276
4. Readhead & Wilkinson 1978, ApJ 223, 25 — https://ui.adsabs.harvard.edu/abs/1978ApJ...223...25R/abstract
5. Shimwell et al. 2022, A&A 659, A1 — https://doi.org/10.1051/0004-6361/202142484
6. Tasse et al. 2021, A&A 648, A1 — https://doi.org/10.1051/0004-6361/202038804
7. Scaife & Heald 2012, MNRAS 423, L30 — https://doi.org/10.1093/mnrasl/sls053
8. Intema et al. 2011, A&A 535, A38 — https://doi.org/10.1051/0004-6361/201014253
9. Mohan & Rafferty 2015, PyBDSF, ASCL 1502.007 — https://ui.adsabs.harvard.edu/abs/2015ascl.soft02007M/abstract
