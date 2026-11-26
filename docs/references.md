# Scientific references

The pipeline is the implementation of a body of radio-interferometry work. This page
records the papers, the invariant each pins, and whether the harness (ADR-0006) can test
it. Two groups: the **baseline** method papers, and the **recent / extended** set
(ddf-pipeline, LoTSS, LOFAR, SKA).

Tiers: **A** = `harness/fundamentals/` (pure numpy, CI), **B** = `harness/synthetic/`
(known-answer on the cluster), **C** = QC oracle (deferred).

## Baseline — the methods the code implements

| Paper | Pins | Layer | Harness value |
| --- | --- | --- | --- |
| **Smirnov & Tasse 2015**, MNRAS 449, 2668 (arXiv:1502.06974) | RIME as a complex optimisation (StEFCal), gain degeneracies | A + B | **high** — add a StEFCal step and a degeneracy test |
| **Tasse et al. 2017**, A&A 611, A87 (arXiv:1712.02078, DDFacet) | facet tesselation, `lmShift`/`lmExtent`, DD deconvolution | A | medium-high — facets tile the field; a source in a facet is recovered |
| **Tasse 2014**, A&A 566, A127 (arXiv:1410.8706) + **arXiv:1403.6308** | KAFCA Kalman filter; Wirtinger (half) Jacobian | A | medium — Kalman update on a gain track; Wirtinger gradient vs finite differences |
| **Tasse et al. 2021**, A&A 648, A1 | LoTSS Deep Fields DD calibration and imaging; production QC | C | high (deferred) |
| **Shimwell et al. 2019, 2022**, A&A 622, A1; A&A 659, A1 | LoTSS DR1/DR2 survey QC (flux, astrometry, completeness) | C | medium |

## Recent / extended — ddf-pipeline, LoTSS, LOFAR, SKA

| Paper | Pins | Layer | Harness value |
| --- | --- | --- | --- |
| **Certenais, Bodin & Morin 2025** (arXiv:2509.03075) | the DDF Pipeline itself: architecture, stages, validation | all | **very high** — the harness should mirror the stages this paper describes |
| **Shimwell, Hardcastle, Tasse et al. 2026** (LoTSS VII, DR3; arXiv:2602.15949) | latest LoTSS; DD calibration in ~45 facets; survey QC | C | high — the current survey reference |
| **van Weeren et al. 2016** (arXiv:1601.05422, LOFAR facet calibration) | the facet-calibration method that DDFacet builds on | A (facets) + C | medium-high |
| **Morabito et al. 2025**, Ap&SS 370, 19 (doi:10.1007/s10509-025-04406-x) | ILT sub-arcsecond DD calibration in practice; ionosphere, smearing, astrometry | C | medium |
| **de Jong et al. 2024** (ELAIS-N1 deep, A&A aa50595-24) | improved DI/DD calibration; deepest wide-field imaging | C | medium |
| **Wu et al. 2025**, AJ 170, 308 (doi:10.3847/1538-3881/ae10a9) | SKA1-LOW **IQA**: completeness, false-detection rate, position, flux; source-finder choice | C | **high** — the layer-C IQA oracle |
| **Xie et al. 2024**, AJ 167, 259 (doi:10.3847/1538-3881/ad3af5) | SKA1-LOW point-source sensitivity vs weighting and tapering | C + A | medium-high — pins the weighting/sensitivity relation |
| **Thekkeppattu, Wayth & Sokolowski 2024** (arXiv:2401.08039) | SKA-Low station calibration by holography | context | low for the pipeline harness (station level) |
| **Hardcastle et al. 2025** (arXiv:2504.09303) | radio AGN selection in LoTSS DR2 | context | low — survey science, not a calibration/imaging invariant |

## What is directly testable

- **Now (layer A)**: Smirnov & Tasse 2015 (StEFCal + degeneracy), Tasse 2017 and
  van Weeren 2016 (facets), Xie 2024 (weighting/sensitivity).
- **Layer C (deferred)**: Wu 2025 (IQA), Tasse 2021, Shimwell 2019/2022/2026,
  de Jong 2024, Morabito 2025.
- **Context only**: Thekkeppattu 2024 (station holography), Hardcastle 2025 (survey science).
- **The spec**: Certenais et al. 2025 describes the pipeline the harness validates.
