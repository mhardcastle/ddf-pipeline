# ddf-pipeline

Direction-dependent self-calibration and imaging of LOFAR data using DDFacet and KillMS.
This glossary fixes the vocabulary for the domain; keep definitions free of implementation
detail.

## Language

### Instruments and surveys

**LOFAR**: The Low-Frequency Array, a radio interferometer.

**LoTSS**: The LOFAR Two-metre Sky Survey; the survey ddf-pipeline reduces.

**DR2**: LoTSS Data Release 2.

**KSP**: Key Science Project.

### Calibration and imaging

**DDF**: Direction-dependent — effects that vary across the field of view.
_Avoid_: DD (unless the abbreviation is defined nearby)

**Self-calibration (selfcal)**: Iteratively calibrating a field using its own data.

**Direction**: A sky position solved independently during self-calibration.

**Facet**: A sub-region of the field of view, imaged and calibrated separately.

**DI calibration**: Direction-independent calibration.

**Solutions**: The calibration gains produced by the solver.

**DDFacet**: The direction-dependent imaging tool (Tasse et al. 2017).

**KillMS**: The solver that produces calibration solutions (Tasse 2014; Smirnov & Tasse 2015).

**Stokes I**: Total intensity.

**Continuum**: Broadband emission, as opposed to spectral line.

### Data and files

**Measurement Set (.ms)**: The CASA format holding visibilities.

**Modeled Measurements Dataset (.ms)**: A measurement set output by the pipeline, with the sky model applied.

**Configuration (.cfg)**: The parameter-set file that drives a run.
_Avoid_: parset file (parset is the parser module, not the file)

**mslist**: A list of measurement sets to process (`mslist.txt`, `big-mslist.txt`).

**FITS**: Flexible Image Transport System; the imaging products.

### Testing

**Smoke test**: A run that proves the plumbing works end to end — `pipeline.py` completes on
the simulated dataset and produces its artefacts. It is not a scientific validation: it says
nothing about flux scale, astrometry, dynamic range or convergence.
_Avoid_: integration test (until the inputs are pinned and a truth-based assertion exists)

**Unit test**: An algorithm-level test with no data, exact assertions, sub-second.

**Acceptance test**: An intermediate-step test (imaging, weighting, ClipCal, solutions)
compared against a stored reference, with a pixel-exact oracle because it runs one process at
a fixed thread count.

**End-to-end test**: A test of the whole workflow including multi-node MPI; its oracle is
statistical, because an iterative, threaded, MPI run is not pixel-reproducible.

**Golden master**: A stored reference output (an image, a solution set) that a run is
compared against.

**Oracle**: The criterion that decides pass from fail: *pixel-exact* (acceptance) or
*statistical* (end-to-end, ADR-0003).

### Quality control

**Calibration residual**: The difference between the corrected visibilities and the sky
model (`colname − DD_PREDICT`). Its RMS, normalised by the visibility weights, is the
primary measure of calibration quality.

**Closure phase**: The sum of visibility phases around an antenna triangle; independent of
all antenna-based gains, so it reveals only baseline-dependent (non-closing) errors.

**Closure amplitude**: The ratio of visibility amplitudes around a four-antenna loop; for a
point source it is unity, so its departure from unity reveals baseline amplitude errors.

**Dynamic range**: The ratio of a source's peak flux to the noise or artefact level in its
vicinity; the image-side signature of residual calibration errors.

**Flux scale**: The absolute scaling between recovered and true flux; an absolute check that
self-calibration cannot verify on its own (the ambiguity is absorbed into the gains).

**Astrometric offset**: The angular offset between a recovered source position and its true
or reference position.

**Completeness**: The fraction of true sources recovered above a flux threshold.
**Purity**: The fraction of detected sources that are real rather than artefacts.

**Source counts**: The differential number of sources per unit area per unit flux; compared
to the Euclidean-normalised reference to flag flux-scale or completeness errors.

**Thermal noise**: The radiometer noise floor; the level a correct calibration and model
should leave in the residual.
