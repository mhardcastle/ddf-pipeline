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
