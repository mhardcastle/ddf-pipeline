# A fundamentals harness with known answers

The pipeline's correctness rests on interferometry fundamentals — the Fourier/imaging chain,
the RIME, calibration, deconvolution — yet the existing tests check the *plumbing* (keys,
files, statistics), not those laws. We add a **harness** that encodes the fundamentals as
executable checks with analytic or known answers, so a physics regression is caught.

`o-smirnov/fundamentals_of_interferometry` (the NASSP course) is the reference for the
concepts. It is pedagogical and old (Python 2.7, external KAT-7 data), so we do **not** reuse
its notebooks or data: we re-derive the concepts as modern, self-contained numpy tests.

## Shape

`harness/` in `ddf-pipeline`, two layers:

- **`harness/fundamentals/` (layer A)** — pure-numpy pytest modules, one per book section that
  matters: Fourier/convolution/sampling (2), visibility and Van Cittert-Zernike (4),
  imaging/PSF/weighting (5), deconvolution (6), Jones/RIME and closure (7), calibration (8).
  Exact/analytic assertions; runs in CI.
- **`harness/synthetic/` (layer B)** — a tiny synthetic Measurement Set (python-casacore) with
  a known sky and known gains; run KillMS/DDFacet and verify recovery with a known-answer
  oracle (statistical tolerance, ADR-0003). Runs on Jean Zay.

A QC-oracle layer (C) is deferred.

## Considered options

- **Reuse the book's notebooks/data** — rejected: pedagogical, Python 2.7, unpinned external
  tarballs; the lessons are the value, not the code.
- **Split the harness across the three repos** — rejected: the fundamentals are cross-cutting;
  duplicating them invites divergence.

## Consequences

- CI gains a physics layer (numpy only); the synthetic layer needs python-casacore and the
  cluster.
- The harness is **generic** interferometry (KAT-7/VLA/SKA examples), not LOFAR/DDF-specific:
  it complements, and does not replace, the pipeline's own tests and the ADR-0005 tiers.
- A fundamentals failure is a **physics** regression; a tier failure is a plumbing/behaviour
  regression.
