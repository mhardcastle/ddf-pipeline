# The smoke oracle compares image statistics, not pixels

Two identical runs of the simulated dataset on Jean Zay produce pixel-different restored
images (RMS difference ≈ 35% of the image RMS), while their global statistics agree to
about 1% (mean, std) and 5% (1st/99th percentiles). The pipeline is therefore not
pixel-reproducible, so the smoke test asserts that the named artefacts exist, that the
logs contain no `Traceback`, and that the restored image's global statistics match
`tests/data/smoke_reference.json` within per-statistic tolerances.

Considered options: a pixel/hash comparison (rejected — it fails on every run because of
threading, FFT and iterative-solver non-determinism); an exit-code-only check (rejected —
too weak to catch a broken image). Consequence: the reference must be re-pinned
deliberately when the pipeline's output changes; the tolerances (mean/std 5%, p1/p99 20%,
max 10%) were set from the two-run spread and should only be widened with evidence.
