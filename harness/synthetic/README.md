# Layer B — synthetic calibration harness

A **known-answer** test of the solver: inject known complex gains into real data,
run KillMS 1GC on the same model, and check what it recovers. Unlike layer A
(`harness/fundamentals/`, pure numpy, CI), this needs casacore and a completed
smoke run, so it runs on **Jean Zay**.

## Idea

A point-like gain on every baseline multiplies the visibilities by
`g_i conj(g_j)`. If KillMS solves the same field, its recovered gains must equal
the injected ones **up to a global complex factor** (the self-cal degeneracy).
The check compares gains normalised by antenna 0 (`check_sols.compare_gains`).

## Procedure (Jean Zay)

1. **Inject gains** into a copy of one sim MS, into the column KillMS reads:

   ```bash
   python harness/synthetic/inject_gains.py \
       "$SCRATCH/ddf-sim-dataset/Sim/0000.p00.MS" /tmp/known.ms --col DATA
   # -> /tmp/known.ms and /tmp/known.ms.true_gains.npy
   ```

2. **Solve** with the same `kMS.py` command the pipeline uses (see
   `killms_data_serial` in `scripts/pipeline.py`), with `--MSName /tmp/known.ms`,
   `--InCol DATA`, and the `--DicoModel` produced by a completed smoke run for the
   same image. Keep the solver settings identical so the test is about recovery,
   not tuning.

3. **Check** the recovered solutions against the injected gains:

   ```bash
   python harness/synthetic/check_sols.py \
       /tmp/known.ms/killMS.<outsols>.sols.npz \
       /tmp/known.ms.true_gains.npy
   # prints "max |recovered - true| ..." then PASS/FAIL
   ```

## Status

- [x] `visibility_model.py` — the gain model (pure numpy, CI-tested).
- [x] `check_sols.py` — `gains_from_G` / `compare_gains` (pure, CI-tested).
- [ ] `inject_gains.py` — casacore; **validate on Jean Zay**.
- [ ] Wire steps 1–3 into `scripts/run_tests_jean_zay.sh`.

The sols format is `Sols.G` of shape `(nt, nch, na, nd, 2, 2)`; scalar mode uses
the `(0,0)` element. If the pinned KillMS changes that, update `gains_from_G`.
