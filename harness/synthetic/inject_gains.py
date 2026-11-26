"""Inject known antenna gains into a copy of a real MS (layer B, Jean Zay).

Multiplies every visibility by g_i conj(g_j) for a known set of complex gains, so
a subsequent KillMS run should recover them (up to the self-cal degeneracy). Uses
an existing MS as a template, so no table schema is created from scratch.

    python harness/synthetic/inject_gains.py \
        "$SCRATCH/ddf-sim-dataset/Sim/0000.p00.MS" /tmp/known.ms --col DATA

Writes /tmp/known.ms and /tmp/known.ms.true_gains.npy.
Casacore is required; run on Jean Zay, not in CI.
"""
import argparse
import shutil

import numpy as np


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("template_ms")
    p.add_argument("out_ms")
    p.add_argument("--col", default="DATA", help="column to scale (the one kMS reads)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--amp-spread", type=float, default=0.1)
    args = p.parse_args(argv)

    from casacore.tables import table

    shutil.copytree(args.template_ms, args.out_ms)
    t = table(args.out_ms, readonly=False, ack=False)
    try:
        a1 = t.getcol("ANTENNA1")
        a2 = t.getcol("ANTENNA2")
        data = t.getcol(args.col)
        n_ant = int(max(a1.max(), a2.max())) + 1

        rng = np.random.default_rng(args.seed)
        amp = 1.0 + args.amp_spread * rng.normal(size=n_ant)
        gains = amp * np.exp(2j * np.pi * rng.random(n_ant))

        factor = (gains[a1] * np.conj(gains[a2]))[:, None, None]
        t.putcol(args.col, data * factor)
    finally:
        t.close()

    np.save(args.out_ms + ".true_gains.npy", gains)
    print("wrote", args.out_ms, "-> true gains in", args.out_ms + ".true_gains.npy")


if __name__ == "__main__":
    main()
