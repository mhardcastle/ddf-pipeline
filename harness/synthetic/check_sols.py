"""Compare KillMS 1GC solutions against the known gains (layer B).

The recovered gains are only defined up to a global complex factor (self-cal
degeneracy), so the check compares the gains normalised by antenna 0.
"""
import argparse

import numpy as np


def gains_from_G(G):
    """Recovered scalar gains (na,) from a killMS Sols.G array (nt, nch, na, nd, 2, 2)."""
    G = np.asarray(G)
    if G.ndim == 6:
        g = G[..., 0, 0][:, :, :, 0]  # (nt, nch, na), direction 0
    elif G.ndim == 4:
        g = G[..., 0]
    else:
        raise ValueError("unexpected Sols.G shape %r" % (G.shape,))
    flat = g.reshape(-1, g.shape[-1])
    return np.median(flat, axis=0)


def compare_gains(g_hat, g_true):
    """Max |g_hat/g_hat[0] - g_true/g_true[0]| (global phase removed)."""
    g_hat = np.asarray(g_hat, dtype=complex)
    g_true = np.asarray(g_true, dtype=complex)
    return float(np.max(np.abs(g_hat / g_hat[0] - g_true / g_true[0])))


def load_solutions(sols_npz):
    d = np.load(sols_npz, allow_pickle=False)
    return gains_from_G(d["Sols"]["G"])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sols_npz")
    p.add_argument("true_gains_npy")
    p.add_argument("--tol", type=float, default=0.1)
    args = p.parse_args(argv)

    g_true = np.load(args.true_gains_npy)
    g_hat = load_solutions(args.sols_npz)
    err = compare_gains(g_hat, g_true)
    print("max |recovered - true| (global phase removed) = %.4f" % err)
    if err > args.tol:
        raise SystemExit("FAIL: gains not recovered (tol %.3f)" % args.tol)
    print("PASS")


if __name__ == "__main__":
    main()
