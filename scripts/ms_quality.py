#!/usr/bin/env python3
"""Quality metrics for a calibrated Measurement Set and its KillMS solutions.

The primary output of ddf-pipeline is the calibrated Measurement Set. This script
computes the cheap, self-consistency metrics on it (and on the killMS solution .npz)
that the smoke test deliberately does not: the normalised calibration residual, flag
fraction, weight statistics, model-vs-data correlation, and per-antenna solution
amplitude/phase/coverage.

Usage:
  ms_quality.py MS_OR_DIR [--colname COL] [--model-col COL] [--sols GLOB]
                 [--out PATH] [--tolerance PATH]

  MS_OR_DIR    a single .ms, or a run directory (its top-level *.ms are processed).
  --colname    calibrated column (default: auto-detect DATA_DI_CORRECTED, then
               CORRECTED_DATA, DATA_SUB, SCALED_DATA). Must exist, else exit 1.
  --model-col  model column (default: DD_PREDICT; falls back to MODEL_DATA/PREDICT_DATA).
  --sols       explicit solution file or glob (default: <ms>/killMS.*.sols.npz).
  --out        write JSON here (default: stdout).
  --tolerance  optional tolerance map (same shape as tests/data/smoke_reference.json).

Exit codes: 0 = metrics computed (and within tolerance if --tolerance given);
            1 = a required column is missing, or a metric is outside tolerance;
            2 = usage / IO error (unopenable MS, corrupt .npz).

This script is a plumbing check: it proves the arithmetic and the column/key
conventions, not the scientific correctness of the run (see docs/quality-control.md).
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

# ---------------------------------------------------------------------------
# Pure-numpy metric functions (unit-tested without python-casacore).
# ---------------------------------------------------------------------------


def _as_complex(a):
    return np.asarray(a, dtype=np.complex128)


def _flagged_mask(flag, shape):
    """Return a bool mask broadcast to `shape`, or an all-False mask."""
    if flag is None:
        return np.zeros(shape, dtype=bool)
    flag = np.asarray(flag, dtype=bool)
    return flag if flag.shape == shape else np.broadcast_to(flag, shape)


def residual_rms(data, model, flag=None):
    """RMS of (data - model) per correlation, and its ratio to the data RMS.

    data, model: complex (nrow, nchan, ncorr). flag: bool same shape, or None.
    Flagged visibilities are *excluded* from the mean (not zeroed, which would
    bias the RMS downward).
    """
    data = _as_complex(data)
    model = _as_complex(model)
    resid = data - model
    mask = _flagged_mask(flag, resid.shape)
    ncorr = resid.shape[-1]
    per_corr = []
    rms_over_data = []
    for c in range(ncorr):
        r = resid[..., c].ravel()
        d = data[..., c].ravel()
        good = (~mask[..., c].ravel()) & np.isfinite(r) & np.isfinite(d)
        r = r[good]
        d = d[good]
        rms = float(np.sqrt(np.mean(np.abs(r) ** 2))) if r.size else float("nan")
        drms = float(np.sqrt(np.mean(np.abs(d) ** 2))) if d.size else 0.0
        per_corr.append(rms)
        rms_over_data.append(rms / drms if drms > 0 else float("nan"))
    return {"rms_abs": per_corr, "rms_over_data": rms_over_data}


def _broadcast_sqrt_weight(weight, shape):
    """Return sqrt(weight) broadcast to `shape` (nrow, nchan, ncorr)."""
    w = np.asarray(weight, dtype=np.float64)
    if w.ndim == 1:
        w = w.reshape((-1, 1, 1))
    elif w.ndim == 2:
        w = w.reshape(w.shape + (1,))
    return np.broadcast_to(np.sqrt(w), shape)


def normalised_residual(data, model, weight, flag=None):
    """std((data - model) * sqrt(weight)) per correlation.

    A correct calibration and model leave this ~ 1 (Re/Im ~ 1/sqrt(2) each).
    Flagged visibilities are excluded from the standard deviation.
    """
    data = _as_complex(data)
    model = _as_complex(model)
    resid = data - model
    mask = _flagged_mask(flag, resid.shape)
    sw = _broadcast_sqrt_weight(weight, resid.shape)
    ncorr = resid.shape[-1]
    per_corr = []
    for c in range(ncorr):
        v = (resid[..., c] * sw[..., c]).ravel()
        v = v[~mask[..., c].ravel()]
        v = v[np.isfinite(v)]
        per_corr.append(float(np.std(v)) if v.size else float("nan"))
    return per_corr


def flag_fraction(flag):
    """Flag fraction: total, per correlation, per-channel median and worst channel."""
    flag = np.asarray(flag, dtype=bool)
    total = float(flag.mean())
    per_corr = [float(flag[..., c].mean()) for c in range(flag.shape[-1])]
    per_channel = flag.mean(axis=(0, 2)) if flag.ndim == 3 else flag.mean(axis=0)
    return {
        "total": total,
        "per_corr": per_corr,
        "per_channel_median": float(np.median(per_channel)),
        "worst_channel": float(np.max(per_channel)),
    }


def weight_stats(w):
    """Mean/min/max and non-positive/non-finite fraction of a weight array."""
    w = np.asarray(w, dtype=np.float64).ravel()
    finite = w[np.isfinite(w)]
    return {
        "mean": float(finite.mean()) if finite.size else float("nan"),
        "min": float(finite.min()) if finite.size else float("nan"),
        "max": float(finite.max()) if finite.size else float("nan"),
        "frac_nonpositive": float((finite <= 0).mean()) if finite.size else float("nan"),
        "frac_nonfinite": float((~np.isfinite(w)).mean()),
    }


def solution_stats(G, t0=None, t1=None):
    """Per-antenna/direction gain statistics from a killMS Sols.G array.

    G: complex (nt, nch, na, nd) scalar gains, or (nt, nch, na, nd, 2, 2) Jones
    matrices. The gain amplitude is the Frobenius norm / sqrt(2), so a scalar
    gain g (G = g.I) maps to |g|. Phase is read from the (0,0) element.
    """
    G = _as_complex(G)
    if G.ndim == 6:
        amp = np.sqrt(np.sum(np.abs(G) ** 2, axis=(-2, -1))) / np.sqrt(2.0)
        phi = np.angle(G[..., 0, 0])
    elif G.ndim == 4:
        amp = np.abs(G)
        phi = np.angle(G)
    else:
        raise ValueError(f"unexpected Sols.G shape {G.shape}")
    finite = np.isfinite(amp) & (amp != 0)
    amp_f = amp[finite]
    phi_f = phi[finite]
    # circular standard deviation of phase (wrap-free)
    if phi_f.size:
        r = abs(np.mean(np.exp(1j * phi_f)))
        circ_std = np.sqrt(-2.0 * np.log(r)) if r > 0 else float("inf")
    else:
        circ_std = float("nan")
    out = {
        "shape": list(G.shape),
        "amp_median": float(np.median(amp_f)) if amp_f.size else float("nan"),
        "amp_rms": float(np.std(amp_f)) if amp_f.size else float("nan"),
        "phase_rms_deg": float(np.degrees(circ_std)),
        "frac_nonfinite": float((~finite).mean()),
        "coverage": float(finite.mean()),
    }
    if t0 is not None and t1 is not None:
        out["t_range_s"] = [float(np.min(t0)), float(np.max(t1))]
    return out


def model_data_correlation(data, model, flag=None):
    """|sum(D . conj(M))| / (||D|| . ||M||) per correlation, a sanity check."""
    data = _as_complex(data)
    model = _as_complex(model)
    if flag is not None:
        flag = np.asarray(flag, dtype=bool)
        if flag.shape != data.shape:
            flag = np.broadcast_to(flag, data.shape)
        data = np.where(flag, 0.0 + 0.0j, data)
        model = np.where(flag, 0.0 + 0.0j, model)
    per_corr = []
    for c in range(data.shape[-1]):
        d = data[..., c].ravel()
        m = model[..., c].ravel()
        d = d[np.isfinite(d)]
        m = m[np.isfinite(m)]
        num = abs(np.vdot(d, m))
        den = np.linalg.norm(d) * np.linalg.norm(m)
        per_corr.append(float(num / den) if den > 0 else float("nan"))
    return per_corr


# ---------------------------------------------------------------------------
# python-casacore IO (pyrap.tables is the legacy alias).
# ---------------------------------------------------------------------------


def _tables():
    try:
        import casacore.tables as pt  # noqa: F401
    except ImportError:
        import pyrap.tables as pt  # noqa: F401
    return pt


def _colnames(table):
    try:
        return list(table.colnames())
    except AttributeError:
        return [c for c in table.keys()]


def _load_solution(filepath):
    """Load a killMS .npz and return (G, t0, t1) or None if schema unsupported."""
    data = dict(np.load(filepath, allow_pickle=False))
    if "Sols" not in data:
        return None
    sols = data["Sols"]
    G = np.asarray(sols["G"]) if hasattr(sols, "dtype") and "G" in sols.dtype.names else None
    if G is None:
        # Some files store the gains at top level.
        G = data.get("Gains", data.get("G"))
    t0 = np.asarray(sols["t0"]) if "t0" in (sols.dtype.names if hasattr(sols, "dtype") else {}) else None
    t1 = np.asarray(sols["t1"]) if "t1" in (sols.dtype.names if hasattr(sols, "dtype") else {}) else None
    return G, t0, t1


def _find_solutions(ms_path, explicit):
    if explicit:
        return sorted(glob.glob(explicit))
    patterns = [
        os.path.join(ms_path, "killMS.*.sols.npz"),
        os.path.join(ms_path, "**", "killMS.*.sols.npz"),
    ]
    found = []
    for pat in patterns:
        found.extend(glob.glob(pat, recursive=True))
    return sorted(set(found))


def _detect_colname(table, requested):
    if requested:
        return requested
    for candidate in ("DATA_DI_CORRECTED", "CORRECTED_DATA", "DATA_SUB", "SCALED_DATA"):
        if candidate in _colnames(table):
            return candidate
    return None


def _detect_model_col(table, requested):
    if requested:
        return requested
    for candidate in ("DD_PREDICT", "MODEL_DATA", "PREDICT_DATA"):
        if candidate in _colnames(table):
            return candidate
    return None


def _read_col(table, name, fallback=None):
    cols = _colnames(table)
    if name not in cols:
        return fallback
    return np.asarray(table.getcol(name))


def _ms_metrics(ms_path, colname_arg, model_arg, sols_arg):
    """Compute the metrics for one MS. Returns a dict, or raises ValueError/OSError."""
    pt = _tables()
    table = pt.table(ms_path, readonly=True, ack=False)
    try:
        cols = _colnames(table)
        colname = _detect_colname(table, colname_arg)
        if colname is None:
            raise ValueError(f"no calibrated column found in {ms_path} (looked for "
                             "DATA_DI_CORRECTED/CORRECTED_DATA/DATA_SUB/SCALED_DATA); "
                             "use --colname")
        model = _detect_model_col(table, model_arg)
        if model is None:
            raise ValueError(f"no model column found in {ms_path} (looked for "
                             "DD_PREDICT/MODEL_DATA/PREDICT_DATA); use --model-col")

        data = _read_col(table, "DATA")
        flag = _read_col(table, "FLAG")
        corrected = _read_col(table, colname)
        predict = _read_col(table, model)
        weight = _read_col(table, "WEIGHT_SPECTRUM", _read_col(table, "WEIGHT"))
        imaging_weight = _read_col(table, "IMAGING_WEIGHT")

        out = {
            "ms": ms_path,
            "columns_present": cols,
            "colname": colname,
            "model_col": model,
            "nrow": int(table.nrows()),
        }
        if data is not None:
            out["nchan"] = int(data.shape[1]) if data.ndim > 1 else 1
            out["ncorr"] = int(data.shape[-1]) if data.ndim > 1 else 1

        if corrected is not None and predict is not None:
            rr = residual_rms(corrected, predict, flag)
            out["residual"] = {
                "rms_abs": rr["rms_abs"],
                "rms_over_data": rr["rms_over_data"],
            }
            if data is not None:
                out["residual"]["rms_over_data"] = residual_rms(
                    corrected, predict, flag)["rms_over_data"]
            if weight is not None:
                out["residual"]["normalised_std"] = normalised_residual(
                    corrected, predict, weight, flag)
            out["model_data_corr"] = model_data_correlation(corrected, predict, flag)
        else:
            out["residual"] = None
            out["model_data_corr"] = None

        if flag is not None:
            out["flag_fraction"] = flag_fraction(flag)

        if imaging_weight is not None:
            out["imaging_weight"] = weight_stats(imaging_weight)
        if weight is not None:
            out["weight_spectrum"] = weight_stats(weight)

        out["solutions"] = []
        for solfile in _find_solutions(ms_path, sols_arg):
            entry = {"file": solfile}
            try:
                loaded = _load_solution(solfile)
            except Exception as exc:  # corrupt npz -> IO error, not a silent skip
                raise OSError(f"could not read {solfile}: {exc}")
            if loaded is None:
                entry["error"] = "unsupported schema (no Sols key)"
            else:
                G, t0, t1 = loaded
                if G is None:
                    entry["error"] = "unsupported schema (no gain array)"
                else:
                    entry.update(solution_stats(G, t0, t1))
            out["solutions"].append(entry)
        return out
    finally:
        try:
            table.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Tolerance application.
# ---------------------------------------------------------------------------


def _check_tolerance(metrics, tol):
    """Return a list of failure strings (empty = pass)."""
    failures = []
    r = metrics.get("residual") or {}
    if r and "rms_over_data" in r and "residual_rms_over_data_max" in tol:
        v = max([x for x in r["rms_over_data"] if np.isfinite(x)] or [float("nan")])
        if v > tol["residual_rms_over_data_max"]:
            failures.append(f"residual rms_over_data max {v:.3f} > {tol['residual_rms_over_data_max']}")
    for sol in metrics.get("solutions", []):
        if "error" in sol:
            continue
        if "amp_median" in sol and "sols_amp_median" in tol:
            lo, hi = tol["sols_amp_median"]
            if not (lo <= sol["amp_median"] <= hi):
                failures.append(f"sols amp_median {sol['amp_median']:.3f} not in [{lo},{hi}]")
        if "frac_nonfinite" in sol and "sols_frac_nonfinite_max" in tol:
            if sol["frac_nonfinite"] > tol["sols_frac_nonfinite_max"]:
                failures.append(f"sols frac_nonfinite {sol['frac_nonfinite']:.3%} > "
                                f"{tol['sols_frac_nonfinite_max']:.3%}")
        if "coverage" in sol and "sols_coverage_min" in tol:
            if sol["coverage"] < tol["sols_coverage_min"]:
                failures.append(f"sols coverage {sol['coverage']:.3%} < {tol['sols_coverage_min']:.3%}")
    ff = metrics.get("flag_fraction")
    if ff and "flag_worst_channel_max" in tol:
        if ff["worst_channel"] > tol["flag_worst_channel_max"]:
            failures.append(f"flag worst_channel {ff['worst_channel']:.3%} > "
                            f"{tol['flag_worst_channel_max']:.3%}")
    return failures


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ms_or_dir")
    parser.add_argument("--colname", default=None)
    parser.add_argument("--model-col", default=None)
    parser.add_argument("--sols", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--tolerance", default=None)
    args = parser.parse_args(argv[1:])

    if os.path.isdir(args.ms_or_dir) and args.ms_or_dir.endswith(".ms"):
        ms_paths = [args.ms_or_dir]
    elif os.path.isdir(args.ms_or_dir):
        ms_paths = sorted(glob.glob(os.path.join(args.ms_or_dir, "*.ms")))
        if not ms_paths:
            print(f"error: no *.ms in {args.ms_or_dir}", file=sys.stderr)
            return 2
    elif os.path.exists(args.ms_or_dir):
        ms_paths = [args.ms_or_dir]
    else:
        print(f"error: {args.ms_or_dir} does not exist", file=sys.stderr)
        return 2

    tol = None
    if args.tolerance:
        with open(args.tolerance) as fh:
            tol = json.load(fh)

    report = {"metrics": []}
    had_error = False
    for ms in ms_paths:
        try:
            report["metrics"].append(
                _ms_metrics(ms, args.colname, args.model_col, args.sols))
        except (ValueError, OSError) as exc:
            report["metrics"].append({"ms": ms, "error": str(exc)})
            had_error = True

    if tol:
        for m in report["metrics"]:
            if "error" in m:
                continue
            fails = _check_tolerance(m, tol)
            m["tolerance_failures"] = fails
            if fails:
                had_error = True

    payload = json.dumps(report, indent=2)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(payload + "\n")
    else:
        print(payload)

    return 1 if had_error else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
