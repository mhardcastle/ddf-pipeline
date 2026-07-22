#!/usr/bin/env python3
"""
generate_ms_set_simms.py

Drives simms (https://github.com/ratt-ru/simms) to produce a set of empty
Measurement Sets (MSs) tiling a frequency band [f0, f1] into nMS sub-bands,
each with the same channel width df, the same total time coverage dT split
into time bins of dt.

Optional diagnostic plots (--plot / --plot-save):
  - Elevation vs time
  - Azimuth vs time
  - UV coverage (km)
  - W-term vs time (km)

No astropy dependency: all astrometry uses pure numpy.

Example
-------
python generate_ms_set_simms.py \
    --nMS 2 --f0 100e6 --f1 101e6 --df 1e5 \
    --dt 8 --dT 21600 \
    --telescope ska1_10stations.tm \
    --ra-deg 170.0 --dec-deg -10.0 \
    --date "UTC,2025/3/31" --start-utc "08:00:00" \
    --outdir ms_output --ms-list ms_list.txt \
    --plot --plot-save obs_summary.png
"""

import argparse
import os
import sys
import tempfile
import numpy as np


# ---------------------------------------------------------------------------
# Frequency / time utilities
# ---------------------------------------------------------------------------

def plan_frequency_chunks(f0, f1, df, nMS):
    """
    Split [f0, f1] into nMS equal chunks each with the same integer number of
    channels.

    Priority: keep f0, f1, and nMS exactly as specified.
    Adapt:    df is adjusted slightly so that nMS * channels_per_ms channels
              span exactly (f1 - f0).  channels_per_ms is chosen as the nearest
              integer to (f1-f0)/df/nMS.

    This guarantees:
      - Every MS has the same number of channels (channels_per_ms).
      - The global channel grid covers [f0, f1] exactly.
      - df may differ from the requested value by at most half a channel width
        divided by nMS (typically sub-Hz to a few hundred Hz).
    """
    if f1 <= f0:
        raise ValueError("f1 must be greater than f0")
    if df <= 0:
        raise ValueError("df must be positive")
    if nMS <= 0:
        raise ValueError("nMS must be a positive integer")

    total_bw         = f1 - f0
    channels_per_ms  = max(1, int(round(total_bw / df / nMS)))
    total_channels   = nMS * channels_per_ms
    df_adapted       = total_bw / total_channels   # exact: covers f1-f0 in total_channels steps

    if abs(df_adapted - df) > df:
        raise ValueError(
            f"Adapted df ({df_adapted/1e3:.3f} kHz) deviates by more than 100% "
            f"from requested df ({df/1e3:.3f} kHz). Check f0/f1/df/nMS."
        )
    if abs(df_adapted - df) > 1.0:
        print(f"  Note: df adjusted from {df:.4f} Hz to {df_adapted:.4f} Hz "
              f"({df_adapted-df:+.4f} Hz, {(df_adapted-df)/df*100:+.4f}%) "
              f"so that {nMS} MSs × {channels_per_ms} channels cover "
              f"{total_bw/1e6:.6f} MHz exactly.")

    # Use df_adapted for all chunk boundaries (single multiply — no float drift)
    return [(f0 + i * channels_per_ms * df_adapted, channels_per_ms)
            for i in range(nMS)], df_adapted


def resolve_task_id_and_count(args):
    task_id = args.task_id
    if task_id is None:
        task_id = int(os.environ.get("SLURM_PROCID", 0))
    num_tasks = args.num_tasks
    if num_tasks is None:
        num_tasks = int(os.environ.get("SLURM_NTASKS", 1))
    if not (0 <= task_id < num_tasks):
        raise ValueError(f"task_id ({task_id}) must be in [0, num_tasks={num_tasks})")
    return task_id, num_tasks


def my_share(items, task_id, num_tasks):
    n = len(items)
    base, rem = divmod(n, num_tasks)
    start = task_id * base + min(task_id, rem)
    count = base + (1 if task_id < rem else 0)
    return items[start:start + count], start


# ---------------------------------------------------------------------------
# OSKAR .tm directory parser
# ---------------------------------------------------------------------------

def is_oskar_tm(path):
    return (
        os.path.isdir(path)
        and os.path.isfile(os.path.join(path, "position.txt"))
        and os.path.isfile(os.path.join(path, "layout.txt"))
    )


def read_oskar_tm(tm_dir):
    pos_file    = os.path.join(tm_dir, "position.txt")
    layout_file = os.path.join(tm_dir, "layout.txt")
    with open(pos_file) as f:
        parts = f.read().split()
    lon_deg = float(parts[0])
    lat_deg = float(parts[1])
    elev_m  = float(parts[2]) if len(parts) > 2 else 0.0
    stations = []
    with open(layout_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            e = float(cols[0])
            n = float(cols[1])
            u = float(cols[2]) if len(cols) > 2 else 0.0
            stations.append((e, n, u))
    return lon_deg, lat_deg, elev_m, stations


def write_simms_ascii(stations, path, dish_diam=87.0, mount="ALT-AZ", prefix="SKA"):
    with open(path, "w") as f:
        f.write("# E N U dish_diam station mount\n")
        for i, (e, n, u) in enumerate(stations):
            f.write(f"{e:.6f} {n:.6f} {u:.6f} {dish_diam:.1f} {prefix}{i:03d} {mount}\n")


# ---------------------------------------------------------------------------
# Astrometry (pure numpy, no astropy)
# ---------------------------------------------------------------------------

def ymdhms_to_jd(year, month, day, hour=0, minute=0, second=0):
    """Convert a calendar date/time to Julian Date."""
    a = (14 - month) // 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    jdn = day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
    return jdn + (hour - 12) / 24.0 + minute / 1440.0 + second / 86400.0


def jd_to_casa_date(jd):
    """
    Convert a Julian Date to a CASA date string 'YYYY/M/D/HH:MM:SS.ss'.
    Used to encode the observation centre time for simms.
    """
    # Meeus algorithm: JD → Gregorian calendar
    jd_shift = jd + 0.5
    z = int(jd_shift)
    f = jd_shift - z
    if z < 2299161:
        a = z
    else:
        alpha = int((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - alpha // 4
    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)
    day   = b - d - int(30.6001 * e)
    month = e - 1 if e < 14 else e - 13
    year  = c - 4716 if month > 2 else c - 4715
    # Fractional day → HH:MM:SS
    total_sec = f * 86400.0
    hh = int(total_sec // 3600);  total_sec -= hh * 3600
    mm = int(total_sec // 60);    ss = total_sec - mm * 60
    return f"{year}/{month}/{day}/{hh:02d}:{mm:02d}:{ss:05.2f}"


def parse_date_start(date_str, start_utc_str="00:00:00"):
    """
    Parse --date "UTC,YYYY/M/D" and --start-utc "HH:MM:SS" to a JD.
    Returns the Julian Date of the observation start.
    """
    date_part = date_str.split(",")[-1]          # "YYYY/M/D"
    yy, mm, dd = [int(x) for x in date_part.split("/")]
    hh, mn, ss = [int(x) for x in start_utc_str.split(":")]
    return ymdhms_to_jd(yy, mm, dd, hh, mn, ss)


def jd_to_gmst_rad(jd):
    """Greenwich Mean Sidereal Time in radians from Julian Date."""
    gmst_deg = 280.46061837 + 360.98564736629 * (jd - 2451545.0)
    return np.radians(gmst_deg % 360.0)


def compute_azel(jd_arr, ra_deg, dec_deg, lon_deg, lat_deg):
    """
    Return (elevation_deg, azimuth_deg, hour_angle_rad) arrays,
    one value per element of jd_arr.

    Azimuth is measured North-through-East (0=N, 90=E).
    """
    ra_rad  = np.radians(ra_deg)
    dec_rad = np.radians(dec_deg)
    lat_rad = np.radians(lat_deg)
    lon_rad = np.radians(lon_deg)

    gmst = jd_to_gmst_rad(jd_arr)
    lst  = gmst + lon_rad
    H    = lst - ra_rad                           # hour angle

    sin_el = (np.sin(lat_rad) * np.sin(dec_rad)
              + np.cos(lat_rad) * np.cos(dec_rad) * np.cos(H))
    el_deg = np.degrees(np.arcsin(np.clip(sin_el, -1.0, 1.0)))

    # Az: North = 0, East = 90 (standard astronomical convention)
    az_y = np.sin(H)
    az_x = np.cos(H) * np.sin(lat_rad) - np.tan(dec_rad) * np.cos(lat_rad)
    az_deg = (np.degrees(np.arctan2(az_y, az_x)) + 180.0) % 360.0

    return el_deg, az_deg, H


def compute_uvw(stations_enu, H_arr, dec_deg, lon_deg, lat_deg):
    """
    Compute UVW (metres) for every baseline and every time step.

    Parameters
    ----------
    stations_enu : (nStat, 3) array, ENU offsets in metres
    H_arr        : (nTime,) hour angles in radians
    dec_deg      : source declination in degrees
    lon_deg, lat_deg : array reference position in degrees

    Returns
    -------
    u, v, w : (nTime, nBase) arrays in metres
    """
    lon_rad = np.radians(lon_deg)
    lat_rad = np.radians(lat_deg)
    dec_rad = np.radians(dec_deg)

    # ENU → ITRF rotation (columns = East, North, Up unit vectors in ITRF)
    sl, cl = np.sin(lon_rad), np.cos(lon_rad)
    sf, cf = np.sin(lat_rad), np.cos(lat_rad)
    M_ENU_to_ITRF = np.array([
        [-sl,  -sf * cl,   cf * cl],
        [ cl,  -sf * sl,   cf * sl],
        [ 0.0,  cf,        sf     ],
    ])  # shape (3, 3), columns = [East, North, Up] in ITRF

    stat_itrf = stations_enu @ M_ENU_to_ITRF.T   # (nStat, 3)

    # All baselines i→j with i < j
    nStat = len(stat_itrf)
    bl = np.array([stat_itrf[i] - stat_itrf[j]
                   for i in range(nStat) for j in range(i + 1, nStat)])  # (nBase, 3)

    Bx, By, Bz = bl[:, 0], bl[:, 1], bl[:, 2]

    # UVW rotation (TMS convention):
    # u =  Bx sin(H) + By cos(H)
    # v = -Bx sin(d)cos(H) + By sin(d)sin(H) + Bz cos(d)
    # w =  Bx cos(d)cos(H) - By cos(d)sin(H) + Bz sin(d)
    sinH, cosH = np.sin(H_arr)[:, None], np.cos(H_arr)[:, None]   # (nTime, 1)
    sd, cd = np.sin(dec_rad), np.cos(dec_rad)

    u =  sinH * Bx + cosH * By
    v = -sd * cosH * Bx + sd * sinH * By + cd * Bz
    w =  cd * cosH * Bx - cd * sinH * By + sd * Bz

    return u, v, w   # (nTime, nBase) in metres


# ---------------------------------------------------------------------------
# Diagnostic plots
# ---------------------------------------------------------------------------

def make_observation_plots(
    stations_enu,
    ra_deg, dec_deg,
    lon_deg, lat_deg,
    jd_start, dT, dt,
    save_path=None,
):
    """
    Four-panel diagnostic figure:
      [0,0] Elevation (deg) vs time (h)
      [0,1] Azimuth   (deg) vs time (h)
      [1,0] UV coverage (km)
      [1,1] |W| (km) per baseline vs time (h)

    Parameters
    ----------
    stations_enu : (nStat, 3) ENU offsets in metres from read_oskar_tm
    ra_deg, dec_deg : source coordinates
    lon_deg, lat_deg : array reference position
    jd_start : Julian Date of observation start
    dT : total duration in seconds
    dt : time step in seconds
    save_path : if not None, save figure to this PNG path
    """
    import matplotlib
    matplotlib.use("Agg" if save_path and not _display_available() else "TkAgg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    # Time axis
    t_sec  = np.arange(0, dT, dt)
    t_hr   = t_sec / 3600.0
    jd_arr = jd_start + t_sec / 86400.0

    el_deg, az_deg, H_arr = compute_azel(jd_arr, ra_deg, dec_deg, lon_deg, lat_deg)

    stat_arr = np.array(stations_enu)           # (nStat, 3)
    u, v, w  = compute_uvw(stat_arr, H_arr, dec_deg, lon_deg, lat_deg)

    u_km = u / 1e3
    v_km = v / 1e3
    w_km = w / 1e3

    # ---- figure layout ----
    fig = plt.figure(figsize=(6.5, 5))
    fig.patch.set_facecolor("#0d1117")
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.48, wspace=0.38,
                            left=0.10, right=0.97, top=0.91, bottom=0.10)
    axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(2)]

    accent   = "#58a6ff"
    accent2  = "#f78166"
    gridcol  = "#21262d"
    textcol  = "#e6edf3"
    subfgcol = "#8b949e"

    for ax in axes:
        ax.set_facecolor("#161b22")
        ax.tick_params(colors=textcol, which="both", labelsize=4.5)
        ax.xaxis.label.set_color(textcol)
        ax.yaxis.label.set_color(textcol)
        ax.title.set_color(textcol)
        for spine in ax.spines.values():
            spine.set_edgecolor(gridcol)
        ax.grid(True, color=gridcol, linewidth=0.4)

    def style_ax(ax, title, xlabel, ylabel):
        ax.set_title(title, fontsize=5.5, fontweight="bold", pad=3.5)
        ax.set_xlabel(xlabel, fontsize=4.5)
        ax.set_ylabel(ylabel, fontsize=4.5)

    # ── [0] Elevation ──────────────────────────────────────────────────────
    ax = axes[0]
    ax.plot(t_hr, el_deg, color=accent, linewidth=0.8)
    ax.axhline(0, color=accent2, linewidth=0.4, linestyle="--", label="Horizon")
    ax.fill_between(t_hr, el_deg, 0, where=(el_deg > 0),
                    alpha=0.15, color=accent)
    style_ax(ax, "Elevation", "Time from start (h)", "Elevation (°)")
    ax.set_ylim(-10, 90)
    el_max = el_deg.max()
    ax.annotate(f"Peak {el_max:.1f}°",
                xy=(t_hr[np.argmax(el_deg)], el_max),
                xytext=(5, -8), textcoords="offset points",
                color=subfgcol, fontsize=4,
                arrowprops=dict(arrowstyle="->", color=subfgcol, lw=0.4))

    # ── [1] Azimuth ────────────────────────────────────────────────────────
    ax = axes[1]
    ax.scatter(t_hr, az_deg, c=el_deg, cmap="plasma",
               s=0.8, vmin=0, vmax=90, rasterized=True)
    sm = plt.cm.ScalarMappable(cmap="plasma",
                               norm=plt.Normalize(vmin=0, vmax=90))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.01, fraction=0.046)
    cb.set_label("Elevation (°)", color=textcol, fontsize=4)
    cb.ax.yaxis.set_tick_params(color=textcol, labelsize=4)
    plt.setp(cb.ax.yaxis.get_ticklabels(), color=textcol, fontsize=4)
    style_ax(ax, "Azimuth (coloured by elevation)",
             "Time from start (h)", "Azimuth (°, N=0 E=90)")
    ax.set_ylim(0, 360)
    ax.set_yticks([0, 90, 180, 270, 360])

    # ── [2] UV coverage ────────────────────────────────────────────────────
    ax = axes[2]
    stride = max(1, len(t_hr) // 300)
    u_s, v_s = u_km[::stride].ravel(), v_km[::stride].ravel()
    ax.scatter( u_s,  v_s, s=0.2, color=accent,  alpha=0.5, rasterized=True)
    ax.scatter(-u_s, -v_s, s=0.2, color=accent2, alpha=0.3, rasterized=True,
               label="Conjugate")
    ax.set_aspect("equal", adjustable="datalim")
    style_ax(ax, "UV coverage", "u (km)", "v (km)")
    ax.axhline(0, color=gridcol, linewidth=0.3)
    ax.axvline(0, color=gridcol, linewidth=0.3)
    max_bl = np.sqrt(u_km**2 + v_km**2).max()
    ax.annotate(f"Max baseline {max_bl:.2f} km",
                xy=(0.03, 0.95), xycoords="axes fraction",
                color=subfgcol, fontsize=4)

    # ── [3] W-term vs time ─────────────────────────────────────────────────
    ax = axes[3]
    w_abs = np.abs(w_km)
    stride_t = max(1, len(t_hr) // 500)
    for b in range(w_abs.shape[1]):
        ax.plot(t_hr[::stride_t], w_abs[::stride_t, b],
                color=accent, alpha=0.15, linewidth=0.3)
    ax.plot(t_hr, w_abs.max(axis=1),  color=accent2, linewidth=0.7,
            label="|W| max")
    ax.plot(t_hr, w_abs.mean(axis=1), color=accent,  linewidth=0.7,
            linestyle="--", label="|W| mean")
    ax.legend(fontsize=4, facecolor="#161b22", edgecolor=gridcol,
              labelcolor=textcol, handlelength=1.5, borderpad=0.5)
    style_ax(ax, "W-term per baseline", "Time from start (h)", "|w| (km)")

    # ── suptitle ───────────────────────────────────────────────────────────
    nStat = len(stations_enu)
    nBase = nStat * (nStat - 1) // 2
    fig.suptitle(
        f"Observation summary   RA={ra_deg:.2f}°  Dec={dec_deg:.2f}°   "
        f"dT={dT/3600:.1f}h  dt={dt:.0f}s   "
        f"{nStat} stations / {nBase} baselines",
        fontsize=5.5, color=textcol, y=0.975,
    )

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        print(f"  Plot saved to {save_path}")

    if _display_available() and not save_path:
        plt.show()
    elif _display_available() and save_path:
        plt.show()

    plt.close(fig)


def _display_available():
    """True if a graphical display is reachable."""
    import subprocess
    if sys.platform == "darwin":
        return True
    display = os.environ.get("DISPLAY", "")
    return bool(display)


# ---------------------------------------------------------------------------
# Flag zeroing
# ---------------------------------------------------------------------------

def zero_flags(ms_path):
    """
    Explicitly set all FLAG and FLAG_ROW values to False in an MS.

    simms creates an empty MS whose FLAG column may contain uninitialised or
    non-zero values.  Downstream tools (everybeam, DDECal, etc.) treat any
    non-zero flag as "flagged data", so we zero them explicitly here.
    """
    from casatools import table as tb_tool
    tb = tb_tool()
    tb.open(ms_path, nomodify=False)
    try:
        nrows = tb.nrows()
        if nrows == 0:
            return
        # FLAG: shape (nCorr, nChan, nRows) — zero everything
        flag = tb.getcol("FLAG")
        if flag.any():
            tb.putcol("FLAG", np.zeros_like(flag, dtype=bool))
        # FLAG_ROW: shape (nRows,)
        flag_row = tb.getcol("FLAG_ROW")
        if flag_row.any():
            tb.putcol("FLAG_ROW", np.zeros_like(flag_row, dtype=bool))
    finally:
        tb.close()
    print(f"  Flags zeroed ({nrows} rows)")


# ---------------------------------------------------------------------------
# Post-creation metadata corrections
# ---------------------------------------------------------------------------

def fix_ms_metadata(ms_path):
    """
    Correct ANTENNA and FEED metadata written by simms to match OSKAR conventions
    expected by everybeam.

    Key corrections:
      ANTENNA/MOUNT        'ALT-AZ'     → 'FIXED'   (SKA-LOW stations are fixed on the ground;
                                                       everybeam uses this to decide whether to
                                                       apply parallactic angle rotation)
      ANTENNA/NAME         'ALT-AZ'     → 's{i:04d} (station{i:03d})'
      ANTENNA/STATION      'P'          → 's{i:04d}'
      FEED/RECEPTOR_ANGLE  [0, 0]       → [π/2, 0]  (X feed rotated 90° relative to Y;
                                                       matches OSKAR's convention for SKA-LOW
                                                       and affects beam polarisation orientation)
      FEED/POL_RESPONSE    identity     → zeros      (OSKAR writes zeros; field is unused by
                                                       everybeam for phased-array beams)
    """
    from casatools import table as tb_tool

    # --- ANTENNA table ---
    tb = tb_tool()
    tb.open(os.path.join(ms_path, "ANTENNA"), nomodify=False)
    n_ant = tb.nrows()
    tb.putcol("MOUNT",   ["FIXED"] * n_ant)
    tb.putcol("NAME",    [f"s{i:04d} (station{i:03d})" for i in range(n_ant)])
    tb.putcol("STATION", [f"s{i:04d}" for i in range(n_ant)])
    tb.putcol("TYPE",    [""] * n_ant)               # OSKAR leaves this empty
    tb.putcol("DISH_DIAMETER", np.ones(n_ant))       # OSKAR uses 1m placeholder for phased arrays
    tb.close()

    # --- FEED table ---
    tb.open(os.path.join(ms_path, "FEED"), nomodify=False)
    n_feed = tb.nrows()
    # RECEPTOR_ANGLE shape: (n_receptors, n_feed) for putcol → (2, n_feed)
    receptor_angle = np.zeros((2, n_feed), dtype=np.float64)
    receptor_angle[0, :] = np.pi / 2      # first receptor (X): 90°
    receptor_angle[1, :] = 0.0            # second receptor (Y): 0°
    tb.putcol("RECEPTOR_ANGLE", receptor_angle)
    # POL_RESPONSE: set to zeros like OSKAR ((2,2) per feed, shape (2,2,n_feed) for putcol)
    pol_response = np.zeros((2, 2, n_feed), dtype=complex)
    tb.putcol("POL_RESPONSE", pol_response)
    # Match OSKAR FEED conventions exactly
    tb.putcol("BEAM_ID",           np.zeros(n_feed, dtype=np.int32))    # 0, not -1
    tb.putcol("SPECTRAL_WINDOW_ID", np.zeros(n_feed, dtype=np.int32))   # 0, not -1
    tb.putcol("INTERVAL",          np.zeros(n_feed, dtype=np.float64))  # 0., not 1e30
    tb.close()

    # --- SPECTRAL_WINDOW: match OSKAR exactly ---
    # NET_SIDEBAND: OSKAR writes 0 (unspecified), simms/CASA writes 1 (upper sideband).
    # everybeam reads CHAN_FREQ directly, but match OSKAR to be safe.
    # NAME and FREQ_GROUP_NAME: OSKAR leaves blank, simms writes '00'/'Group 1'.
    tb.open(os.path.join(ms_path, "SPECTRAL_WINDOW"), nomodify=False)
    n_spw = tb.nrows()
    tb.putcol("NET_SIDEBAND",   np.zeros(n_spw, dtype=np.int32))
    tb.putcol("NAME",           [""] * n_spw)
    tb.putcol("FREQ_GROUP_NAME",[""] * n_spw)
    tb.close()

    # --- FIELD: clear NAME (OSKAR leaves it blank) ---
    tb.open(os.path.join(ms_path, "FIELD"), nomodify=False)
    tb.putcol("NAME", [""] * tb.nrows())
    tb.close()

    # --- SOURCE: delete the subtable entirely (OSKAR doesn't write one) ---
    # simms always calls sm.setsource() which writes a SOURCE subtable;
    # everybeam checks nrows()>0 to decide whether to use it.
    # NOTE: we use shutil (not tabledelete from python-casacore) because mixing
    # casatools and python-casacore in the same process causes C-level heap
    # corruption (they each link their own casacore C++ runtime).
    import shutil as _shutil
    source_path = os.path.join(ms_path, "SOURCE")
    if os.path.exists(source_path):
        # 1. Unlink keyword so the main table no longer references SOURCE
        tb.open(ms_path, nomodify=False)
        if "SOURCE" in tb.keywordnames():
            tb.removekeyword("SOURCE")
        tb.close()
        # 2. Delete the subtable directory from disk
        _shutil.rmtree(source_path)

    # --- OBSERVATION: clear CASA-specific fields to match OSKAR ---
    tb.open(os.path.join(ms_path, "OBSERVATION"), nomodify=False)
    tb.putcol("RELEASE_DATE", np.zeros(1, dtype=np.float64))
    tb.putcol("OBSERVER", [""])
    tb.putcol("PROJECT",  [""])
    tb.close()

    # --- STATE: fix OBS_MODE separator ('.' → '#') ---
    # simms writes 'OBSERVE_TARGET.ON_SOURCE'; OSKAR writes 'OBSERVE_TARGET#ON_SOURCE'
    tb.open(os.path.join(ms_path, "STATE"), nomodify=False)
    obs_modes = tb.getcol("OBS_MODE")
    tb.putcol("OBS_MODE", [m.replace(".", "#") for m in obs_modes])
    tb.close()

    # --- POINTING table: empty it preserving the valid MSPointing structure ---
    # Problem: a generic tb.create() produces a plain table that fails
    # casacore's MSPointing(const Table&) validator because Measure columns
    # (DIRECTION, TARGET, TIME) are missing their MEASINFO keywords.
    # Solution: copy the original simms POINTING (which IS a valid MSPointing)
    # with norows=True — this preserves all column definitions and MEASINFO
    # keywords, giving a valid empty subtable identical in structure to OSKAR's.
    import shutil as _shutil

    pointing_path    = os.path.join(ms_path, "POINTING")
    pointing_tmp     = pointing_path + "._empty"

    # 1. Copy structure only (0 rows) from the valid simms-created POINTING
    pt = tb_tool()
    pt.open(pointing_path)
    pt.copy(pointing_tmp, norows=True, returnobject=False)
    pt.close()

    # 2. Replace the original with the empty copy
    _shutil.rmtree(pointing_path)
    os.rename(pointing_tmp, pointing_path)

    print(f"  Metadata fixed: MOUNT=FIXED, NAME/STATION corrected, "
          f"RECEPTOR_ANGLE=[π/2,0], POINTING cleared ({n_ant} antennas)")


# ---------------------------------------------------------------------------
# Checkpoint / done-file helpers
# ---------------------------------------------------------------------------

def _done_file(ms_path):
    """Path of the sentinel file that marks a completed MS chunk."""
    return ms_path + ".done"


def _is_done(ms_path):
    return os.path.exists(_done_file(ms_path))


def _mark_done(ms_path):
    """Touch <ms_path>.done to record successful completion."""
    import datetime
    ts = datetime.datetime.utcnow().isoformat() + "Z"
    with open(_done_file(ms_path), "w") as fh:
        fh.write("completed " + ts + "\n")


# ---------------------------------------------------------------------------
# Optional column removal
# ---------------------------------------------------------------------------

def _remove_optional_cols(ms_path, cols):
    """
    Remove columns from the main MS table if they exist.
    Silently skips columns that are not present.
    Useful for stripping MODEL_DATA / CORRECTED_DATA that simms may add,
    before using the MS as a copy template (saves disk space on all copies).
    """
    from casatools import table as tb_tool
    tb = tb_tool()
    tb.open(ms_path, nomodify=False)
    existing = tb.colnames()
    to_remove = [c for c in cols if c in existing]
    if to_remove:
        tb.removecols(to_remove)
        print(f"  Removed columns from main table: {to_remove}")
    tb.close()


# ---------------------------------------------------------------------------
# SPW frequency update (used when copying a template MS)
# ---------------------------------------------------------------------------

def update_spw_frequencies(ms_path, start_freq_hz, num_channels, df_hz):
    """
    Update SPECTRAL_WINDOW/CHAN_FREQ and REF_FREQUENCY for a copied MS.

    Only these two columns carry the per-chunk frequency information;
    everything else (CHAN_WIDTH, EFFECTIVE_BW, RESOLUTION, TOTAL_BANDWIDTH,
    NUM_CHAN) is identical across chunks and is already correct in the copy.

    Parameters
    ----------
    ms_path      : path to the copied MS
    start_freq_hz: centre frequency of the first channel in this chunk (Hz)
    num_channels : number of channels (same for every chunk)
    df_hz        : channel width (Hz)
    """
    from casatools import table as tb_tool
    import numpy as np

    chan_freqs = np.array(
        [start_freq_hz + k * df_hz for k in range(num_channels)],
        dtype=np.float64
    )
    ref_freq = chan_freqs[0]   # same convention as simms

    tb = tb_tool()
    tb.open(os.path.join(ms_path, "SPECTRAL_WINDOW"), nomodify=False)
    # CHAN_FREQ shape per row: (nchan,) — putcol expects (nchan, nSPW)
    tb.putcol("CHAN_FREQ",     chan_freqs[:, None])   # (nchan, 1)
    tb.putcol("REF_FREQUENCY", np.array([ref_freq]))
    tb.close()


# ---------------------------------------------------------------------------
# simms helpers
# ---------------------------------------------------------------------------

def hz_to_simms_str(freq_hz):
    return f"{freq_hz / 1e6:.9f}MHz"


def deg_to_direction(ra_deg, dec_deg, frame="J2000"):
    return f"{frame},{ra_deg}deg,{dec_deg}deg"


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate empty MSs with simms tiling a frequency band.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--nMS",  type=int,   required=True)
    p.add_argument("--f0",   type=float, required=True, help="Band start, Hz.")
    p.add_argument("--f1",   type=float, required=True, help="Band end, Hz.")
    p.add_argument("--df",   type=float, default=1e5,   help="Channel width, Hz.")
    p.add_argument("--dt",   type=float, default=None,
                   help="Dump (integration) time in seconds. "
                        "Provide either --dt or --num-time-steps, not both.")
    p.add_argument("--num-time-steps", type=int, default=None,
                   help="Number of time steps. dt is computed as dT/num-time-steps. "
                        "Alternative to --dt; matches OSKAR's num_time_steps convention.")
    p.add_argument("--dT",   type=float, required=True,
                   help="Total integration time per MS, seconds.")

    p.add_argument("--telescope", type=str, required=True,
                   help="OSKAR .tm dir, simms ASCII file, or CASA table path.")
    p.add_argument("--pos-type", type=str, default=None,
                   choices=["CASA", "ascii"],
                   help="Force position file type. Auto-detected when not set.")
    p.add_argument("--tel-name", type=str, default="OSKAR",
                   help="Telescope label.")
    p.add_argument("--lon-lat", type=str, default=None,
                   help="'lon_deg,lat_deg[,elev_m]'. Auto-read from .tm for OSKAR dirs.")
    p.add_argument("--dish-diam", type=float, default=87.0,
                   help="Station diameter in metres (for ASCII generation).")
    p.add_argument("--coords", type=str, default=None,
                   choices=["itrf", "enu", "wgs84"],
                   help="Coordinate system for ASCII file.")

    p.add_argument("--ra-deg",  type=float, default=170.0)
    p.add_argument("--dec-deg", type=float, default=-10.0)
    p.add_argument("--date",    type=str,   default="UTC,2025/3/31",
                   help="Observation date, CASA format 'UTC,YYYY/M/D'.")
    p.add_argument("--start-utc", type=str, default="00:00:00",
                   help="Observation start time within --date, as HH:MM:SS.")
    p.add_argument("--stokes",  type=str,   default="XX XY YX YY",
                   help="Stokes parameters. Use XX XY YX YY (linear, default) for SKA-LOW/everybeam.")
    p.add_argument("--feed",    type=str,   default="perfect X Y",
                   help="Feed type. Use perfect X Y (linear, default) for SKA-LOW/everybeam.")

    p.add_argument("--outdir",    type=str, default="ms_output")
    p.add_argument("--ms-prefix", type=str, default="chunk")
    p.add_argument("--ms-list",   type=str, default="ms_list.txt")

    p.add_argument("--task-id",   type=int, default=None,
                   help="0-based task index. Defaults to $SLURM_PROCID.")
    p.add_argument("--num-tasks", type=int, default=None,
                   help="Total parallel tasks. Defaults to $SLURM_NTASKS.")

    p.add_argument("--n-stations", type=int, default=None, metavar="N",
                   help="Randomly subsample N stations from the telescope model. "
                        "Default: use all stations.")
    p.add_argument("--station-seed", type=int, default=42,
                   help="Random seed for station subsampling. Default: 42.")

    p.add_argument("--plot",      action="store_true",
                   help="Show diagnostic observation plots.")
    p.add_argument("--plot-save", type=str, default=None, metavar="FILE.png",
                   help="Save the diagnostic plot to this PNG file "
                        "(implies --plot; does not require a display).")

    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    do_plot = args.plot or (args.plot_save is not None)

    # --- Time ---
    if args.dt is None and args.num_time_steps is None:
        raise ValueError("Provide either --dt or --num-time-steps.")
    if args.dt is not None and args.num_time_steps is not None:
        raise ValueError("Provide --dt or --num-time-steps, not both.")
    if args.num_time_steps is not None:
        if args.dT % args.num_time_steps != 0 and abs(args.dT / args.num_time_steps - round(args.dT / args.num_time_steps)) > 1e-6:
            raise ValueError(f"dT ({args.dT}) is not evenly divisible by --num-time-steps ({args.num_time_steps}).")
        args.dt = args.dT / args.num_time_steps
    num_time_steps_f = args.dT / args.dt
    if abs(num_time_steps_f - round(num_time_steps_f)) > 1e-6:
        raise ValueError(f"dT / dt = {num_time_steps_f} is not an integer.")
    num_time_steps  = int(round(num_time_steps_f))
    synthesis_hours = args.dT / 3600.0

    # simms uses usehourangle=True internally: the date parameter is the
    # MERIDIAN TRANSIT time (HA=0), and the synthesis window is placed
    # symmetrically as [transit - half_sidereal, transit + half_sidereal].
    #
    # Because simms counts synthesis in SIDEREAL hours but OSKAR uses SOLAR
    # hours, a sidereal/solar correction is needed to align start times:
    #   1 sidereal hour = 86164.0905/24 = 3590.171 solar seconds
    #   vs 1 solar hour = 3600 solar seconds
    #   correction per half-synthesis = synthesis_hours × (3600 - 3590.171) / 2
    #                                 ≈ 4.91 solar seconds per hour of synthesis
    # For 8h synthesis this is ~39.3 s — enough to shift the TIME column
    # by that amount compared to OSKAR if not corrected.
    # simms uses usehourangle=True internally: it places the synthesis window
    # as [ref - half_HA, ref + half_HA] where half_HA is in SIDEREAL seconds.
    #   1 sidereal hour = 86164.0905 / 24 = 3590.171 solar seconds  (not 3600)
    #
    # To get the solar observation START at the user-specified time:
    #   ref = obs_start + (synthesis_hours × sidereal_seconds_per_hour / 2)
    # This ensures simms_start = ref - half_sidereal = obs_start (in solar time).
    _SIDEREAL_S_PER_H = 86164.0905 / 24.0    # solar seconds per sidereal hour = 3590.171
    _half_sidereal_s  = synthesis_hours * _SIDEREAL_S_PER_H / 2.0

    jd_obs_start = parse_date_start(args.date, args.start_utc)
    jd_simms_ref = jd_obs_start + _half_sidereal_s / 86400.0
    simms_date   = f"UTC,{jd_to_casa_date(jd_simms_ref)}"
    print(f"  Simms date (HA=0 reference) : {simms_date}  "
          f"(half-synthesis sidereal: {_half_sidereal_s:.1f} s = {_half_sidereal_s/3600:.4f} h)")
    print(f"  Expected obs start (UTC)    : {args.date.split(',')[-1]} {args.start_utc}")

    # --- Antenna table: auto-detect OSKAR .tm ---
    tel_path      = args.telescope
    pos_type      = args.pos_type
    coords        = args.coords
    lon_lat       = args.lon_lat
    _tmp_ascii    = None
    stations_enu  = None   # (nStat, 3) for plotting
    lon_deg = lat_deg = None

    if is_oskar_tm(tel_path):
        print(f"Detected OSKAR .tm directory: {tel_path}")
        lon_deg, lat_deg, elev_m, stations = read_oskar_tm(tel_path)
        stations_enu = np.array(stations)
        print(f"  Reference position : lon={lon_deg}, lat={lat_deg}, elev={elev_m} m")
        print(f"  Stations read      : {len(stations)}")

        # --- Optional station subsampling ---
        selected_indices = list(range(len(stations)))
        if args.n_stations is not None:
            n_avail = len(stations)
            if args.n_stations > n_avail:
                raise ValueError(
                    f"--n-stations {args.n_stations} exceeds the number of "
                    f"available stations ({n_avail})."
                )
            rng = np.random.default_rng(args.station_seed)
            selected_indices = sorted(
                rng.choice(n_avail, size=args.n_stations, replace=False).tolist()
            )
            stations = [stations[i] for i in selected_indices]
            stations_enu = np.array(stations)
            print(f"  Subsampled {args.n_stations}/{n_avail} stations "
                  f"(seed={args.station_seed}): indices {selected_indices}")
        else:
            print(f"  Using all {len(stations)} stations")

        _tmp_ascii = tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", prefix="ska_antennas_", delete=False
        )
        write_simms_ascii(stations, _tmp_ascii.name,
                          dish_diam=args.dish_diam, prefix="SKA")
        _tmp_ascii.close()
        tel_path = _tmp_ascii.name

        pos_type = pos_type or "ascii"
        coords   = coords   or "enu"
        lon_lat  = lon_lat  or f"{lon_deg},{lat_deg},{elev_m}"
        print(f"  ASCII antenna file : {tel_path}")
        print(f"  lon-lat            : {lon_lat}")
    else:
        pos_type = pos_type or "CASA"
        coords   = coords   or "itrf"
        tel_path = os.path.abspath(tel_path)

    # If plotting but not an OSKAR tm, try to parse lon/lat from --lon-lat
    if do_plot and stations_enu is None and lon_lat is not None:
        ll = [float(x) for x in lon_lat.split(",")]
        lon_deg, lat_deg = ll[0], ll[1]

    # --- Diagnostic plots (before running simms, so visible on dry-run too) ---
    if do_plot:
        if stations_enu is None:
            print("Warning: --plot requires an OSKAR .tm dir or --lon-lat "
                  "with known station positions; skipping plot.")
        else:
            jd_start = parse_date_start(args.date, args.start_utc)
            print(f"Generating observation plots "
                  f"(JD start={jd_start:.5f}, start_utc={args.start_utc}) ...")
            make_observation_plots(
                stations_enu=stations_enu,
                ra_deg=args.ra_deg, dec_deg=args.dec_deg,
                lon_deg=lon_deg, lat_deg=lat_deg,
                jd_start=jd_start, dT=args.dT, dt=args.dt,
                save_path=args.plot_save,
            )

    # --- Frequency chunking ---
    chunks, df_adapted = plan_frequency_chunks(args.f0, args.f1, args.df, args.nMS)
    args.df = df_adapted   # use the (possibly adjusted) channel width throughout

    # --- Task slice ---
    task_id, num_tasks = resolve_task_id_and_count(args)
    my_chunks, offset  = my_share(chunks, task_id, num_tasks)

    direction = deg_to_direction(args.ra_deg, args.dec_deg)
    df_str    = hz_to_simms_str(args.df)

    print(f"Planned {args.nMS} MS(s) across {num_tasks} task(s); "
          f"task {task_id} handles {len(my_chunks)} MS(s).")
    print(f"  Band : {args.f0/1e6:.6f} - {args.f1/1e6:.6f} MHz, "
          f"df = {args.df/1e3:.3f} kHz, {chunks[0][1]} channels/MS")
    print(f"  Time : dT={args.dT} s ({synthesis_hours:.4f} h), "
          f"dt={args.dt} s, steps={num_time_steps}")
    print(f"  Direction : {direction}, date : {args.date} {args.start_utc}")
    print(f"  pos_type={pos_type}, coords={coords}, lon_lat={lon_lat}")

    if not args.dry_run:
        from simms.casasm import makems

    os.makedirs(args.outdir, exist_ok=True)
    ms_paths      = []
    template_path = None   # path of the first (template) MS, reused by later chunks

    try:
        for local_i, (start_freq, num_channels) in enumerate(my_chunks):
            i        = offset + local_i
            ms_name  = f"{args.ms_prefix}_{i:04d}.MS"
            ms_path  = os.path.join(args.outdir, ms_name)
            ms_paths.append(os.path.abspath(ms_path))

            freq0_str = hz_to_simms_str(start_freq)
            print(
                f"[task {task_id}] [{local_i+1}/{len(my_chunks)}] "
                f"(global {i+1}/{args.nMS}) {ms_path}: "
                f"freq0={freq0_str}, nchan={num_channels}, df={df_str}"
            )

            if args.dry_run:
                continue

            if _is_done(ms_path):
                print(f"  Already done (found {_done_file(ms_path)}), skipping.")
                # Still need template_path for subsequent chunks
                if local_i == 0:
                    template_path = ms_path
                continue

            if local_i == 0:
                # --- First chunk: run simms once, apply all metadata fixes ---
                makems(
                    msname=ms_path,
                    tel=args.tel_name,
                    pos=tel_path,
                    pos_type=pos_type,
                    coords=coords,
                    lon_lat=lon_lat,
                    direction=[direction],
                    synthesis=synthesis_hours,
                    dtime=int(args.dt),
                    freq0=[freq0_str],
                    dfreq=[df_str],
                    nchan=[num_channels],
                    nbands=1,
                    stokes=args.stokes,
                    feed=args.feed,
                    date=simms_date,
                    outdir=None,
                )
                zero_flags(ms_path)
                fix_ms_metadata(ms_path)
                # Remove MODEL_DATA and CORRECTED_DATA if simms created them —
                # they waste disk space and are not needed in the template copies.
                _remove_optional_cols(ms_path, ["MODEL_DATA", "CORRECTED_DATA"])

                # Add PHASED_ARRAY subtable (only needed on the template;
                # all copytree copies inherit it automatically).
                # Run as a subprocess to avoid casatools/python-casacore ABI clash.
                if _oskar_tm_dir is not None:
                    import subprocess as _sp
                    _script = os.path.join(
                        os.path.dirname(os.path.abspath(__file__)),
                        "add_phased_array.py")
                    print(f"  Adding PHASED_ARRAY table to template MS ...")
                    _sp.run(
                        [sys.executable, _script,
                         ms_path, "--tm", _oskar_tm_dir, "--overwrite"],
                        check=True)

                template_path = ms_path

            else:
                # --- Subsequent chunks: copy template + patch SPECTRAL_WINDOW only ---
                # UVW, FLAG, TIME, ANTENNA, FEED … are identical across chunks.
                # Only CHAN_FREQ and REF_FREQUENCY differ between frequency chunks.
                import shutil as _shutil
                if os.path.exists(ms_path):
                    _shutil.rmtree(ms_path)
                _shutil.copytree(template_path, ms_path)
                update_spw_frequencies(ms_path, start_freq, num_channels, args.df)
                print(f"  Copied from template, SPW updated to {freq0_str}")

            _mark_done(ms_path)

    finally:
        if _tmp_ascii is not None and os.path.exists(_tmp_ascii.name):
            os.unlink(_tmp_ascii.name)

    # --- Write per-task MS list ---
    if num_tasks > 1:
        root, ext = os.path.splitext(args.ms_list)
        list_path = f"{root}.task{task_id:04d}{ext}"
    else:
        list_path = args.ms_list

    with open(list_path, "w") as f:
        for p in ms_paths:
            f.write(p + "\n")

    print(f"\n[task {task_id}] Wrote {len(ms_paths)} MS path(s) to {list_path}")


if __name__ == "__main__":
    main()
