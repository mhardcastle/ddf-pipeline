#!/usr/bin/env python3
"""
add_phased_array.py

Adds a PHASED_ARRAY subtable to one or more Measurement Sets, reading element
positions from an OSKAR telescope model (.tm) directory.

Uses python-casacore (casacore.tables) throughout — NOT casatools — so that
the putcol convention matches what OSKAR writes:
    putcol(col, data (nRows, d1, d2, ...)) → cell for row s = data[s]

This avoids both the transpose issue (casatools transposes cells internally,
giving the wrong shape for non-square array columns) and the ABI crash caused
by mixing casatools and python-casacore in the same process.

ELEMENT_OFFSET shape per cell: (3, nElem)   → row = ITRF component, col = element
ELEMENT_FLAG   shape per cell: (2, nElem)   → row = polarisation
COORDINATE_AXES shape per cell: (3, 3)      → rows = [East, North, Up] in ITRF
"""

import argparse
import glob
import os
import sys
import numpy as np


# ---------------------------------------------------------------------------
# Astrometry helpers
# ---------------------------------------------------------------------------

def geodetic_to_ecef(lon_deg, lat_deg, alt_m=0.0):
    """Convert WGS84 geodetic to ITRF/ECEF XYZ (metres)."""
    a  = 6378137.0
    f  = 1.0 / 298.257223563
    e2 = 2.0 * f - f**2
    lon = np.radians(lon_deg)
    lat = np.radians(lat_deg)
    N = a / np.sqrt(1.0 - e2 * np.sin(lat)**2)
    x = (N + alt_m) * np.cos(lat) * np.cos(lon)
    y = (N + alt_m) * np.cos(lat) * np.sin(lon)
    z = (N * (1.0 - e2) + alt_m) * np.sin(lat)
    return np.array([x, y, z])


def itrf_to_lonlat(xyz):
    """
    ITRF XYZ → (lon_rad, geodetic_lat_rad) via WGS84 Bowring iteration.
    Uses geodetic (not geocentric) latitude — same as OSKAR.
    Geocentric latitude differs by ~9 arcmin at SKA-LOW, corrupting COORDINATE_AXES.
    """
    x, y, z = xyz
    lon = np.arctan2(y, x)
    a  = 6378137.0
    f  = 1.0 / 298.257223563
    e2 = 2.0 * f - f**2
    p  = np.sqrt(x**2 + y**2)
    lat = np.arctan2(z, p * (1.0 - e2))
    for _ in range(10):
        N = a / np.sqrt(1.0 - e2 * np.sin(lat)**2)
        lat_new = np.arctan2(z + e2 * N * np.sin(lat), p)
        if abs(lat_new - lat) < 1e-14:
            break
        lat = lat_new
    return lon, lat_new


def enu_to_itrf_matrix(xyz):
    """
    3×3 rotation matrix whose rows are [East, North, Up] unit vectors in ITRF.
    M @ v_ITRF = v_ENU  (projects ITRF vector onto local ENU axes)
    """
    lon, lat = itrf_to_lonlat(xyz)
    sl, cl = np.sin(lon), np.cos(lon)
    sf, cf = np.sin(lat), np.cos(lat)
    east  = np.array([-sl,         cl,         0.0])
    north = np.array([-sf * cl,   -sf * sl,    cf ])
    up    = np.array([ cf * cl,    cf * sl,    sf ])
    return np.array([east, north, up])


# ---------------------------------------------------------------------------
# Station position matching
# ---------------------------------------------------------------------------

def tm_station_itrf_positions(tm_dir):
    """
    Compute ITRF XYZ for every station in an OSKAR .tm directory.
    Reads position.txt (lon/lat/elev) and layout.txt (ENU offsets).
    Returns (positions (nStat, 3), station_dirs list).
    """
    with open(os.path.join(tm_dir, "position.txt")) as f:
        parts = f.read().split()
    lon_deg = float(parts[0])
    lat_deg = float(parts[1])
    elev_m  = float(parts[2]) if len(parts) > 2 else 0.0
    ref_itrf = geodetic_to_ecef(lon_deg, lat_deg, elev_m)

    lon_r = np.radians(lon_deg)
    lat_r = np.radians(lat_deg)
    sl, cl = np.sin(lon_r), np.cos(lon_r)
    sf, cf = np.sin(lat_r), np.cos(lat_r)
    M = np.array([[-sl, -sf*cl,  cf*cl],
                  [ cl, -sf*sl,  cf*sl],
                  [ 0.,  cf,     sf   ]])

    station_enu = []
    with open(os.path.join(tm_dir, "layout.txt")) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            station_enu.append([float(cols[0]), float(cols[1]),
                                 float(cols[2]) if len(cols) > 2 else 0.0])
    station_enu = np.array(station_enu)
    positions   = ref_itrf + station_enu @ M.T

    station_dirs = sorted(glob.glob(os.path.join(tm_dir, "station[0-9]*")))
    if len(station_dirs) != len(positions):
        raise ValueError(
            f"layout.txt has {len(positions)} stations but found "
            f"{len(station_dirs)} station subdirectories in {tm_dir!r}.")
    return positions, station_dirs


def match_ms_to_tm_stations(ms_positions, tm_dir):
    """
    Match each MS ANTENNA to its .tm station by nearest ITRF position.
    Injectivity check: no two MS stations may map to the same .tm station.
    """
    tm_positions, station_dirs = tm_station_itrf_positions(tm_dir)
    diff        = ms_positions[:, None, :] - tm_positions[None, :, :]
    dist_matrix = np.linalg.norm(diff, axis=2)
    matched_indices = dist_matrix.argmin(axis=1).tolist()
    match_dists     = dist_matrix[np.arange(len(ms_positions)), matched_indices]

    seen = {}
    for ms_i, tm_j in enumerate(matched_indices):
        if tm_j in seen:
            raise ValueError(
                f"Duplicate match: MS stations {seen[tm_j]} and {ms_i} "
                f"both map to .tm station {tm_j}. Wrong --tm?")
        seen[tm_j] = ms_i

    print(f"  Positional offset MS vs .tm ITRF: "
          f"max={match_dists.max():.2f} m  mean={match_dists.mean():.2f} m")
    return [station_dirs[j] for j in matched_indices], matched_indices


# ---------------------------------------------------------------------------
# Station element layout reader
# ---------------------------------------------------------------------------

def read_station_layouts(station_dirs):
    """
    Read per-station element ENU offsets from a list of station<NNN>/ dirs.
    Each layout.txt has E N [U] columns in metres; U defaults to 0 if absent.
    Returns list of np.ndarray, each shape (n_elements_i, 3).
    """
    layouts = []
    for sdir in station_dirs:
        layout_file = os.path.join(sdir, "layout.txt")
        if not os.path.isfile(layout_file):
            raise FileNotFoundError(f"Missing: {layout_file}")
        elements = []
        with open(layout_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                cols = line.split()
                elements.append([float(cols[0]), float(cols[1]),
                                  float(cols[2]) if len(cols) >= 3 else 0.0])
        layouts.append(np.array(elements, dtype=np.float64))
        print(f"  {os.path.basename(sdir)}: {len(elements)} elements")
    return layouts


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------

def add_phased_array_table(ms_path, tm_dir, overwrite=False,
                            oskar_element_order=False):
    """
    Create and populate a PHASED_ARRAY subtable inside an existing MS.

    Uses python-casacore (casacore.tables) — NOT casatools — so that
    putcol convention matches OSKAR:  data[s] is the cell for station s.

    Element layout convention
    -------------------------
    Default (everybeam):   per-cell shape (3, nElem) and (2, nElem)
                           row = ITRF component / polarisation, col = element
    --oskar-element-order: per-cell shape (nElem, 3) and (nElem, 2)  [legacy]
    """
    from casacore.tables import (
        table as pt, makescacoldesc, makearrcoldesc, maketabdesc
    )

    pa_path = os.path.join(ms_path, "PHASED_ARRAY")

    if os.path.exists(pa_path):
        if overwrite:
            import shutil
            shutil.rmtree(pa_path)
            print(f"  Removed existing PHASED_ARRAY at {pa_path}")
        else:
            print(f"  PHASED_ARRAY already exists in {ms_path} (use --overwrite).")
            return

    # --- Read station ITRF positions from ANTENNA subtable ---
    ant_tb   = pt(os.path.join(ms_path, "ANTENNA"))
    # casacore putcol/getcol convention: data[s] = cell for row s
    # ANTENNA/POSITION: getcol returns (nStat, 3)
    positions = ant_tb.getcol("POSITION").astype(np.float64)   # (nStat, 3)
    ant_tb.close()
    n_stations = positions.shape[0]
    print(f"  Stations in ANTENNA table : {n_stations}")

    # --- Match MS stations to .tm stations by position ---
    print(f"  Matching {n_stations} MS stations → .tm stations in {tm_dir}")
    matched_dirs, matched_indices = match_ms_to_tm_stations(positions, tm_dir)
    print(f"  Matched (first 5): {matched_indices[:5]} ...")

    # --- Read element layouts for matched stations ---
    print(f"  Reading element layouts")
    layouts    = read_station_layouts(matched_dirs)
    n_elements = max(len(l) for l in layouts)
    print(f"  Max elements per station  : {n_elements}")

    # --- COORDINATE_AXES: ENU→ITRF rotation per station ---
    # rows = [East, North, Up] unit vectors in ITRF  →  cell shape (3, 3)
    coordinate_axes = np.array(
        [enu_to_itrf_matrix(positions[i]) for i in range(n_stations)],
        dtype=np.float64)   # (nStat, 3, 3)

    # --- ELEMENT_OFFSET: element positions in ITRF frame ---
    # Stored as (3, nElem) per cell [everybeam convention]:
    #   cell[comp, elem] = ITRF component comp of element elem
    # Convert ENU layout → ITRF:  itrf = layout @ M   (nElem,3)@(3,3) → (nElem,3)
    # Then transpose to (3, nElem) for the cell.
    element_offset_3xN = np.zeros((n_stations, 3, n_elements), dtype=np.float64)
    element_flag_2xN   = np.ones ((n_stations, 2, n_elements), dtype=bool)

    for i, layout in enumerate(layouts):
        ne     = len(layout)
        M      = coordinate_axes[i]
        itrf   = layout @ M                   # (ne, 3) in ITRF
        if oskar_element_order:
            # (nElem, 3) per cell [legacy]: transpose to (3, nElem) only if needed
            element_offset_3xN[i, :,  :ne] = itrf.T        # stored (3, ne)
        else:
            element_offset_3xN[i, :,  :ne] = itrf.T        # (3, ne)
        element_flag_2xN[i, :, :ne]       = False          # valid elements unflagged

    if oskar_element_order:
        # Legacy (nElem, 3/2): transpose the per-cell matrix
        off_data  = element_offset_3xN.transpose(0, 2, 1)   # (nStat, nElem, 3)
        flag_data = element_flag_2xN.transpose(0, 2, 1)     # (nStat, nElem, 2)
        off_shape  = [n_elements, 3]
        flag_shape = [n_elements, 2]
        print(f"  ELEMENT_OFFSET per-cell shape : ({n_elements}, 3)  [OSKAR/legacy]")
    else:
        off_data  = element_offset_3xN                       # (nStat, 3, nElem)
        flag_data = element_flag_2xN                         # (nStat, 2, nElem)
        off_shape  = [3, n_elements]
        flag_shape = [2, n_elements]
        print(f"  ELEMENT_OFFSET per-cell shape : (3, {n_elements})  [everybeam default]")

    # --- Create the PHASED_ARRAY table using python-casacore ---
    # makearrcoldesc / maketabdesc / table() from casacore.tables
    # putcol convention: data[s] = cell for row s  →  matches OSKAR exactly
    col_defs = maketabdesc([
        makescacoldesc("ANTENNA_ID",      0,     comment="Antenna/station index"),
        makearrcoldesc("POSITION",        0.0,   ndim=1, shape=[3],
                       comment="Station ITRF XYZ (m)"),
        makearrcoldesc("COORDINATE_AXES", 0.0,   ndim=2, shape=[3, 3],
                       comment="ENU→ITRF rotation (rows: E, N, U)"),
        makearrcoldesc("ELEMENT_OFFSET",  0.0,   ndim=2, shape=off_shape,
                       comment="Element offsets in ITRF (m)"),
        makearrcoldesc("ELEMENT_FLAG",    False,  ndim=2, shape=flag_shape,
                       comment="Element flag per polarisation"),
    ])

    pa_tb = pt(pa_path, col_defs, nrow=n_stations)
    pa_tb.putcol("ANTENNA_ID",      np.arange(n_stations, dtype=np.int32))
    pa_tb.putcol("POSITION",        positions)         # (nStat, 3)
    pa_tb.putcol("COORDINATE_AXES", coordinate_axes)  # (nStat, 3, 3)
    pa_tb.putcol("ELEMENT_OFFSET",  off_data)          # (nStat, 3|nElem, nElem|3)
    pa_tb.putcol("ELEMENT_FLAG",    flag_data)         # (nStat, 2|nElem, nElem|2)
    pa_tb.close()

    # Register the subtable as a keyword on the main MS table
    ms_tb = pt(ms_path, readonly=False)
    ms_tb.putkeyword("PHASED_ARRAY", f"Table: {os.path.abspath(pa_path)}")
    ms_tb.close()

    print(f"  PHASED_ARRAY written to {pa_path}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Add a PHASED_ARRAY subtable to MS(s) from an OSKAR .tm dir.")
    p.add_argument("ms", nargs="+", help="One or more Measurement Set paths.")
    p.add_argument("--tm", required=True,
                   help="OSKAR telescope model directory.")
    p.add_argument("--overwrite", action="store_true",
                   help="Overwrite existing PHASED_ARRAY tables.")
    p.add_argument("--oskar-element-order", action="store_true",
                   help="Write ELEMENT_OFFSET/FLAG with per-cell shape (nElem, 3/2) "
                        "[legacy]. Default (everybeam) is (3/2, nElem).")
    args = p.parse_args()

    for ms_path in args.ms:
        print(f"\nProcessing {ms_path}")
        add_phased_array_table(ms_path, args.tm, overwrite=args.overwrite,
                               oskar_element_order=args.oskar_element_order)
    print("\nDone.")


if __name__ == "__main__":
    main()
