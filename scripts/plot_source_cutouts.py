#!/usr/bin/env python3
"""
plot_source_cutouts.py - grid of cutouts around sources, one column per image.

Rows    = sources, taken from a ds9 region file of circles (centre = source,
          diameter = cutout field of view, times --pad).
Columns = images, e.g. the app.restored images of the successive
          self-calibration steps of ddf-stream / pipeline_noWriteCol.

Colour scale (--scale):
  global (default)  one scale for every panel, as dsm.py does in ds9: sigma = std
                    of 10000 random pixels of the LAST image, linear from
                    -10 sigma to +40 sigma (--vmin-sigma / --vmax-sigma)
  row               one scale per row (source), from -nsigma x rms to the row peak
  panel             one scale per panel
Each panel shows its peak and its rms (estimated outside the circle), and the
restoring beam when the header has BMAJ/BMIN/BPA.

Facets (--tessel, default <stem of the last image>.tessel.reg when it exists): the
DDFacet tessel region file (ds9 lines + labelled points) is drawn on the cutouts
(instead of the circles; --circle to keep them) and on an overview PNG of the
whole image (<out>_tessel.png, image --overview, default the last one).

Examples
--------
# Default steps of the pipeline, looked up in the current directory:
plot_source_cutouts.py sources.reg -o cutouts.png

# Same, in another directory, with the MFS restored images:
plot_source_cutouts.py sources.reg --dir /path/to/run --suffix .app.restored_mfs.fits

# Explicit images with labels (label=path):
plot_source_cutouts.py sources.reg \
    "DI=image_dirin_SSD_m_c.app.restored.fits" \
    "DI-corr=image_dirin_SSD_m_c_di_m.app.restored.fits" \
    "phase1=image_phase1.app.restored.fits" \
    "ampphase1=image_ampphase1.app.restored.fits" -o cutouts.png

# Several pages (6 sources per PNG) and one PNG per source as well:
plot_source_cutouts.py sources.reg --per-page 6 --each -o cutouts.png

# Facets from a given tessel file, overview of the phase1 image:
plot_source_cutouts.py sources.reg --suffix .app.restored_mfs.fits \
    --tessel image_phase1.tessel.reg --overview 2 -o cutouts.png

Region file
-----------
ds9 format, circles only (other shapes are ignored), e.g.

    # Region file format: DS9 version 4.1
    fk5
    circle(10:20:30.5,-30:12:05.2,45") # text={S1}
    circle(155.31,-29.95,0.02)          # degrees when no unit (fk5)
    circle(183.2,-31.1,1.5')

Sexagesimal or decimal degrees; radius in " ' d (bare number = degrees in a
sky frame, pixels in 'image' frame). text={...} gives the row label.

Requires numpy, matplotlib, astropy (cmocean optional, for --cmap cmo.<name>).
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.visualization import ImageNormalize, AsinhStretch, LinearStretch

DEFAULT_STEPS = [
    ("DI", "image_dirin_SSD_m_c"),
    ("DI-corrected", "image_dirin_SSD_m_c_di_m"),
    ("phase1 (DDS0)", "image_phase1"),
    ("ampphase1 (Fast+Slow)", "image_ampphase1"),
]

SKY_FRAMES = ("fk5", "icrs", "fk4", "galactic", "ecliptic", "j2000")

# colormaps with a dark background: overlays in light colours
DARK_CMAPS = ("afmhot", "gist_heat", "hot", "inferno", "magma", "plasma", "viridis",
              "cividis", "cubehelix", "gnuplot", "gnuplot2", "cmo.thermal", "cmo.solar",
              "cmo.matter_r", "cmo.amp_r", "cmo.haline", "gray", "bone")


# --------------------------------------------------------------------------
# Region file
# --------------------------------------------------------------------------
def _parse_radius(tok, frame):
    tok = tok.strip()
    if tok.endswith('"'):
        return float(tok[:-1]) * u.arcsec
    if tok.endswith("'"):
        return float(tok[:-1]) * u.arcmin
    if tok.endswith("d"):
        return float(tok[:-1]) * u.deg
    if tok.endswith("r"):
        return (float(tok[:-1]) * u.rad).to(u.deg)
    if tok.endswith("p") or frame == "image":
        return float(tok.rstrip("p")) * u.pix
    return float(tok) * u.deg


def read_ds9_circles(path):
    """Return a list of dicts: name, coord (SkyCoord) or pix (x,y), radius (Quantity)."""
    frame = "fk5"
    out = []
    circ = re.compile(r"^\s*(?:(fk5|icrs|fk4|galactic|ecliptic|image|j2000)\s*;\s*)?"
                      r"circle\s*\(([^)]*)\)(.*)$", re.I)
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        low = s.lower()
        if low in SKY_FRAMES + ("image", "physical"):
            frame = low
            continue
        m = circ.match(s)
        if not m:
            continue
        fr = (m.group(1) or frame).lower()
        args = [a.strip() for a in m.group(2).split(",")]
        if len(args) < 3:
            continue
        name = None
        t = re.search(r"text\s*=\s*[{\"']([^}\"']*)[}\"']", m.group(3))
        if t:
            name = t.group(1).strip()
        rad = _parse_radius(args[2], fr)
        src = {"name": name or "S%i" % (len(out) + 1), "radius": rad}
        if fr == "image":
            src["pix"] = (float(args[0]) - 1, float(args[1]) - 1)  # ds9 is 1-based
        else:
            sexa = ":" in args[0] or "h" in args[0].lower()
            unit = (u.hourangle, u.deg) if sexa else (u.deg, u.deg)
            frm = "fk5" if fr == "j2000" else fr
            src["coord"] = SkyCoord(args[0], args[1], unit=unit, frame=frm)
        out.append(src)
    if not out:
        sys.exit("No circle found in %s" % path)
    return out


def read_tessel_reg(path):
    """DDFacet .tessel.reg (ds9, fk5 degrees): facet edges as line(ra0,dec0,ra1,dec1)
    and labelled points point(ra,dec) # text={...}. Returns (lines (N,4) deg, [(ra,dec,text)])"""
    lines, labels = [], []
    rl = re.compile(r"^\s*line\s*\(([^)]*)\)", re.I)
    rp = re.compile(r"^\s*point\s*\(([^)]*)\)(.*)$", re.I)
    for line in open(path):
        m = rl.match(line)
        if m:
            v = [float(x) for x in m.group(1).split(",")[:4]]
            if len(v) == 4:
                lines.append(v)
            continue
        m = rp.match(line)
        if m:
            v = [float(x) for x in m.group(1).split(",")[:2]]
            t = re.search(r"text\s*=\s*[{\"']([^}\"']*)[}\"']", m.group(2))
            labels.append((v[0], v[1], t.group(1) if t else ""))
    return np.array(lines, np.float64).reshape((-1, 4)), labels


def draw_tessel(ax, wcs, tessel, color, fontsize, shape=None):
    """Facet edges and labels of `tessel` on axis `ax` whose pixel frame is given by `wcs`"""
    lines, labels = tessel
    if lines.size:
        x0, y0 = wcs.wcs_world2pix(lines[:, 0], lines[:, 1], 0)
        x1, y1 = wcs.wcs_world2pix(lines[:, 2], lines[:, 3], 0)
        segs = np.stack([np.stack([x0, y0], -1), np.stack([x1, y1], -1)], 1)
        ok = np.all(np.isfinite(segs.reshape(len(segs), -1)), axis=1)
        from matplotlib.collections import LineCollection
        ax.add_collection(LineCollection(segs[ok], colors=color, linewidths=0.7,
                                         linestyles="--", alpha=0.85))
    for ra, dec, txt in labels:
        if not txt:
            continue
        x, y = wcs.wcs_world2pix([ra], [dec], 0)
        x, y = float(x[0]), float(y[0])
        if not (np.isfinite(x) and np.isfinite(y)):
            continue
        if shape is not None and not (0 <= x < shape[1] and 0 <= y < shape[0]):
            continue
        ax.plot([x], [y], "o", ms=2.5, mfc="none", mec=color, mew=0.7)
        ax.text(x, y, " " + txt, color=color, fontsize=fontsize, va="bottom", ha="left", clip_on=True)


# --------------------------------------------------------------------------
# Images
# --------------------------------------------------------------------------
class Image:
    """Lazily-read 2D plane of a (possibly 4D) FITS image."""

    def __init__(self, path, label):
        self.path, self.label = path, label
        self.hdul = fits.open(path, memmap=True)
        h = self.hdul[0].header
        self.header = h
        self.wcs = WCS(h).celestial
        d = self.hdul[0].data
        while d.ndim > 2:            # first frequency / Stokes plane
            d = d[0]
        self.data = d                # memmap view: only the cutouts are read
        self.pixscale = proj_plane_pixel_scales(self.wcs) * u.deg  # deg/pix
        self.bunit = h.get("BUNIT", "Jy/beam")
        self.beam = None
        if "BMAJ" in h and "BMIN" in h:
            self.beam = (h["BMAJ"] * u.deg, h["BMIN"] * u.deg, h.get("BPA", 0.0))

    def cutout(self, src, pad):
        if "pix" in src:
            pos = src["pix"]
            rad_pix = src["radius"].value
        else:
            pos = src["coord"]
            rad_pix = (src["radius"] / self.pixscale[1]).decompose().value
        size = int(np.ceil(2 * rad_pix * pad)) | 1
        try:
            c = Cutout2D(self.data, pos, (size, size), wcs=self.wcs,
                         mode="partial", fill_value=np.nan, copy=True)
        except Exception:            # source outside this image
            return None, rad_pix
        return c, rad_pix

    def binned(self, maxpix=2000):
        """Whole image binned by an integer factor to <= maxpix per side (max per block,
        so that sources stay visible); read in strips to bound memory"""
        ny, nx = self.data.shape
        f = max(1, int(np.ceil(max(nx, ny) / float(maxpix))))
        by, bx = ny // f, nx // f
        out = np.empty((by, bx), np.float32)
        for i in range(by):
            blk = np.asarray(self.data[i * f:(i + 1) * f, :bx * f], dtype=np.float32)
            blk = blk.reshape(f, bx, f)
            with np.errstate(all="ignore"):
                out[i] = np.nanmax(np.nanmax(blk, axis=2), axis=0)
        return out, f

    def dsm_sigma(self, n=10000, seed=0):
        """As dsm.py: std of n random pixels of the whole image (NaNs dropped)"""
        rng = np.random.default_rng(seed)
        flat = self.data.reshape(-1)
        ind = np.sort(rng.integers(0, flat.size, n))
        a = np.asarray(flat[ind], dtype=np.float64)
        a = a[np.isfinite(a)]
        return float(np.std(a)) if a.size else np.nan


def robust_stats(arr, centre, rad_pix):
    """Peak inside the circle, rms (MAD) outside it."""
    a = np.asarray(arr, dtype=np.float64)
    yy, xx = np.indices(a.shape)
    r = np.hypot(xx - centre[0], yy - centre[1])
    inside = (r <= rad_pix) & np.isfinite(a)
    outside = (r > rad_pix) & np.isfinite(a)
    peak = np.nanmax(a[inside]) if inside.any() else np.nan
    ref = a[outside] if outside.sum() > 50 else a[np.isfinite(a)]
    if ref.size == 0:
        return peak, np.nan
    rms = 1.4826 * np.median(np.abs(ref - np.median(ref)))
    return peak, rms


def fmt_flux(v, unit="Jy"):
    if not np.isfinite(v):
        return "n/a"
    a = abs(v)
    if a >= 1:
        return "%.2f %s" % (v, unit)
    if a >= 1e-3:
        return "%.2f m%s" % (v * 1e3, unit)
    return "%.1f µ%s" % (v * 1e6, unit)


def fmt_angle(q):
    q = q.to(u.arcsec)
    if q.value >= 120:
        return "%.1f′" % q.to(u.arcmin).value
    return "%.0f″" % q.value


def get_cmap(name):
    if name.startswith("cmo."):
        try:
            import cmocean  # noqa: F401  (registers the cmo.* colormaps)
        except ImportError:
            print("WARNING: cmocean not installed, using afmhot instead of %s" % name)
            return "afmhot"
    return name


# --------------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------------
def make_norm(vmin, vmax, args):
    stretch = LinearStretch() if args.stretch == "linear" else AsinhStretch(args.asinh)
    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)


def plot_page(sources, images, args, outname, global_lims=None):
    nr, nc = len(sources), len(images)
    pw = args.panel
    fig, axes = plt.subplots(nr, nc, figsize=(pw * nc + 0.9, pw * nr + 0.5),
                             squeeze=False,
                             gridspec_kw=dict(wspace=0.04, hspace=0.06))
    H = pw * nr + 0.5 + (0.35 if args.title else 0)
    fig.set_figheight(H)
    fig.subplots_adjust(left=0.07, right=0.90,
                        top=1 - (0.45 + (0.35 if args.title else 0)) / H, bottom=0.02)
    textbg = "white" if not args.dark else "black"

    mappable = None
    for i, src in enumerate(sources):
        cuts = []
        for im in images:
            c, rad_pix = im.cutout(src, args.pad)
            if c is None:
                cuts.append((None, None, np.nan, np.nan, rad_pix))
                continue
            ctr = c.position_cutout
            peak, rms = robust_stats(c.data, ctr, rad_pix)
            cuts.append((c, ctr, peak, rms, rad_pix))

        peaks = [x[2] for x in cuts if np.isfinite(x[2])]
        rmss = [x[3] for x in cuts if np.isfinite(x[3])]
        row_vmax = max(peaks) if peaks else 1.0
        row_rms = np.nanmin(rmss) if rmss else row_vmax / 100.0

        row_mappable = None
        for j, (im, (c, ctr, peak, rms, rad_pix)) in enumerate(zip(images, cuts)):
            ax = axes[i, j]
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_linewidth(0.5); sp.set_color("0.6")
            if c is None:
                ax.set_facecolor("0.92")
                ax.text(0.5, 0.5, "outside image", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8, color="0.4")
            else:
                if args.scale == "global":
                    vmin, vmax = global_lims
                elif args.scale == "panel" and np.isfinite(peak) and np.isfinite(rms):
                    vmin, vmax = -args.nsigma * rms, peak
                else:
                    vmin, vmax = -args.nsigma * row_rms, row_vmax
                if args.vmax_frac and args.scale != "global":
                    vmax *= args.vmax_frac
                row_mappable = ax.imshow(c.data, origin="lower", cmap=args.cmap,
                                         norm=make_norm(vmin, vmax, args),
                                         interpolation="nearest")
                mappable = row_mappable
                if args.show_circle:
                    ax.add_patch(Circle(ctr, rad_pix, fill=False, ls="--", lw=0.8,
                                        ec=args.overlay, alpha=0.8))
                if args.tessel_data is not None:
                    draw_tessel(ax, c.wcs, args.tessel_data, args.tessel_color,
                                args.fontsize - 2.5, shape=c.data.shape)
                if im.beam is not None:
                    bmaj = (im.beam[0] / im.pixscale[1]).decompose().value
                    bmin = (im.beam[1] / im.pixscale[1]).decompose().value
                    ny, nx = c.data.shape
                    off = max(bmaj, 0.06 * nx) * 0.8
                    ax.add_patch(Ellipse((off, off), bmaj, bmin,
                                         angle=90 + im.beam[2], fc="none",
                                         ec=args.overlay, lw=0.9))
                ax.text(0.03, 0.97, "peak %s\nrms %s" % (fmt_flux(peak),
                                                       fmt_flux(rms)),
                        transform=ax.transAxes, va="top", ha="left",
                        fontsize=args.fontsize - 1.5, color=args.overlay,
                        family="monospace",
                        bbox=dict(boxstyle="round,pad=0.2", lw=0, alpha=0.55, fc=textbg))
                ax.set_xlim(-0.5, c.data.shape[1] - 0.5)
                ax.set_ylim(-0.5, c.data.shape[0] - 0.5)
            if i == 0:
                ax.set_title(im.label, fontsize=args.fontsize)
            if j == 0:
                fov = 2 * src["radius"] * args.pad
                fovs = ("%.0f px" % fov.value) if fov.unit == u.pix else fmt_angle(fov)
                ax.set_ylabel("%s\n%s" % (src["name"], fovs), fontsize=args.fontsize)

        if row_mappable is not None and args.scale == "row":
            p0, p1 = axes[i, 0].get_position(), axes[i, -1].get_position()
            cax = fig.add_axes([p1.x1 + 0.008, p1.y0, 0.012, p1.height])
            cb = fig.colorbar(row_mappable, cax=cax)
            cb.ax.tick_params(labelsize=args.fontsize - 2.5, length=2)
            cb.set_label(images[0].bunit, fontsize=args.fontsize - 2.5)

    if mappable is not None and args.scale == "global":
        ptop, pbot = axes[0, -1].get_position(), axes[-1, -1].get_position()
        cax = fig.add_axes([ptop.x1 + 0.008, pbot.y0, 0.012, ptop.y1 - pbot.y0])
        cb = fig.colorbar(mappable, cax=cax)
        cb.ax.tick_params(labelsize=args.fontsize - 2, length=2)
        cb.set_label("%s  (%g to %g sigma, sigma = %s)" % (images[0].bunit, -args.vmin_sigma,
                                                           args.vmax_sigma, fmt_flux(args.sigma_used)),
                     fontsize=args.fontsize - 2)

    if args.title:
        fig.suptitle(args.title, fontsize=args.fontsize + 1,
                     y=1 - 0.12 / H)
    fig.savefig(outname, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", outname)


def plot_overview(im, args, outname, lims, sources):
    data, f = im.binned(args.overview_pix)
    w = im.wcs.deepcopy()
    # pixel frame of the binned image: block (i,j) covers original pixels [i*f, (i+1)*f)
    w.wcs.crpix = (w.wcs.crpix - 0.5) / f + 0.5
    if w.wcs.has_cd():
        w.wcs.cd = w.wcs.cd * f
    else:
        w.wcs.cdelt = w.wcs.cdelt * f
    if lims is None:
        sig = im.dsm_sigma()
        lims = (-args.vmin_sigma * sig, args.vmax_sigma * sig)
    fig, ax = plt.subplots(figsize=(args.overview_size, args.overview_size))
    mp = ax.imshow(data, origin="lower", cmap=args.cmap, interpolation="nearest",
                   norm=make_norm(lims[0], lims[1], args))
    if args.tessel_data is not None:
        draw_tessel(ax, w, args.tessel_data, args.tessel_color, args.fontsize - 1, shape=data.shape)
    for src in sources:   # where the cutouts are
        if "coord" in src:
            x, y = w.wcs_world2pix([src["coord"].ra.deg], [src["coord"].dec.deg], 0)
            ax.plot(x, y, "s", ms=7, mfc="none", mec=args.overlay, mew=1.0)
            ax.text(float(x[0]), float(y[0]), "  " + src["name"], color=args.overlay,
                    fontsize=args.fontsize, va="center")
    ax.set_xlim(-0.5, data.shape[1] - 0.5); ax.set_ylim(-0.5, data.shape[0] - 0.5)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("%s  (binned x%i, max per block)" % (im.label, f), fontsize=args.fontsize + 1)
    cb = fig.colorbar(mp, ax=ax, fraction=0.035, pad=0.01)
    cb.set_label(im.bunit, fontsize=args.fontsize)
    fig.savefig(outname, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", outname)


def build_images(args):
    images = []
    if args.images:
        for spec in args.images:
            if "=" in spec and not os.path.exists(spec):
                label, path = spec.split("=", 1)
            else:
                path = spec
                label = os.path.basename(spec)
                for sfx in (".app.restored_mfs.fits", ".app.restored.fits",
                            ".int.restored.fits", ".fits"):
                    if label.endswith(sfx):
                        label = label[:-len(sfx)]
                        break
            if not os.path.exists(path):
                print("WARNING: %s not found, skipped" % path)
                continue
            images.append(Image(path, label))
    else:
        for label, stem in DEFAULT_STEPS:
            path = os.path.join(args.dir, stem + args.suffix)
            if os.path.exists(path):
                images.append(Image(path, label))
            else:
                print("WARNING: %s not found, skipped" % path)
    if not images:
        sys.exit("No image to plot")
    return images


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("region", help="ds9 region file with circles")
    p.add_argument("images", nargs="*",
                   help="FITS images, optionally 'label=path'. Default: the "
                        "pipeline steps found in --dir")
    p.add_argument("--dir", default=".", help="directory of the default steps")
    p.add_argument("--suffix", default=".app.restored.fits",
                   help="suffix of the default steps (e.g. .app.restored_mfs.fits)")
    p.add_argument("-o", "--out", default="cutouts.png", help="output PNG")
    p.add_argument("--pad", type=float, default=1.5,
                   help="cutout side = pad x circle diameter (default 1.5)")
    p.add_argument("--scale", choices=("global", "row", "panel"), default="global",
                   help="global: one scale for all panels as dsm.py (default); "
                        "row: per source; panel: per panel")
    p.add_argument("--vmin-sigma", type=float, default=10.,
                   help="global scale: lower limit = -vmin_sigma x sigma (dsm.py: 10)")
    p.add_argument("--vmax-sigma", type=float, default=40.,
                   help="global scale: upper limit = vmax_sigma x sigma (dsm.py: 40)")
    p.add_argument("--sigma-image", type=int, default=-1,
                   help="global scale: index of the image sigma is measured on "
                        "(dsm.py: the last one, -1)")
    p.add_argument("--stretch", choices=("linear", "asinh"), default=None,
                   help="default: linear for --scale global (as dsm.py/ds9), asinh otherwise")
    p.add_argument("--nsigma", type=float, default=3.0,
                   help="row/panel scales: lower colour limit = -nsigma x rms (default 3)")
    p.add_argument("--vmax-frac", type=float, default=None,
                   help="row/panel scales: saturate at this fraction of the peak (e.g. 0.3)")
    p.add_argument("--asinh", type=float, default=0.05,
                   help="asinh stretch parameter (smaller = more contrast on faint)")
    p.add_argument("--cmap", default="afmhot",
                   help="colormap (default afmhot, close to ds9's 'bb' used by dsm.py; "
                        "e.g. inferno, magma, gist_heat, cmo.thermal with cmocean)")
    p.add_argument("--overlay", default=None,
                   help="colour of circles/text (default: auto for the colormap)")
    p.add_argument("--per-page", type=int, default=0,
                   help="sources per PNG (0 = all in one)")
    p.add_argument("--each", action="store_true",
                   help="also write one PNG per source (<out>_<name>.png)")
    p.add_argument("--panel", type=float, default=2.2, help="panel size in inches")
    p.add_argument("--fontsize", type=float, default=9)
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--title", default=None)
    p.add_argument("--tessel", default=None,
                   help="DDFacet .tessel.reg drawn on the cutouts and the overview "
                        "(default: <stem of the last image>.tessel.reg if it exists; 'none' to disable)")
    p.add_argument("--tessel-color", default="lime", help="colour of the facets (default lime)")
    p.add_argument("--circle", dest="show_circle", action="store_true",
                   help="also draw the region circles when facets are drawn")
    p.add_argument("--overview", type=int, default=-1,
                   help="index of the image shown whole with the facets in <out>_tessel.png "
                        "(default: the last one)")
    p.add_argument("--no-overview", action="store_true", help="do not write <out>_tessel.png")
    p.add_argument("--overview-pix", type=int, default=2000,
                   help="max size of the binned overview image in pixels (default 2000)")
    p.add_argument("--overview-size", type=float, default=12., help="overview size in inches")
    args = p.parse_args()
    args.cmap = get_cmap(args.cmap)
    if args.stretch is None:
        args.stretch = "linear" if args.scale == "global" else "asinh"
    args.dark = not args.cmap.endswith("_r") and (args.cmap in DARK_CMAPS or args.cmap == "afmhot")
    if args.overlay is None:
        args.overlay = "cyan" if args.dark else "tab:red"

    sources = read_ds9_circles(args.region)
    images = build_images(args)
    tessel = args.tessel
    if tessel is None:
        stem = images[-1].path
        for sfx in (".app.restored_mfs.fits", ".app.restored.fits", ".int.restored.fits",
                    ".app.residual.fits", ".fits"):
            if stem.endswith(sfx):
                stem = stem[:-len(sfx)]
                break
        tessel = stem + ".tessel.reg" if os.path.exists(stem + ".tessel.reg") else None
    elif tessel.lower() == "none":
        tessel = None
    args.tessel_data = None
    if tessel is not None:
        if not os.path.exists(tessel):
            sys.exit("tessel file %s not found" % tessel)
        args.tessel_data = read_tessel_reg(tessel)
        print("facets from %s: %i edges, %i labels" % (tessel, len(args.tessel_data[0]),
                                                      len(args.tessel_data[1])))
    else:
        args.show_circle = True
    keep = []
    for src in sources:
        if any(im.cutout(src, args.pad)[0] is not None for im in images):
            keep.append(src)
        else:
            print("WARNING: %s is outside every image, skipped" % src["name"])
    if not keep:
        sys.exit("No source inside the images")
    sources = keep
    print("%i source(s) x %i image(s): %s" % (len(sources), len(images),
                                           ", ".join(im.label for im in images)))

    global_lims = None
    args.sigma_used = np.nan
    if args.scale == "global":
        ref = images[args.sigma_image]
        sig = ref.dsm_sigma()
        args.sigma_used = sig
        global_lims = (-args.vmin_sigma * sig, args.vmax_sigma * sig)
        print("global scale (as dsm.py): sigma = %.3g on %s -> [%.3g, %.3g]"
              % (sig, ref.label, global_lims[0], global_lims[1]))

    base, ext = os.path.splitext(args.out)
    ext = ext or ".png"
    if not args.no_overview:
        plot_overview(images[args.overview], args, "%s_tessel%s" % (base, ext), global_lims, sources)
    n = args.per_page if args.per_page > 0 else len(sources)
    pages = [sources[k:k + n] for k in range(0, len(sources), n)]
    for k, page in enumerate(pages):
        name = args.out if len(pages) == 1 else "%s_p%i%s" % (base, k + 1, ext)
        plot_page(page, images, args, name, global_lims)
    if args.each:
        for src in sources:
            safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", src["name"])
            plot_page([src], images, args, "%s_%s%s" % (base, safe, ext), global_lims)


if __name__ == "__main__":
    main()
