#!/usr/bin/env python3
"""
Generate a 3D-printable STL heightmap from a PSRCHIVE archive.

The dynamic spectrum (frequency × time) is read WITHOUT dedispersion so the
dispersion sweep diagonal is preserved as a dramatic ridge across the model.

Each pixel of the dynamic spectrum becomes a point on a height-map surface.
The model has:
  - A solid flat base (base_mm thick)
  - The signal extruded upward: height ∝ intensity (clipped to ≥0)
  - Border walls so the base forms a complete enclosed solid

Output is a binary STL file ready to slice and print.

Usage:
    python ar_to_stl.py <archive> [options]

Examples:
    python ar_to_stl.py frb010724/FRB010724/ar_files/6.0
    python ar_to_stl.py burst.ar --width 120 --depth 80 --signal-height 20
    python ar_to_stl.py burst.ar -f 2 -t 4 --smooth 1.5
"""

import sys, os, argparse, struct
import numpy as np
from scipy.ndimage import gaussian_filter, uniform_filter

try:
    import psrchive
except ImportError:
    print("ERROR: psrchive Python bindings not found.")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Binary STL writer  (no external library needed)
# ---------------------------------------------------------------------------

def write_stl(filename, triangles, normals=None):
    """
    Write a binary STL file.

    Parameters
    ----------
    filename  : str
    triangles : (N, 3, 3) float32  – N triangles, each with 3 vertices (x,y,z)
    normals   : (N, 3) float32 or None  – per-triangle normals (computed if None)
    """
    n = len(triangles)
    if normals is None:
        v0 = triangles[:, 0, :]
        v1 = triangles[:, 1, :]
        v2 = triangles[:, 2, :]
        e1 = v1 - v0
        e2 = v2 - v0
        normals = np.cross(e1, e2)
        norms   = np.linalg.norm(normals, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        normals = normals / norms

    with open(filename, 'wb') as f:
        f.write(b'\x00' * 80)          # 80-byte header
        f.write(struct.pack('<I', n))   # triangle count
        for i in range(n):
            f.write(struct.pack('<3f', *normals[i]))
            for v in range(3):
                f.write(struct.pack('<3f', *triangles[i, v]))
            f.write(b'\x00\x00')       # attribute byte count
    print(f"  Written {n:,} triangles → {filename}")


# ---------------------------------------------------------------------------
# Heightmap → solid mesh
# ---------------------------------------------------------------------------

def heightmap_to_solid_stl(Z, x_mm, y_mm, base_mm=2.0):
    """
    Convert a 2-D height array Z[row, col] into a solid (watertight) STL mesh.

    Z       : (nrows, ncols) float  – height ABOVE base top surface in mm
    x_mm    : (ncols,) float        – x coordinates (time axis)
    y_mm    : (nrows,) float        – y coordinates (frequency axis)
    base_mm : float                 – thickness of flat base below Z=0
        The solid has:
            top surface   – the heightmap
      bottom face   – flat at z = -base_mm
      four side walls connecting top edge to bottom edge
    """
    nrows, ncols = Z.shape
    tris = []

    # ---- top surface (2 triangles per quad) ----
    for r in range(nrows - 1):
        for c in range(ncols - 1):
            x0, x1 = x_mm[c],   x_mm[c+1]
            y0, y1 = y_mm[r],   y_mm[r+1]
            z00 = Z[r,   c]
            z10 = Z[r+1, c]
            z01 = Z[r,   c+1]
            z11 = Z[r+1, c+1]
            # two triangles per grid cell (CCW from above = normal pointing up)
            tris.append([[x0,y0,z00],[x1,y0,z01],[x0,y1,z10]])
            tris.append([[x1,y0,z01],[x1,y1,z11],[x0,y1,z10]])

    # ---- bottom face (flat at z = -base_mm, normal pointing down) ----
    xmin, xmax = x_mm[0],  x_mm[-1]
    ymin, ymax = y_mm[0],  y_mm[-1]
    zb = -base_mm
    tris.append([[xmin,ymin,zb],[xmax,ymin,zb],[xmin,ymax,zb]])
    tris.append([[xmax,ymin,zb],[xmax,ymax,zb],[xmin,ymax,zb]])

    # ---- four side walls ----
    # Front wall  (y = ymin, c varies 0..ncols-1)
    for c in range(ncols - 1):
        x0, x1 = x_mm[c], x_mm[c+1]
        zt0, zt1 = Z[0, c], Z[0, c+1]
        tris.append([[x0,ymin,zb],  [x1,ymin,zb],  [x0,ymin,zt0]])
        tris.append([[x1,ymin,zb],  [x1,ymin,zt1], [x0,ymin,zt0]])

    # Back wall   (y = ymax)
    for c in range(ncols - 1):
        x0, x1 = x_mm[c], x_mm[c+1]
        zt0, zt1 = Z[-1, c], Z[-1, c+1]
        tris.append([[x1,ymax,zb],  [x0,ymax,zb],  [x0,ymax,zt0]])
        tris.append([[x1,ymax,zb],  [x1,ymax,zt1], [x0,ymax,zt0]])

    # Left wall   (x = xmin, r varies)
    for r in range(nrows - 1):
        y0, y1 = y_mm[r], y_mm[r+1]
        zt0, zt1 = Z[r, 0], Z[r+1, 0]
        tris.append([[xmin,y1,zb],  [xmin,y0,zb],  [xmin,y0,zt0]])
        tris.append([[xmin,y1,zb],  [xmin,y0,zt0], [xmin,y1,zt1]])

    # Right wall  (x = xmax)
    for r in range(nrows - 1):
        y0, y1 = y_mm[r], y_mm[r+1]
        zt0, zt1 = Z[r, -1], Z[r+1, -1]
        tris.append([[xmax,y0,zb],  [xmax,y1,zb],  [xmax,y0,zt0]])
        tris.append([[xmax,y0,zb],  [xmax,y0,zt0], [xmax,y1,zt1]])

    return np.array(tris, dtype=np.float32)


# ---------------------------------------------------------------------------
# Archive reader  (no dedispersion — keep the sweep)
# ---------------------------------------------------------------------------

def read_ar_nodedisp(filename):
    ar = psrchive.Archive_load(filename)
    ar.remove_baseline()
    ar.pscrunch()        # collapse to total intensity

    nsubint = ar.get_nsubint()
    nchan   = ar.get_nchan()
    nbin    = ar.get_nbin()

    all_data = []
    all_dur  = []
    for isub in range(nsubint):
        sub  = ar.get_Integration(isub)
        dur  = sub.get_duration()
        sd   = np.zeros((nchan, nbin), dtype=np.float32)
        for ic in range(nchan):
            prof = sub.get_Profile(0, ic)
            amps = np.array(prof.get_amps(), dtype=np.float32)
            if sub.get_weight(ic) == 0.0:
                amps[:] = 0.0
            sd[ic, :] = amps
        all_data.append(sd)
        all_dur.append(dur)

    data  = np.concatenate(all_data, axis=1) if nsubint > 1 else all_data[0]
    tsamp = all_dur[0] / nbin

    sub0  = ar.get_Integration(0)
    freqs = np.array([sub0.get_centre_frequency(ic) for ic in range(nchan)])
    times = np.arange(data.shape[1]) * tsamp * 1000  # ms

    weights = np.array([ar.get_Integration(0).get_weight(ic) for ic in range(nchan)])
    bad     = weights == 0.0

    print(f"  nsubint={nsubint}, nchan={nchan}, nbin={nbin}")
    print(f"  Freq range : {freqs.min():.1f} – {freqs.max():.1f} MHz")
    print(f"  Duration   : {times[-1]:.1f} ms")
    print(f"  Zero-weight channels: {bad.sum()}")

    return data, freqs, times, bad, tsamp


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='PSRCHIVE archive → 3D-printable STL heightmap (undedispersed)'
    )
    parser.add_argument('filename')
    parser.add_argument('-f', '--freq-downsample', type=int, default=1,
                        help='Frequency downsample factor (default: 1)')
    parser.add_argument('-t', '--time-downsample', type=int, default=1,
                        help='Time downsample factor (default: 1)')
    parser.add_argument('--flag-freq', type=str, default=None,
                        help='Flag bands e.g. "1540-1560,1200-1220"')
    parser.add_argument('--smooth', type=float, default=1.5,
                        help='Gaussian smoothing σ in pixels before extrusion (default: 1.5)')
    parser.add_argument('--noise-floor', type=float, default=1.0,
                        help='Clip normalised intensities below this σ to zero before '
                             'extrusion — flattens background noise (default: 1.0)')
    parser.add_argument('--width',  type=float, default=100.0,
                        help='Model width  (time axis) in mm (default: 100)')
    parser.add_argument('--depth',  type=float, default=80.0,
                        help='Model depth  (freq axis) in mm (default: 80)')
    parser.add_argument('--signal-height', type=float, default=15.0,
                        help='Max signal extrusion height in mm (default: 15)')
    parser.add_argument('--base',   type=float, default=2.0,
                        help='Base thickness in mm (default: 2)')
    parser.add_argument('--clip-low',  type=float, default=0.0,
                        help='Clip intensities below this percentile to 0 (default: 0)')
    parser.add_argument('--clip-high', type=float, default=99.5,
                        help='Clip intensities above this percentile (default: 99.5)')
    parser.add_argument('--time-start', type=float, default=None,
                        help='Start time in ms for pulse extraction (default: beginning)')
    parser.add_argument('--time-end', type=float, default=None,
                        help='End time in ms for pulse extraction (default: end)')
    parser.add_argument('-o', '--output', type=str, default=None)
    args = parser.parse_args()

    if not os.path.exists(args.filename):
        print(f"Error: '{args.filename}' not found"); sys.exit(1)

    FREQ_DS = args.freq_downsample
    TIME_DS = args.time_downsample

    # ------------------------------------------------------------------
    # 1. Load archive (NO dedispersion)
    # ------------------------------------------------------------------
    print(f"\nReading (no dedispersion): {args.filename}")
    data, freqs, times, bad_chan, tsamp = read_ar_nodedisp(args.filename)

    # ------------------------------------------------------------------
    # 1.5. Time window extraction
    # ------------------------------------------------------------------
    if args.time_start is not None or args.time_end is not None:
        t_start = args.time_start if args.time_start is not None else times[0]
        t_end = args.time_end if args.time_end is not None else times[-1]
        mask = (times >= t_start) & (times <= t_end)
        time_indices = np.where(mask)[0]
        if len(time_indices) > 0:
            t_idx_min, t_idx_max = time_indices[0], time_indices[-1] + 1
            data = data[:, t_idx_min:t_idx_max]
            times = times[t_idx_min:t_idx_max]
            print(f"  Extracted time window: {t_start:.3f} – {t_end:.3f} ms ({len(time_indices)} samples)")
        else:
            print(f"  Warning: no data in time range {t_start:.3f} – {t_end:.3f} ms")

    # ------------------------------------------------------------------
    # 2. Manual RFI flag
    # ------------------------------------------------------------------
    rfi_flag = bad_chan.copy()
    if args.flag_freq:
        for band in args.flag_freq.split(','):
            parts = band.strip().split('-')
            if len(parts) == 2:
                fl, fh = float(parts[0]), float(parts[1])
                if fl > fh: fl, fh = fh, fl
                m = (freqs >= fl) & (freqs <= fh)
                rfi_flag |= m
                print(f"  Flagged {fl}-{fh} MHz: {m.sum()} channels")
    data[rfi_flag, :] = 0.0

    # ------------------------------------------------------------------
    # 3. Downsample
    # ------------------------------------------------------------------
    def ds2d(d, ff, tf):
        nf, nt = d.shape
        nf2 = (nf // ff) * ff;  nt2 = (nt // tf) * tf
        return d[:nf2, :nt2].reshape(nf2//ff, ff, nt2//tf, tf).mean(axis=(1,3))
    def ds1d(a, f):
        n = len(a); n2 = (n//f)*f
        return a[:n2].reshape(n2//f, f).mean(axis=1)

    if FREQ_DS > 1 or TIME_DS > 1:
        print(f"\nDownsampling {data.shape} → ", end='')
        data  = ds2d(data, FREQ_DS, TIME_DS)
        freqs = ds1d(freqs, FREQ_DS)
        times = ds1d(times, TIME_DS)
        nf_trim = (len(rfi_flag) // FREQ_DS) * FREQ_DS
        rfi_flag = rfi_flag[:nf_trim].reshape(-1, FREQ_DS).any(axis=1)
        data[rfi_flag, :] = 0.0
        print(data.shape)

    # ------------------------------------------------------------------
    # 4. Per-channel normalisation
    #    Use the median of the lowest-intensity quartile of each channel
    #    as the baseline, so the dispersion sweep (which only occupies a
    #    small fraction of any given channel's time axis) does not inflate
    #    the baseline estimate and wash out the ridge.
    # ------------------------------------------------------------------
    nrows, ncols = data.shape   # (nchan, ntime)
    data_norm = np.zeros_like(data)
    for i in range(nrows):
        if rfi_flag[i]:
            continue
        ch   = data[i, :]
        # Estimate baseline from the quietest 50% of samples
        # (robust against the sweep brightening a few samples in this channel)
        sorted_ch = np.sort(ch)
        half      = max(2, len(sorted_ch) // 2)
        base      = np.median(sorted_ch[:half])
        off       = sorted_ch[:half] - base
        mad       = np.median(np.abs(off))
        sigma     = 1.4826 * mad if mad > 0 else (np.std(off) if np.std(off) > 0 else 1.0)
        data_norm[i, :] = (ch - base) / sigma

    # Zero out remaining bad channels
    data_norm[rfi_flag, :] = 0.0

    print(f"\nNormalised data: min={data_norm.min():.2f}  max={data_norm.max():.2f}")

    # ------------------------------------------------------------------
    # 5. Build height map Z[freq_idx, time_idx]
    #    - clip negatives to 0 (we only extrude positive signal)
    #    - percentile-scale to [0, signal_height_mm]
    # ------------------------------------------------------------------
    Z = data_norm.copy()

    # Apply noise floor BEFORE smoothing: zero out anything below threshold.
    # This eliminates the spiky 1-bit noise background so smoothing blends
    # genuine zeros rather than spreading spike energy into neighbours.
    if args.noise_floor > 0:
        Z[Z < args.noise_floor] = 0.0
        frac_zeroed = (Z == 0).mean() * 100
        print(f"  Noise floor {args.noise_floor}σ: zeroed {frac_zeroed:.1f}% of pixels")

    # Optional Gaussian smoothing (helps print quality)
    if args.smooth > 0:
        Z = gaussian_filter(Z, sigma=args.smooth)
        print(f"  Applied Gaussian smooth σ={args.smooth} px")

    # Clip to [clip_low_pct, clip_high_pct] of non-zero values
    nonzero = Z[Z != 0]
    if len(nonzero) > 0:
        lo = np.percentile(nonzero, args.clip_low)
        hi = np.percentile(nonzero, args.clip_high)
    else:
        lo, hi = 0.0, 1.0
    Z = np.clip(Z, lo, hi)

    # Shift so minimum is 0
    Z -= Z.min()

    # Scale to [0, signal_height_mm]
    zmax = Z.max()
    if zmax > 0:
        Z = Z / zmax * args.signal_height
    print(f"  Height map: 0 – {Z.max():.2f} mm  shape={Z.shape}")

    # ------------------------------------------------------------------
    # 6. Physical coordinate axes
    # freqs may be descending (foff < 0 for Parkes); normalise to ascending y
    freq_sorted = np.sort(freqs)   # ascending
    if freqs[0] > freqs[-1]:
        Z = Z[::-1, :]             # flip freq axis to match sorted order

    x_mm = np.linspace(0, args.width, ncols)   # time → x
    y_mm = np.linspace(0, args.depth, nrows)   # freq → y

    print(f"  Physical size: {args.width:.0f} × {args.depth:.0f} × "
          f"{args.base + args.signal_height:.0f} mm  (W × D × H)")

    # ------------------------------------------------------------------
    # 7. Build STL mesh
    # ------------------------------------------------------------------
    print(f"\nBuilding mesh (~{2*(nrows-1)*(ncols-1) + 2 + 4*(nrows+ncols-2):,} triangles) ...")
    tris = heightmap_to_solid_stl(Z, x_mm, y_mm, base_mm=args.base)

    # ------------------------------------------------------------------
    # 8. Write STL
    # ------------------------------------------------------------------
    base_out = os.path.splitext(os.path.basename(args.filename))[0]
    outname  = args.output or f"{base_out}_frb_heightmap.stl"
    print(f"\nWriting STL: {outname}")
    write_stl(outname, tris)

    print(f"\nDone.")
    print(f"  Model dimensions : {args.width:.0f} × {args.depth:.0f} × "
          f"{args.base + args.signal_height:.0f} mm")
    print(f"  Base thickness   : {args.base:.1f} mm")
    print(f"  Signal height    : {args.signal_height:.1f} mm")
    print(f"  Downsample       : freq×{FREQ_DS}  time×{TIME_DS}")
    print(f"  Grid resolution  : {nrows} × {ncols} points")
    print(f"\nPrint tips:")
    print(f"  - Layer height 0.15–0.2 mm works well for heightmaps")
    print(f"  - No supports needed (solid base)")
    print(f"  - Try 15–20% infill for the base; 100% for top 2 mm if you want crispness")
    print(f"  - If dispersion sweep detail is lost, increase --signal-height or reduce downsample")
    print(f"  Background noise level : --noise-floor (default 1.0σ, raise to 2.0 for flatter bg)")
    print(f"  Surface smoothness     : --smooth (default 1.5 px, raise to 2.0-3.0 for rounder ridges)")


if __name__ == '__main__':
    main()
