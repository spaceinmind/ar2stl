# ar2stl 🌌📡

**Convert a PSRCHIVE pulsar/FRB archive into a 3D-printable STL heightmap.**

The dynamic spectrum (frequency × time) is read **without dedispersion**, so the dispersion sweep diagonal is preserved as a dramatic physical ridge across the printed model. Each pixel of the dynamic spectrum becomes a point on the height-map surface, extruded upward in proportion to its intensity.


---

## Features

- Reads any PSRCHIVE-compatible archive (`.ar`)
- Preserves the dispersion sweep as a printable ridge — no dedispersion applied
- Per-channel robust normalisation (median/MAD, using the quietest 50 % of samples) keeps the sweep visible without baseline inflation
- Optional Gaussian smoothing for cleaner print surfaces
- Configurable noise floor to suppress background speckle
- Manual RFI flagging by frequency band
- Frequency and time downsampling
- Time-window extraction to focus on a specific pulse
- Produces a **watertight solid** STL (base + top surface + four walls) — no supports required when printing
- Pure Python + NumPy mesh writer; no external mesh library needed

---

## Requirements

See [`requirements.txt`](requirements.txt). The key non-standard dependency is the [PSRCHIVE](http://psrchive.sourceforge.net/) Python bindings, which must be built and installed separately (see below).

---

## Installation

### 1. Install PSRCHIVE with Python bindings

PSRCHIVE is not on PyPI. Follow the [official build instructions](http://psrchive.sourceforge.net/download.shtml), making sure to enable the Python bindings:

```bash
./configure --enable-shared PYTHON=$(which python3)
make
make install
```

### 2. Install Python dependencies

```bash
pip install -r requirements.txt
```

---

## Usage

```
python ar2stl.py <archive> [options]
```

### Quick examples

```bash
# Basic usage — all defaults
python ar2stl.py burst.ar

# Custom physical size and smoother surface
python ar2stl.py burst.ar --width 120 --depth 80 --signal-height 20 --smooth 2.0

# Extract a pulse window, downsample, flag RFI band
python ar2stl.py burst.ar \
    --time-start 42.0 --time-end 55.0 \
    -f 2 -t 4 \
    --flag-freq "1540-1560,1200-1220" \
    --noise-floor 2.0

# Specify output filename
python ar2stl.py frb010724/FRB010724/ar_files/6.0 -o frb010724_model.stl
```

---

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `-f`, `--freq-downsample` | `1` | Frequency downsample factor |
| `-t`, `--time-downsample` | `1` | Time downsample factor |
| `--flag-freq` | — | Flag RFI bands, e.g. `"1540-1560,1200-1220"` |
| `--smooth` | `1.5` | Gaussian smoothing σ in pixels before extrusion |
| `--noise-floor` | `1.0` | Zero out pixels below this σ level (flattens background noise) |
| `--width` | `100.0` | Model width along the time axis (mm) |
| `--depth` | `80.0` | Model depth along the frequency axis (mm) |
| `--signal-height` | `15.0` | Maximum signal extrusion height (mm) |
| `--base` | `2.0` | Base plate thickness (mm) |
| `--clip-low` | `0.0` | Clip intensities below this percentile to 0 |
| `--clip-high` | `99.5` | Clip intensities above this percentile |
| `--time-start` | — | Start of time window to extract (ms) |
| `--time-end` | — | End of time window to extract (ms) |
| `-o`, `--output` | `<archive>_frb_heightmap.stl` | Output STL filename |

---

## Output

A binary STL file containing a solid, watertight mesh:

| Component | Description |
|-----------|-------------|
| **Top surface** | Heightmap — intensity extruded upward |
| **Base plate** | Flat bottom at `z = -base_mm` |
| **Side walls** | Four walls connecting top surface edges to the base |

The mesh requires **no supports** when sliced with the base face down.

---

## 3D Printing Tips

- **Layer height**: 0.15–0.20 mm works well for heightmap detail
- **Infill**: 15–20 % for the base; 100 % for the top 2 mm if you want crisp ridge detail
- **No supports needed** — the solid base is self-supporting
- **Filament**: PLA or PETG both work; translucent filaments with a backlight look great for science outreach

### Tuning the model

| Problem | Solution |
|---------|----------|
| Dispersion ridge detail lost | Increase `--signal-height` or reduce `-f`/`-t` downsampling |
| Background too noisy / spiky | Raise `--noise-floor` (try `2.0`) |
| Surface too jagged | Raise `--smooth` (try `2.0`–`3.0`) |
| File too large to slice | Increase `-f` and/or `-t` |
| Ridge too flat | Lower `--noise-floor` or `--clip-high` |

---

## How it works

1. **Load** — Archive is loaded with `psrchive`, baseline-removed, and Stokes-I collapsed. No dedispersion is applied.
2. **RFI masking** — Zero-weight channels and any manually specified bands are blanked.
3. **Downsample** — Optional binning in frequency and time.
4. **Normalise** — Each channel is normalised using the median and MAD of its quietest 50 % of samples, making the per-channel baseline robust against the bright sweep.
5. **Noise floor** — Pixels below the threshold are zeroed before smoothing to prevent noise spikes from spreading.
6. **Smooth** — Gaussian filter sharpens the mesh aesthetically for printing.
7. **Scale** — Heights are clipped by percentile and linearly scaled to `[0, signal_height_mm]`.
8. **Mesh** — A watertight solid mesh is built (top surface + base + walls) and written as a binary STL.

---

## License

MIT

---

## Citation / Acknowledgements

If you use this tool in outreach or publications, a mention would be appreciated. Built on [PSRCHIVE](http://psrchive.sourceforge.net/) and NumPy.
