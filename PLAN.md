# Plan

Status: **in progress.** `astro:background` (first version),
`astro:stretch`, `astro:star-mask`, `astro:star-reduce` and
`astro:star-color` work; ops 6 and 7 and the plug-ins are design. Open questions 1 to 11 have
defaults chosen 2026-11-18 while David was away (end of this file). The sources
behind each method and licence note are in `docs/research.md`.

## Ground rules (from the research)

- **Single-input GEGL ops only.** GIMP 3.2.6 forces any filter with an
  `aux` pad to be destructive (`app/tools/gimpfiltertool.c` lines 431,
  593, 1583). So no op takes a star mask or a starless image as a second
  input: each op that needs a star mask computes it inside, with the
  same shared code and the same mask settings.
- **Menu:** all ops under `<Image>/Filters/Astronomy` (Georg Hennig's
  path for his GIMP 2 plug-ins). GIMP 3.2.6 creates the submenu for a GEGL
  op's `gimp:menu-path` (`app/widgets/gimpmenumodel.c` line 1241); **seen
  working** in GIMP 3.2.6 on Broadway (`docs/astronomy-menu.png`).
- **No neural nets** (StarNet's weights have no licence; GraXpert's AI
  mode is skipped). Classical methods only: morphology, robust
  statistics, surface fits. All older than 20 years in their basic form
  (`docs/research.md` §5).
- **Deterministic:** no random numbers in the ops; subsamples on fixed
  lattices; the same pixels give the same bits wherever the layer sits
  (checked for `astro:background`).
- **Linear light:** the ops work on `RGBA float` in the image's own RGB
  space (linear), as light adds up there. See open question 7.
- **Licence:** GPL-3.0-or-later, which takes both Hennig's code
  (GPL-2.0-or-later) and Siril's (GPL-3.0-or-later) with their copyright
  lines kept.
- **Property names:** never `input`, `output` or `aux` (GIMP's dialog
  leaves such a property out; found with the first op, now a check).

## The operations

### 1. `astro:background`: Remove Sky Background (first version done)

Filters > Astronomy > Remove Sky Background...

| Property | Type, range, default | Meaning |
|---|---|---|
| `correction` | subtract / divide, subtract | sky glow and light pollution add light; vignetting multiplies |
| `degree` | int 0 to 4, 2 | total degree of the polynomial |
| `samples` | int 4 to 64, 16 | cells along the longer side (square cells) |
| `tolerance` | 0.2 to 10, 2.0 | cells brighter than the median cell by more than this many robust sigmas are left out |
| `exclude-stars` | bool, false | leave the star mask's pixels (op 2, default settings) out of the cells |
| `keep-level` | bool, true | add the surface's mean back (or multiply by it) |
| `show` | corrected / model, corrected | the result or the fitted surface |

**Method (now):** per cell and channel a 3 sigma clipped median (3
rounds, sigma = 1.4826 MAD) of a lattice subsample of at most 64 x 64
pixels (Hennig's boxes); cells with under half opaque pixels are skipped;
bright cells dropped by MAD (as Siril); least squares on the cells in
coordinates scaled to -1 to 1 (Gaussian elimination with pivoting, in
double); up to 4 rounds dropping cells more than 3 robust sigmas off the
surface; the degree drops when too few cells are left. Then `in - B +
mean(B)` or `in / B * mean(B)` per channel; alpha kept.

**Why it is whole-image and how it stays fast:** the surface depends on
every cell, so `get_required_for_output`, `get_cached_region` and
`get_invalidated_by_change` return the whole input bounding box and
`threaded = FALSE` (the gegl-presence `dehaze.c` pattern): GEGL computes
the result once and caches it; a change anywhere recomputes everything.
The fit itself runs on a fixed-lattice subsample (the "downscaled copy":
at most 64 x 64 values per cell, so about 1 M values at the default 16
cells per side, however big the image) and on at most 64 x 64 cells; applying
it is one Horner evaluation per pixel and channel. No random numbers.

**Next for this op:**

1. **Done 2026-11-18:** speed on 24 MP (`tests/bench.sh`). Only the
   lattice rows of each cell are read, medians by quickselect (the same
   values as the sort, so every test number stayed the same), the
   correction applied in parallel (`gegl_parallel_distribute_area`).
   Defaults 0.467 s to 0.122 s; 32 cells, degree 4 1.436 s to 0.373 s.
2. **RBF surface: tried 2026-11-18, not shipped** (branch
   `spline-experiment`, with the probes in `tests/probes`). A smoothing
   thin plate spline (r^2 log r, terms 1, x, y, the smoothing on
   Siril's scale) through the cells, at most 32 per side, cells dropped
   when far from a plane through their 3 x 3 neighbours (a spline
   leaves no residuals to reject by), evaluated on an 8 px grid. The
   spline matched an independent numpy one (2.02 % at a corner of a
   noise-free quadratic, where cells start half a cell in: corners are
   extrapolated). Against the degree 4 polynomial on synthetic scenes
   with noise 0.003 (max and mean model error, share of the sky's
   amplitude):

   | Scene | Polynomial 4 | Spline | Both with tolerance 10 |
   |---|---|---|---|
   | horizon band (exp, 60 px) | 11.6 / 0.61 % | 41.3 / 2.77 % | 2.64 / 0.46 % vs 14.7 / 2.85 % |
   | nebula (sigma 55 px, as the gradient) | 6.46 / 2.08 % | 28.7 / 4.50 % | the same |

   The local rejection drops strongly curved cells (the band), and
   without it the spline follows nebulae. It needs better samples
   (user points, item 3, or samples placed away from bright areas, as
   Siril's automatic mode) before it can win. Left for later.

   **Found on the way, for the polynomial:** the bright cell rejection
   (`tolerance` 2) drops a strong light dome's own cells. A dome from
   beyond a corner (falling to 1/e in 0.3 widths): degree 4 within
   19.8 % (max) at tolerance 2, 7.7 % at 5, 4.0 % at 10; the horizon
   band 11.6 % at 2, 2.64 % at 5. The nebula (degree 2) is the same at
   tolerance 2 and 5 (2.61 %). So for strong light pollution a higher
   tolerance with degree 3 or 4 works better; the dialog's text could
   say so, or the default could move to about 5 (not changed: decide
   with real images).
3. **User sample points:** a `points` string property (`"u,v;u,v"` in 0
   to 1 coordinates, saved in the XCF) that replaces the grid, and a
   small plug-in that fills it from GIMP's sample points (Pointer dialog)
   and opens the filter. Open question 4.
4. **Done 2026-11-18:** `exclude-stars` leaves the pixels under the
   shared star mask (op 2, its default settings, no extra properties)
   out of each cell's sample. On a band of 9600 stars across the field
   (51 % of the pixels above the noise, like the Milky Way) the model is
   within 3.17 % of the gradient with the stars left out, 21.44 %
   without. Costs the mask: 1.07 s instead of 0.12 s at 24 MP. A cell
   counts when half of it is opaque and 16 values are left (stars no
   longer count against it).
5. An `amount` (0 to 100 %) for a partial removal.

### 2. `astro:star-mask`: Star Mask (done 2026-11-18)

Filters > Astronomy > Star Mask...

A mask of the stars, grey (white = star), to look at, to copy into a
layer mask, or to guide other GIMP tools. The same code (the shared
`operations/astro-common.h`) makes the mask inside ops 3 and 4, with the
same five settings (one macro, `ASTRO_STAR_MASK_PROPERTIES`).

| Property | Type, range, default | Meaning |
|---|---|---|
| `threshold` | 1 to 100, 4.0 | detection level in noise sigmas of the detection filter |
| `scale` | 0.5 to 10 px, 1.5 | size of the smallest stars (sigma of the detection filter) |
| `max-size` | 2 to 100 px, 15 | the largest star's FWHM; wider structures are not stars |
| `grow` | 0 to 20 px, 1.0 | dilation of the mask by a disc |
| `feather` | 0 to 20 px, 1.5 | Gaussian softening of the mask's edge |
| `show` | mask / overlay, mask | the mask, or the image with the mask in red |

**Method (as built; `astro-common.h` has the details):** luminance of
the image's RGB space; a difference of Gaussians `G(s) * L - G(4 s) * L`
(smooth sky, gradients and nebulae wider than a few s drop out, and
stars on a nebula stay separate); its noise by median and MAD on a
lattice of at most 512 x 512 points; pixels above `threshold` sigmas, in
8-connected regions; a region is a star when its FWHM (the area at or
above half its peak in the light above the local sky, as a disc) is at
most `max-size`; grow, feather. The local sky: L blurred at `max-size`,
L clipped at that plus 3 sigma, blurred again (so star light hardly
counts). Blurs are separable, threaded by GEGL, with mirrored edges.

**Changes from the first design, each found by a test:**

- **Whole-image**, not an area filter with a margin: a region's size is
  only known when all of it is seen. `get_required_for_output`,
  `get_cached_region`, `get_invalidated_by_change` return the whole
  layer, `threaded = FALSE` (the blurs thread themselves). Memory: three
  float planes of the layer at the peak (about 290 MB at 24 MP).
- **No separate median sky for detection:** the difference of Gaussians
  needs none.
- **Size as FWHM, not a second moment in the detection filter:** the
  filter limits what it sees to about 4 s, so every wide object looked
  star sized (a nebula knot of sigma 7 px passed). `max-size` is now the
  largest star's FWHM, default 15 px (was 25 px "structures"). A
  Gaussian knot of FWHM 18.8 px is left out; one of FWHM 1.1 x
  `max-size` reads about 10 % smaller than it is (the sky blur at
  `max-size` takes some of its light), so it is the borderline.
- **Edges:** with repeated edges a single noisy edge pixel dominated the
  blur, and false stars lined the image's border (most of the false
  detections in the first noise test). Now mirrored edges, and the
  detection filter divided by its own noise factor near the edges,
  which is exact for white noise and computed from the two kernels.

**Measured (tests/check.c and tests/check-stars.c, 600 x 400 fields with
a known star list):** at the default threshold 4, stars with a peak of
12, 20 and 50 noise sigmas and the bright ones: 343 of 343 found; at 8
sigmas 69 of 86, at 5 sigmas 29 of 86, at 3 sigmas 0 of 85. No false
regions among 600 stars; on pure noise 4.2 false regions per megapixel
(about 100 in a 24 MP frame; 171 per megapixel at threshold 3). The
share of a star's light under the mask, stars at 20 sigmas and up:
0.893. On a nebula (sigma 60 px, peak 0.15): mean mask away from stars
0.017; the knot (FWHM 18.8 px) not marked; 25 of 25 stars on the nebula
found. Transparent pixels never marked; bit identical for a moved
buffer. The blur matches a plain double sum with mirrored edges to
5.6e-7 (16 sizes and sigmas, kernels longer than the plane too).

### 3. `astro:star-reduce`: Reduce Stars (done 2026-11-18)

Filters > Astronomy > Reduce Stars...

| Property | Type, range, default | Meaning |
|---|---|---|
| `amount` | 0 to 100 %, 70 | blend towards the reduced stars |
| `radius` | 0.5 to 8 px, 1.5 | radius of the round structuring element |
| `method` | erosion / selection, selection | minimum, or the value at rank `selection` |
| `selection` | 0 to 1, 0.25 | rank quantile: 0 = erosion, 0.5 = median (PixInsight's "Morphological Selection" as published) |
| `iterations` | 1 to 4, 1 | |
| `luminance-only` | bool, true | "Keep star colour": shrink the brightness, keep each pixel's colour |
| star mask settings | as op 2 | `threshold`, `scale`, `max-size`, `grow`, `feather` |

**Method:** star mask m (op 2); per connected region of m > 0, the
region plus a margin of `radius` x `iterations` read from the input;
rank filter over a disc of `radius`, `iterations` times; `out = in +
amount * m * (rank - in)`, where `rank` is never taken below the local
sky of op 2, and only darker values are taken. In luminance-only mode
the rank filter runs on the luminance and the channels are scaled by the
new over the old luminance. Outside the mask nothing changes, bit for
bit.

**Changes from the first design, each found by a test:**

- **The opening is left out:** it flattens a star's peak, so the star
  got fainter but wider at half maximum (FWHM x 1.30, flux x 0.82 at
  amount 100).
- **The sky floor:** erosion of the sky's noise under the mask's edge
  pulled the ring around each star down (a dark ring): the minimum of
  the 9 noise values in the disc is on average 1.49 noise sigmas below
  their mean (computed). Now the result is never below the local sky.
- **FWHM from the half maximum, not second moments:** the second
  moments over a disc of 4 star sigmas weight the edge by r^2, so a sky
  offset of a few thousandths decides them (the darkened noise of an
  erosion made stars look 25 % smaller; a sky floor 0.0013 above the
  true sky made them look larger).

**Measured (bright stars and stars at 50 noise sigmas, median):**

| Settings | FWHM after / before | Flux after / before |
|---|---|---|
| default: selection 0.25, amount 70, radius 1.5 | 0.866 | 0.653 |
| erosion, amount 100, radius 1.5 | 0.661 | 0.328 |
| erosion, amount 100, radius 1.5, 2 iterations | 0.542 | 0.109 |

Outside the mask: no pixel changed (nebula flux outside the mask changed
by 0 %, limit 0.1 %); away from stars, under the mask's edges too, the
nebula's flux changed by -0.076 %. No dark rings: the median of the ring
2 to 6 star sigmas out after erosion at 100 % lies at most 0.30 noise
sigmas below the sky (noise 1e-4) and 0.24 (noise 0.003), limit 3.
Keep star colour: R / G of changed pixels kept to 1.2e-7. Amount 0
changes nothing; bit identical for a moved buffer.

### 4. `astro:star-color`: Star Colour (done 2026-11-18)

Filters > Astronomy > Star Colour...

| Property | Type, range, default | Meaning |
|---|---|---|
| `saturation` | 0 to 300 %, 50 | chroma boost inside the star mask |
| `background` | 0 to 10, 1.0 | Siril's background factor: pixels darker than `(median + sigma) * background` stay as they are (0: everything under the mask) |
| `fill-cores` | bool, true | clipped white cores take the colour of the star's halo |
| star mask settings | as op 2 | |

**Method:** each pixel's chroma in linear RGB, `c - Y` (Y its luminance
in the image's RGB space), times `1 + saturation * m * b`, so Y stays
exactly; the factor is limited so that no channel goes below 0. b: 0
below `(median + sigma) * background` (median and robust sigma of the
luminance on a lattice), rising to 1 over one more sigma (Siril's
saturation switches hard; `src/filters/saturation.c` read, written anew
here). The open question "linear ratio space or CIE LCh" is settled as
linear RGB about the luminance (the luminance is kept by construction).
`fill-cores` (our own idea, needs a look on real data): each region
under the mask with all channels at 0.95 or above takes the mean
chromaticity `c / Y` of the unclipped pixels within 2 px of it, scaled
so that its largest channel keeps its value, blended by the mask. Strips
of 128 rows in parallel; the clipped regions in a first pass;
whole-image (mask, statistics).

**Measured (tests/check-stars.c):** 171 star centres (20 sigmas and up)
get their chroma times `1 + 50 % x mask` to 2.1e-7; the luminance
changes by at most 4.4e-8 (without the core fill); the sky's chroma away
from stars changes by 0.209 % of the noise on average (limit 0.5 %;
under the feathered mask's edge noise above the background level is
boosted a little); outside the mask nothing changes; saturation 0
changes nothing; bit identical for a moved buffer. A star clipped to
white in its core (R / B 1.01) with a red halo (unclipped halo R / B
1.98): the core goes to R / B 1.81.

### 5. `astro:stretch`: Astro Stretch (done 2026-11-18)

Filters > Astronomy > Stretch... (title "Astro Stretch")

| Property | Type, range, default | Meaning |
|---|---|---|
| `curve` | ghs / arcsinh, ghs | the generalised hyperbolic stretch or the arcsinh stretch |
| `black-point` | 0 to 0.99, 0 | taken off first: `x' = max (0, (x - bp) / (1 - bp))`, clipped to 0 to 1 |
| `stretch` | 0 to 10, 0 | GHS ln(D + 1) (Siril's range); 0 leaves the image as it is |
| `local-intensity` | -5 to 15, 0 | GHS b (Siril's range) |
| `symmetry-point` | 0 to 1, 0 | GHS SP |
| `protect-shadows` | 0 to 1, 0 | GHS LP, kept at most SP |
| `protect-highlights` | 0 to 1, 1 | GHS HP, kept at least SP |
| `inverse` | bool, false | GHS inverse transform (undoes the same stretch) |
| `beta` | 0 to 1000, 10 | arcsinh strength (Lupton et al. 2004) |
| `colour` | preserve / per-channel, preserve | one factor from the luminance of the image's RGB space for all channels, or each channel alone |
| `auto-points` | bool, false | black point and symmetry point from the image (`auto` is a C keyword) |
| `encoded` | bool, false | "Work on encoded values (as Siril)": the curve on the image's encoded values, not on linear light (open question 7) |

The GHS sliders show only with GHS chosen, `beta` only with arcsinh
(`ui_meta "visible"`, as gegl-depth-blur).

**Method:** a point operation. GHS: Siril's `ght.c` (`GHTsetup`, `GHT`;
master `499fd209`, 2026-11-09) ported to double: the forward and the
inverse transform, all five b cases (b = -1, b < 0, b = 0, b > 0). The
stretch property is Siril's command value, `D = expm1 (stretch)`
(`command.c` line 3534). Arcsinh: `asinh (beta x) / asinh (beta)`.
Preserve colour: curve (Y) / Y times each channel (Lupton); a pixel
whose largest channel would pass 1 is scaled down whole (Siril's
"rescale" clip mode), so its hue stays. Alpha kept.

**A bug found in Siril while porting:** in the inverse transform for
b = -1, Siril sets the transformed symmetry point `SPT = q0 * q`; the
forward curve at SP gives `-q0 * q` (q0 is negative there), so Siril's
inverse never takes its lower branch and is wrong below SP (0.066 off at
0 for b = -1, ln(D+1) = 3, SP = 0.1). Fixed here (`-q0 * q`) and noted in
the source; found by the check against the published equations. Worth
reporting to Siril (not done: David's call).

**Whole-image part:** only `auto-points`: the median and robust sigma
(1.4826 MAD) of the luminance (preserve) or of each channel
(per-channel), on a fixed lattice of at most 256 x 256 points of the
whole layer, give `bp = max (0, median - 3 sigma)` and SP where the
median lands. Then `get_required_for_output` and
`get_invalidated_by_change` return the whole input; without it the op
asks only for the region it draws. The factor 3 is our own choice, not
taken from Siril (Siril's AutoGHS sets SP = median + k sigma with its
own k and leaves the black point).

**Measured (tests/check.c, tests/stretch-reference.py):** the ramp
outputs of 20 GHS cases (b from -5 to 15, ln(D+1) up to 10, LP = SP =
HP, SP = 1, LP and HP out of order, a black point, 5 inverse cases) and
4 arcsinh cases against an independent Python implementation of the
published equations (the inverse by bisection of the forward
equations): largest difference 9.4e-7 (b = 15, ln(D+1) = 10, at SP),
limit 1e-6; most cases 3e-8. 0 and 1 kept to 3e-17, all 24 curves
rising; forward then inverse back to 1.8e-6 where the curve is not
flattened below float resolution; preserve colour: channel ratios kept
to 1.6e-7, the luminance follows the curve to 1.5e-7, nothing above 1;
the automatic points agree with the same points set by hand (from the
whole image's median and MAD) to 0.9 % in the stretched sky median, and
0.08 % of the pixels go below the black point (Gaussian 3 sigma:
0.135 %); bit identical for a moved buffer. In GIMP (headless): every
precision and grayscale, settings through XCF, a layer with offsets, the
same pixels as the gegl command line.

**Encoded values (added 2026-11-18, PLAN question 7):** with `encoded`
the operation takes and gives `R'G'B'A float` in the image's space, so
the curve, the luminance weights and the automatic points all work on
the encoded numbers, which is what Siril does with the numbers of a
file. Checked (tests/check.c, `stretch_encoded_matches_numbers`): the
result read encoded equals the linear mode run on the encoded numbers
handed over as if they were linear, to 5.4e-7, for GHS preserve and per
channel, with a black point, arcsinh and the automatic points; off by
default (`stretch_encoded_off_by_default`); in GIMP headless (a
`stretch_encoded` case) the same pixels as plain GEGL and the gegl
command line on every precision.

What it matches: Siril shows a file's numbers as they are. An image
GIMP holds with the perceptual (sRGB) curve (Image > Precision:
"perceptual gamma", as an 8 bit JPEG or PNG loads) also shows its numbers as they are, so there the same settings give
Siril's picture. An image in "linear light" precision is already shown
brighter than Siril shows the same numbers (GIMP puts them through the
sRGB curve for the screen); for it, the Siril look would need "linear
in, encoded out", which is not built (decide when someone needs it).

### 6. `astro:star-round`: Round Stars (Hennig port, later)

Hennig's Star Rounding: radial (`a3 r^3 + a2 r^2 + a1 r`) and linear
(length, angle) elongation corrected by a directional min (or max) with
bilinear weights, above a threshold, blended by a percentage. Numbers
only, so a clean single-input area op; the port keeps © Georg Hennig.

### 7. `astro:hot-pixels`: Remove Hot Pixels (later, maybe)

Siril's cosmetic correction: deviation from the 5 x 5 median against k
image sigmas, replaced by the 3 x 3 mean (hot) or median (cold). Area op
with radius 2 plus a whole-image sigma. Mostly a raw frame task for Siril.

### Plug-ins (later, only if wanted)

- **Background points:** writes `astro:background`'s `points` from GIMP's
  sample points (op 1, next 3).
- **Stack layers:** Hennig's Merge Layers (median, sigma-clipped mean) for
  a few already aligned layers. Siril does stacking far better; open
  question 9.

## Test strategy

Synthetic data only in the repository (generated, fixed seeds). Real
images only CC0 or public domain with the licence checked on the source
page per image (NASA's are usually public domain, but each is checked;
never the user's photos).

- **Generator** (`tests/check.c`): Gaussian stars (sigma 0.8 to 2.5 px,
  peaks with many faint and few bright, colours red to blue), later
  Moffat PSFs with a known FWHM and beta; gradients: quadratic per
  channel, vignetting, later an exponential light dome and a
  non-polynomial one; a nebula (wide Gaussian); Gaussian noise.
- **Recovery numbers**, per op:
  - background (now): model within 0.54 % of the gradient's amplitude
    (limit 2 %), sky left after correction spans 0.29 % (limit 4 %), star
    flux kept to 0.05 % (limit 1 %), degree 2 to 4 within 0.5 % (limit
    3 %), a nebula of the gradient's size leaks 2.2 % (limit 6 %),
    vignetting divided out to 0.57 % (limit 3 %); flat sky untouched;
    transparent pixels ignored; bit identical for a moved buffer. Next:
    the share of a nebula's flux removed, against the naive fit.
  - star-mask (done, op 2): completeness and false detections against
    the known star list per signal to noise ratio; the share of each
    star's flux inside the mask.
  - star-reduce (done, op 3): FWHM before and after of isolated stars
    (from the half maximum; second moments proved unusable, op 3);
    nebula flux outside the mask unchanged (< 0.1 %); no dark rings (the
    median of the ring around a star not below the local sky minus 3
    noise sigma; the single darkest pixel, as first planned, was below
    that for 101 of 400 stars after an erosion of the noise, before the
    sky floor).
  - star-color (done, op 4): star chroma up by the set share; sky chroma
    changed by < 0.5 % of its noise.
  - stretch (done): closed-form values, monotonic, 0 and 1 fixed, GHS
    against an independent implementation of the published equations to
    1e-6 (measured: at most 9.4e-7; op 5 above).
- **Every op:** GIMP 3.2.6 headless (`tests/gimp-check.sh`): NDE filter,
  every precision and grayscale, settings through XCF, the same pixels as
  plain GEGL and the gegl command line, a layer with offsets. ASan and
  UBSan (`tests/check.sh`). The dialog on Broadway (`tests/gui`).
  Everything through `tests/isolate.sh` with the snapshot of the user's
  folders.
- Timing on 24 MP (`tests/bench.sh`, a synthetic star field of 6000 x
  4000 with 30000 stars, float, 24 GEGL threads on this machine), with a
  target of under 1 s per preview for the point ops and under 3 s for
  the whole-image ones. Measured 2026-11-18:

  | Operation | Time |
  |---|---|
  | gegl:gaussian-blur 40 (for scale) | 0.400 s |
  | Remove Sky Background, defaults | 0.122 s |
  | Remove Sky Background, 32 cells, degree 4 | 0.373 s |
  | Remove Sky Background, stars left out | 1.065 s |
  | Astro Stretch, GHS / arcsinh / automatic points | 0.073 / 0.079 / 0.099 s |
  | Star Mask, defaults | 1.023 s |
  | Reduce Stars, defaults | 1.192 s (3.315 s before the strips) |
  | Reduce Stars, erosion, 2 iterations | 1.314 s |
  | Star Colour, defaults | 1.089 s |

  Reduce Stars runs in strips of 128 rows in parallel, each reading its
  rows plus the disc's reach from the input and filtering only where the
  mask needs it; the result is bit identical to the first version (one
  region at a time). The star mask's time is mostly its blurs (two at
  sigma `max-size`) and the single-threaded region pass.

## Order (proposal)

1. `astro:background` next steps 1, 2, 4 (speed, RBF, star exclusion):
   1 and 4 done 2026-11-18; the RBF tried and not shipped (it lost to
   the polynomial on every synthetic scene; item 2 above).
2. `astro:star-mask` with the shared code, then `astro:star-reduce`:
   done 2026-11-18.
3. `astro:stretch` (small, point op): done 2026-11-18 (built first).
4. `astro:star-color`: done 2026-11-18.
5. Hennig's Star Rounding; the points plug-in; hot pixels if asked for.

## Open questions for you (defaults chosen while David was away)

1. **Menu:** everything under Filters > Astronomy (Hennig's path, one
   place for astro work), or the stretch under Colors next to Curves?
   **Default chosen 2026-11-18 while David was away; easy to change:** everything under Filters > Astronomy, the stretch too.
2. **Scope of the first round:** background plus star mask and star
   reduction, as proposed, or something else first?
   **Default chosen 2026-11-18 while David was away; easy to change:** the proposed first round (background, star mask, star
   reduction) plus the stretch.
3. **Star mask as its own op:** each op computes its own mask (GIMP makes
   aux inputs destructive), so op 2 is only for looking and for layer
   masks. Keep it as a visible op, or only inside ops 3 and 4?
   **Default chosen 2026-11-18 while David was away; easy to change:** keep `astro:star-mask` as a visible op.
4. **Sample points:** automatic grid only, or also user-placed points (a
   string property plus a small plug-in that reads GIMP's sample points)?
   **Default chosen 2026-11-18 while David was away; easy to change:** the automatic grid only for now.
5. **RBF:** is a polynomial of degree up to 4 enough for your images, or
   is the thin plate spline (for local gradients, e.g. near the horizon)
   worth doing soon?
   **Default chosen 2026-11-18 while David was away; easy to change:** the polynomial now, the RBF later.
   (Tried anyway: it lost to the polynomial on synthetic scenes, see op 1
   item 2; it needs better sampling first.)
6. **GHS:** decided 2026-11-18: port Siril's `ght.c` (GPL-3.0-or-later,
   compatible), credited in the source and the README.
7. **Linear or as encoded:** the ops work in linear light. Most images
   that reach GIMP are already stretched (non-linear). Subtracting a
   gradient is physically right only on linear data; on a stretched
   image a fit in the encoded values may look better. Offer a
   `work-in` choice (linear / as encoded), or keep linear only?
   **Default chosen 2026-11-18 while David was away; easy to change:** linear only for now. The trade-off: the background fit and
   the star ops are right on linear data, and wrong in a subtle way on
   an already stretched image (a gradient there is no longer additive).
   For `astro:stretch` it also means: the curve maps linear values to
   linear values, and GIMP then shows them through the image's TRC
   (sRGB), so the same Siril numbers give a brighter picture here than
   in Siril, which writes the stretched numbers out as the encoded
   values. A Siril-identical look needs either an "output as encoded"
   option on the stretch (input read linear, output written as R'G'B')
   or lower settings. Not built; decide with real data.
   **Built 2026-11-18 (evening), off by default:** the stretch's
   `encoded` option works on the encoded values (op 5). The background
   fit and the star ops stay linear only.
8. **Real test images:** which ones may be used (CC0 or public domain,
   checked per image), e.g. a few NASA frames, or synthetic only?
   **Default chosen 2026-11-18 while David was away; easy to change:** synthetic only in the repository.
9. **Hennig's Align and Merge Layers:** port as GIMP 3 plug-ins, or leave
   alignment and stacking to Siril?
   **Default chosen 2026-11-18 while David was away; easy to change:** leave alignment and stacking to Siril.
10. **Names:** "Remove Sky Background" (now), or Siril's "Background
    Extraction", or "Gradient Removal"?
    **Default chosen 2026-11-18 while David was away; easy to change:** keep "Remove Sky Background".
11. **Upstream and hosting:** keep as a separate GEGL module repository
    (as gegl-presence), and later offer ops to GEGL, or aim for GEGL from
    the start?
    **Default chosen 2026-11-18 while David was away; easy to change:** a separate repository (github.com/sandbranch/gegl-astro).
