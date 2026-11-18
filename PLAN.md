# Plan

Status: **skeleton, for review.** One operation exists as a first
version (`astro:background`, below); the rest is design. The sources
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

1. Speed on 24 MP (bench as gegl-presence `tests/bench.sh`). Each cell is
   read whole and then subsampled; reading only the lattice rows would
   cut the reads by the step factor.
2. **RBF surface** (`model` = polynomial / thin plate spline, and a
   `smoothing` property): Siril's thin plate kernel `r^2 log r` with a
   diagonal smoothing term (`docs/research.md` §2.3); solve a system of
   the size of the used cells (at most 4096; dense Cholesky is fine up to
   about 1000, so a coarser grid or a smaller system for the RBF).
3. **User sample points:** a `points` string property (`"u,v;u,v"` in 0
   to 1 coordinates, saved in the XCF) that replaces the grid, and a
   small plug-in that fills it from GIMP's sample points (Pointer dialog)
   and opens the filter. Open question 4.
4. Leave out stars before the per-cell median with the shared star mask
   (op 2), for dense Milky Way fields where most cells hold stars.
5. An `amount` (0 to 100 %) for a partial removal.

### 2. `astro:star-mask`: Star Mask

Filters > Astronomy > Star Mask...

A mask of the stars, grey (white = star), to look at, to copy into a
layer mask, or to guide other GIMP tools. The same code (a shared
`operations/astro-common.h`) makes the mask inside ops 3 and 4.

| Property | Default | Meaning |
|---|---|---|
| `threshold` | 4.0 | detection level in noise sigmas above the local sky |
| `scale` | 1.5 px | size of the smallest stars (sigma of the detection filter) |
| `max-size` | 25 px | structures larger than this are not stars (nebula knots, galaxy cores) |
| `grow` | 1.0 px | dilation of the mask |
| `feather` | 1.5 px | Gaussian softening of the mask's edge |
| `show` | mask / overlay | the mask, or the image with the mask in red |

**Method:** luminance; local sky by a median (or the `astro:background`
surface) over a radius well above `max-size`; noise sigma by MAD of the
sky residual on a fixed lattice of the whole image; a difference of
Gaussians at `scale` against the residual; pixels above `threshold`
sigma; connected regions larger than `max-size` removed; grow, feather.
**Whole-image part:** only the noise statistic (whole bounding box for
the statistic, cached; the rest is an area filter with a margin of the
sky radius). Hennig has no star mask; the approach follows the published
PixInsight StarMask description (multiscale, threshold, growth,
smoothness) without its code.

### 3. `astro:star-reduce`: Reduce Stars

Filters > Astronomy > Reduce Stars...

| Property | Default | Meaning |
|---|---|---|
| `amount` | 70 % | blend towards the reduced stars |
| `radius` | 1.5 px | radius of the round structuring element |
| `method` | selection | erosion (min) / selection (rank at `selection`) / opening |
| `selection` | 0.25 | rank quantile: 0 = erosion, 0.5 = median, 1 = dilation (PixInsight's "Morphological Selection" as published) |
| `iterations` | 1 | 1 to 4 |
| `luminance-only` | true | shrink the brightness, keep each star's colour |
| star mask settings | as op 2 | `threshold`, `scale`, `max-size`, `grow`, `feather` |

**Method:** star mask m (op 2); rank filter over a disc of `radius`,
`iterations` times; `out = in + amount * m * (rank - in)`; in
luminance-only mode the channels are scaled by the new over the old
luminance (colour kept). Area filter; the mask statistic is
whole-image as in op 2. The directional variant (Hennig's Star Rounding)
is op 6.

### 4. `astro:star-color`: Star Colour

Filters > Astronomy > Star Colour...

| Property | Default | Meaning |
|---|---|---|
| `saturation` | 50 % | chroma boost inside the star mask |
| `background` | 1.0 | Siril's background factor: pixels darker than `(median + sigma) * background` stay as they are |
| `fill-cores` | true | clipped white cores take the colour of the star's halo |
| star mask settings | as op 2 | |

**Method:** chroma scaled (in the linear ratio space of gegl-presence, or
CIE LCh: open) by `1 + saturation * m`; `fill-cores`: for pixels at or
near clip in all channels inside a star, the chroma of the ring just
outside the clipped region, blended by the mask. The core fill is our
own idea (`docs/research.md` §2.4) and needs a look on real data.

### 5. `astro:stretch`: Astro Stretch

Filters > Astronomy > Stretch... (or Colors, open question 1)

| Property | Default | Meaning |
|---|---|---|
| `curve` | ghs | arcsinh / ghs |
| `black-point` | 0.0 | subtracted and rescaled first |
| `beta` | 10 | arcsinh strength (Lupton et al. 2004) |
| `stretch` (ln(D + 1)), `local-intensity` (b), `symmetry-point` (SP), `protect-shadows` (LP), `protect-highlights` (HP) | 0, 0, 0, 0, 1 | the five GHS parameters |
| `colour` | preserve | preserve (one factor from luminance for all channels) / per channel |
| `auto` | false | black point and symmetry point from the image's statistics |

**Method:** point operation; arcsinh `x' = x * asinh(beta Y) / (Y
asinh(beta))` on luminance Y; GHS from the published equations (all five
b cases, mirrored about SP, linear beyond LP and HP, normalised) or
Siril's `ght.c` (GPL-3.0-or-later): open question 6. **Whole-image part:**
only `auto` (median and MAD on a lattice, cached); without it a plain
point filter. Test: against an independent Python implementation of the
published equations.

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
  - star-mask: completeness and false detections against the known star
    list per signal to noise ratio; the share of each star's flux inside
    the mask.
  - star-reduce: FWHM before and after from second moments of isolated
    stars (target ratio per `radius`, `amount`); nebula flux outside the
    mask unchanged (< 0.1 %); no dark rings (the minimum around a star
    not below the local sky minus 3 noise sigma).
  - star-color: star chroma up by the set share; sky chroma changed by
    < 0.5 % of its noise.
  - stretch: closed-form values, monotonic, 0 and 1 fixed, GHS against an
    independent implementation of the published equations to 1e-6.
- **Every op:** GIMP 3.2.6 headless (`tests/gimp-check.sh`): NDE filter,
  every precision and grayscale, settings through XCF, the same pixels as
  plain GEGL and the gegl command line, a layer with offsets. ASan and
  UBSan (`tests/check.sh`). The dialog on Broadway (`tests/gui`).
  Everything through `tests/isolate.sh` with the snapshot of the user's
  folders.
- Timing on 24 MP (bench), with a target of under 1 s per preview for
  the point ops and under 3 s for the whole-image ones (open).

## Order (proposal)

1. `astro:background` next steps 1, 2, 4 (speed, RBF, star exclusion).
2. `astro:star-mask` with the shared code, then `astro:star-reduce`.
3. `astro:stretch` (small, point op).
4. `astro:star-color`.
5. Hennig's Star Rounding; the points plug-in; hot pixels if asked for.

## Open questions for you

1. **Menu:** everything under Filters > Astronomy (Hennig's path, one
   place for astro work), or the stretch under Colors next to Curves?
2. **Scope of the first round:** background plus star mask and star
   reduction, as proposed, or something else first?
3. **Star mask as its own op:** each op computes its own mask (GIMP makes
   aux inputs destructive), so op 2 is only for looking and for layer
   masks. Keep it as a visible op, or only inside ops 3 and 4?
4. **Sample points:** automatic grid only, or also user-placed points (a
   string property plus a small plug-in that reads GIMP's sample points)?
5. **RBF:** is a polynomial of degree up to 4 enough for your images, or
   is the thin plate spline (for local gradients, e.g. near the horizon)
   worth doing soon?
6. **GHS:** decided 2026-11-18: port Siril's `ght.c` (GPL-3.0-or-later,
   compatible), credited in the source and the README.
7. **Linear or as encoded:** the ops work in linear light. Most images
   that reach GIMP are already stretched (non-linear). Subtracting a
   gradient is physically right only on linear data; on a stretched
   image a fit in the encoded values may look better. Offer a
   `work-in` choice (linear / as encoded), or keep linear only?
8. **Real test images:** which ones may be used (CC0 or public domain,
   checked per image), e.g. a few NASA frames, or synthetic only?
9. **Hennig's Align and Merge Layers:** port as GIMP 3 plug-ins, or leave
   alignment and stacking to Siril?
10. **Names:** "Remove Sky Background" (now), or Siril's "Background
    Extraction", or "Gradient Removal"?
11. **Upstream and hosting:** keep as a separate GEGL module repository
    (as gegl-presence), and later offer ops to GEGL, or aim for GEGL from
    the start? (Nothing is pushed yet.)
