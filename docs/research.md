# Research: an astrophotography toolkit for GIMP 3

Checked 2026-09-28. **Read** means the source was read (the commit or
tarball is named); **page** means a web page says so and the code was not
read; **not checked** says itself. Nothing here is legal advice; the
licence notes are what the files say.

Starting points: gimp-devtools `docs/photoshop-gaps.md` item 11 and
`docs/prior-art-candidates.md` ("Astro toolkit"): nothing for GIMP 3, no
G'MIC filter and no LinuxBeaver op for star reduction or gradient removal;
Siril and GraXpert do it outside GIMP.

## 1. Georg Hennig's gimp-plugin-astronomy

### Where the code is

- Author's page: <http://www.hennigbuam.de/georg/gimp.html> (HTTP only;
  HTTPS fails the TLS handshake). Latest release 0.10, source
  `downloads/gimp-plugin-astronomy-0.10.tar.bz2`; downloaded, MD5
  `5a70dcd8699790497ff45451c9c4071f`, the same as the page states. Built
  for GIMP 2.10.
- [gimp-plugins-justice/gimp-plugin-astronomy](https://github.com/gimp-plugins-justice/gimp-plugin-astronomy)
  (git, last push 2023-10-21, head `1fd3af2`): 0.10 plus build fixes by
  Alexandre Prokoudine and Joe Da Silva (autotools, gettext, a Mageia gcc 7
  patch). **Read.** `diff` of `src/*.c` against the tarball: only
  `%` escaped as `(%)` in a few help strings and one `static`; the
  plug-ins, parameters and algorithms are the same.
- Forks: JoesCat (2024-02-19, still `gimp-2.0` in `configure.ac`),
  funxiun (2019), OpenMandrivaAssociation (packaging, 2026-09-16). All
  GIMP 2: they use `gimp_pixel_rgn_*`, `gimp_install_procedure` and
  `gimp_plugin_menu_register`, which GIMP 3 removed.

### Licence

Every `src/*.c` and `scripts/*.scm` file: "Copyright (C) 2006-2018 by
Georg Hennig" (years vary per file), "either version 2 of the License, or
(at your option) any later version". `COPYING` is the GPL version 2 text.
So GPL-2.0-or-later, which may be taken under GPL-3.0-or-later: ported
code can go into this repository with its copyright line kept. The
plug-ins link GSL (GPL-3.0-or-later since GSL 1.12, not checked in this
session) and FFTW3 (GPL-2.0-or-later); a port would not need either.

### Every tool (0.10)

C plug-ins, all under `<Image>/Filters/Astronomy` except one:

| Tool | File | What it does |
|---|---|---|
| Align Layers | `alignment.c` | Moves and rotates layers onto the first one. Methods: Gauss fit, centre of brightness (search radius), cross correlation, cross correlation by FFT (search area in % of the selection, subpixel option). Two-star alignment for de-rotation; trim to the overlap; scale before processing; visible layers only; remove invisible layers. The user selects the star or object. |
| Align images with autopano | `alignment_autopano.c` | Under `<Toolbox>/Xtns`: opens images as layers and aligns them with the external autopano program. |
| Merge Layers | `merge.c` | Stacks the layers into a new image: arithmetic mean, geometric mean (both alpha weighted), median, mean of the values within median ± n sigma, sigma median in one or two passes. Multiplier and additive constant per pixel first. A preview. |
| Background Gradient | `background_gradient.c` | An artificial flat field. Splits the layer into boxes (width, height); per box and per channel sorts the values (alpha 0 skipped), optionally keeps only the darkest or brightest z %, then c iterations of clipping at median ± n sigma, and the mean of what is left if at least m values remain. The box means are the points of a least squares fit (GSL `gsl_multifit_linear`) of a **tensor product polynomial of degree 4 in x and 4 in y** (25 terms `x^i y^j`, i, j ≤ 4) in raw pixel coordinates. The polynomial is divided by its maximum (found along rows with a cubic solve) and written to a new 8 bit layer in Divide mode. The documentation: cut nebulae and galaxies out of a copy first; "images of the milky way obviously give bad results". |
| Star Rounding | `star_rounding.c` | Rounds stars stretched by guiding errors or field curvature: a radial pass, length `P(r) = a3 r^3 + a2 r^2 + a1 r` from the centre, and a linear pass (length L, angle phi). Per pixel above (or below) a brightness threshold, the minimum (or maximum) of the pixel and the pixel(s) at the end of the vector (L, phi), bilinear weights among the four neighbours, blended by a percentage. Brightness 0.30 R + 0.59 G + 0.11 B. A directional erosion. |
| Artificial Stars | `artificial_stars.c` | Renders a synthetic star distribution (distribution options on one side, drawing options on the other, sample presets). |
| Artificial Galaxy | `artificial_galaxy.c` | Renders a synthetic galaxy; only elliptical galaxies ("quite difficult to implement"). |

Script-Fu (`scripts/`): `background-gradient-batch` (the gradient plug-in
over a file glob), `border-information` (a frame with up to two colours and
text lines), `brightness-contrast-batch`, `dark-subtraction` (subtract a
dark frame from every layer), `flat-division` (divide every layer by a
flat), `mode-batch` (one layer mode for all layers), `normalize-batch`.

### What is worth porting

- **Background Gradient:** the per box statistics (sort, keep darkest z %,
  sigma clipping, minimum count) are a good sampling front end and small.
  The fit is not worth copying as is: raw pixel coordinates to the 8th
  total power are badly conditioned (x^4 y^4 at 4000 px is 2.6e28), and
  degree 4 per axis oscillates at the edges. Normalised coordinates and a
  total degree polynomial (as Siril) or a thin plate spline are better.
  The 8 bit Divide layer output is replaced by a float op.
- **Star Rounding:** fits as a single-input op (`astro:star-round`, later):
  properties are numbers only.
- **Merge Layers** and **Align Layers:** multi-layer tools, so plug-ins,
  and Siril already stacks and registers far better. Low priority; a
  sigma-clipped layer merge may still be a small, useful plug-in.
- **Artificial stars/galaxy:** useful as **test generators** here (known
  PSFs), not as user tools; our generator is written fresh in the tests.
- The Script-Fu batch scripts: GIMP's batch mode and BIMP cover them.

## 2. Methods

### 2.1 Star mask

The common base of star reduction and star colour. Classical, no neural
net:

- A luminance, a large-scale background estimate (median or blur over a
  radius well above the star size) subtracted, the noise sigma estimated
  robustly (MAD of the residual), pixels above k sigma are star
  candidates; optional a small-scale filter first (difference of Gaussians
  or the first à trous wavelet layers, as Siril's star finder and
  PixInsight StarMask do; PixInsight StarMask is described on
  [deepskycolors.com](https://deepskycolors.com/tools-tutorials/star-size-reduction-via-morphological-transformations/)
  as "a multiscale algorithm", **page**).
- Size limits (reject large connected structures: nebula knots, galaxy
  cores), then grow (dilate) and feather (Gaussian) the mask.
- The deepskycolors tutorial's settings (threshold 0.15, small scale, low
  growth, smoothness 5, "contours" to protect only star edges) show the
  controls users expect.

### 2.2 Star reduction (shrink)

- **PixInsight MorphologicalTransformation** under a star mask (**page**,
  [deepskycolors.com](https://deepskycolors.com/tools-tutorials/star-size-reduction-via-morphological-transformations/)):
  erosion, or "Morphological Selection", which "acts as a blend between
  the erosion and dilation methods ... more erosion the closer you are to
  a value of 0, and more dilation as you get closer to a value of 1"; a
  round 5 x 5 structuring element; amount and iterations. That is a rank
  (order statistic) filter at quantile q over a disc: q = 0 is erosion
  (min), q = 1 dilation (max), q = 0.5 the median. PixInsight is
  proprietary; nothing of it was read or is used, only this published
  description of the method.
- **Siril does not have a morphological star reduction.** Its
  [Automated Star Reduction Script](https://siril.org/tutorials/star-reduction/)
  (**page**) runs StarNet++ once and blends the starless and original
  images with an MTF pixel math formula (value 0.2 default, 0.5 no
  effect). The task's brief assumed a morphological method in Siril; that
  is not what Siril does. dmoulton/siril-scripts has an MTF-based "[DSA]
  Star Reduction-MTF.ssf" (**page**, needs a starless image too).
- **Noel Carboni's "Make Stars Smaller"** (Photoshop actions, $21.95) is
  a closed action set; its steps were not checked.
- **Hennig's Star Rounding** is the directional variant (§1).
- **Plan for `astro:star-reduce`:** star mask (§2.1) inside the op; rank
  filter at quantile q over a disc of radius r (or erosion then a small
  dilation, an "opening"), repeated n times; blend `out = in + m (rank -
  in) * amount` with the feathered mask m; optionally in luminance only,
  colour kept (the ratio method of the presence ops). All classical
  morphology (Matheron and Serra, 1960s and 1970s).
- **StarNet:** [nekitmm/starnet](https://github.com/nekitmm/starnet) v1:
  `LICENSE.md` is MIT and says "Code is Copyright (c) 2018-2019 Nikita
  Misiura"; the weights are separate Dropbox downloads
  (`wherearemyweights.txt`) with **no licence stated**, so no permission
  to redistribute them. StarNet v2 ([starnetastro.com](https://www.starnetastro.com/))
  shows "© 2026 Mikita Misiura. All rights reserved" and no licence on the
  home page (**page**). Not used, as the brief says.

### 2.3 Background gradient removal

- **Siril Background Extraction** (Siril master `f6a19ca`, 2026-09-28,
  `src/algos/background_extraction.c`, **read**; GPL-3.0-or-later):
  - samples of 25 x 25 px (`SAMPLE_SIZE`), on a grid of `nb_per_line`
    per row (`generate_samples`, line 897) or random in darker areas of
    each quadrant after an auto stretch, with an optional gradient descent
    to the local dark minimum (`generate_samples_random`, "inspired by
    SetiAstro's AutoDBE script ... GPL-v3.0-or-later");
  - per sample the median; a sample is rejected when its median is above
    `median + MAD * tolerance` over all samples (line 966), or contains a
    zero pixel;
  - fit: a polynomial of total degree 1 to 4 (3, 6, 10, 15 terms,
    `poly_1` to `poly_4`) by GSL multifit, or an **RBF**: thin plate
    kernel `r^2 log r`, a smoothing term on the diagonal
    (`1e-4 * 10^((s - 0.5) * 3)` times the mean kernel value) and the
    sample medians as values (`computeBackground_RBF`, line 122);
  - correction by subtraction or division (`BACKGROUND_CORRECTION_SUBTRACT`,
    `_DIVIDE`), per channel.
  - New on master (listed under 1.5.0 in the ChangeLog, so after the
    1.4.4 release): an
    automatic sample-free model: an iterative robust rejection of
    structures (`AG_HIGH_K 2`, `AG_LOW_K 4`, 20 iterations) and a
    multiscale smooth surface, optionally after a stiff polynomial.
- **Hennig's Background Gradient** (§1): boxes, sigma clipping,
  tensor degree 4 polynomial, divide.
- **PixInsight DBE / ABE:** proprietary; not read. Only known from
  tutorials as user or auto placed samples with a tolerance and a spline
  surface.
- **GraXpert:** [Steffenhir/GraXpert](https://github.com/Steffenhir/GraXpert)
  is GPL-3.0 (GitHub licence field, last push 2026-01-12); its AI mode is
  skipped as the brief says (its model files' terms **not checked**). Its
  classical modes (RBF, splines, kriging) are **not checked** here either.
- **Plan for `astro:background`:** grid samples with per cell robust
  statistics (Hennig), MAD rejection of bright cells (Siril), iterative
  rejection of cells with large residuals after each fit, a polynomial of
  total degree 0 to 4 in normalised coordinates, later an RBF (thin plate
  with smoothing), subtract or divide, the model level added back so the
  sky keeps its brightness. Whole-image: §3.

### 2.4 Star colour

- **Siril Saturation** (`src/filters/saturation.c`, **read**): HSL
  saturation with a hue range and a "background factor": pixels with
  lightness below `(median + sigma) * factor` are left alone, so the sky
  noise does not get coloured.
- Photoshop "Increase Star Color" (Carboni) is closed, not checked.
- **Plan for `astro:star-color`:** saturation (chroma in a perceptual or
  the linear ratio space) boosted under the star mask (§2.1), with the
  background threshold of Siril; optionally the colour of a star's halo
  carried into a white, clipped core (cores of bright stars clip to
  white; the ring just outside holds the colour). The halo fill is our
  own idea and needs a test on real data.

### 2.5 Hot pixels

- **Siril Cosmetic Correction** (`src/filters/cosmetic_correction.c`,
  GPL-3.0-or-later, **read**): a pixel is hot (cold) when it deviates from
  the median of its 5 x 5 neighbours (without itself) by more than k sigma
  of the image; a hot pixel gets the mean of its 3 x 3 neighbours, a cold
  one the median; CFA aware (step 2).
- Fits as a single-input area filter (fixed radius 2); the image sigma is
  a global statistic, taken from a subsample as in `astro:background`.
  Mostly a raw frame task that Siril does before stacking; low priority
  for GIMP.

### 2.6 Stretching

- **Arcsinh stretch:** Lupton et al., "Preparing Red-Green-Blue (RGB)
  Images from CCD Data", [arXiv:astro-ph/0312483](https://arxiv.org/abs/astro-ph/0312483)
  (2004): an asinh stretch "allows us to show faint objects while
  simultaneously preserving the structure of brighter objects". Colour is
  kept by scaling all channels by the same factor. Siril's version
  (`src/filters/asinh.c`, **read**): `x' = max(0, (x - offset) / (1 -
  offset))`, factor `k = asinh(beta Y) / (Y asinh(beta))` on the
  luminance Y, all channels times k. One formula, published: trivial to
  write fresh.
- **Generalised Hyperbolic Stretch (GHS):** Dave Payne (equations) and
  Mike Cranfield (code), 2021 to 2022. The equations are published in
  the [GHS documentation](https://www.ghsastro.co.uk/doc/tools/GeneralizedHyperbolicStretch/GeneralizedHyperbolicStretch.html)
  (**page**): `D = e^(stretch factor) - 1`; base transforms by local
  intensity b: `ln(1 + D x)` (b = -1), `(1 - (1 - b D x)^((b+1)/b)) /
  (D (b+1))` (b < 0), `1 - e^(-D x)` (b = 0), `1 - (1 + D x)^-1` (b = 1),
  `1 - (1 + b D x)^(-1/b)` (b > 0); mirrored about the symmetry point SP,
  linear below the shadow protection LP and above the highlight
  protection HP, then normalised. Licences: the PixInsight script
  [mikec1485/GHS](https://github.com/mikec1485/GHS) is GPL-3.0 (GitHub
  licence field; the documentation page says "version 3 of the License",
  without "or later"); Siril's implementation `src/filters/ght.c`
  (**read**, `GHTsetup` line 59, all five cases with precomputed
  coefficients) is GPL-3.0-or-later ("either version 3 ... or (at your
  option) any later version"). So: write it from the published equations,
  or port Siril's `ght.c` (licence compatible with this repository);
  do not copy Cranfield's GPL-3.0-only script into a GPL-3.0-or-later
  repository without asking. Equations themselves are maths, not code.
- **Fit:** both are point operations: ideal NDE GEGL ops. An "auto"
  (black point, midtones balance from statistics, as Siril's autostretch
  `find_unlinked_midtones_balance_default`) needs whole-image statistics
  (§3).

## 3. What fits as a GEGL op in GIMP 3.2

Checked in GIMP 3.2.6 (tag `GIMP_3_2_6`, `e101dd19`, 2026-09-10):

- **Ops with an aux input are always destructive.**
  `app/tools/gimpfiltertool.c` lines 431, 593 and 1583: when the operation
  has an `aux` pad, the "merge filter" option is forced on ("TODO: Once we
  can serialize GimpDrawable, remove so that filters with aux nodes can be
  non-destructive"). So every op here has one input and takes its
  settings as properties; a star mask is computed inside the op, not
  passed in as a second image.
- **Menus:** `gimp:menu-path` is read in `app/menus/filters-menu.c` (line
  ~77); the root must be `<Image>` or another named manager (line 97),
  and a missing submenu is created on demand
  (`app/widgets/gimpmenumodel.c` line 1241, `gimp_menu_model_new_submenu`).
  So `<Image>/Filters/Astronomy` (Hennig's path) works for GEGL ops as it
  does for plug-ins. **Seen working** in GIMP 3.2.6 on Broadway
  (`docs/astronomy-menu.png`, `tests/gui/look.sh`).
- **Property names:** a property named `output` (like the output pad) is
  left out of GIMP's filter dialog, though it can still be set from a
  script (found with the first version of `astro:background`; its
  property is now `show`).
- **Whole-image ops:** a gradient fit needs the whole layer. The house
  pattern (gegl-presence `dehaze.c`): `get_required_for_output`,
  `get_cached_region` and `get_invalidated_by_change` all return the input
  bounding box, `threaded = FALSE`, so GEGL computes the result once and
  caches it, and a changed pixel anywhere recomputes everything.
- **Sample points as properties:** a string property (e.g. `"x,y;x,y"` in
  0 to 1 coordinates) survives the XCF save of a filter; GIMP has no
  multi-point property widget, so user placed samples would come from a
  small plug-in that writes the string (or from GIMP's sample points).
- **Plug-in only:** alignment and stacking (many layers in), anything with
  a second image (flat, dark), and a starless/stars split that makes new
  layers.

## 4. Licences summary

| Source | Licence (from the files) | Use here |
|---|---|---|
| Hennig gimp-plugin-astronomy 0.10 | GPL-2.0-or-later (headers), GPLv2 COPYING | Port under GPL-3.0-or-later, keep © Georg Hennig |
| Siril (master `f6a19ca`) | GPL-3.0-or-later (README; headers of `background_extraction.c`, `ght.c`, `cosmetic_correction.c`, `star_finder.c`; `asinh.c` and `saturation.c` have no header, covered by the project licence) | Portable; keep © team free-astro |
| SetiAstro AutoDBE (via Siril's comment) | GPL-3.0-or-later (per Siril's comment, **not checked** at the source) | Idea only |
| GHS script, Cranfield | GPL-3.0 (only, per its docs page) | Equations only |
| GraXpert | GPL-3.0 (GitHub field) | Not used |
| StarNet v1 code / weights | MIT code / weights without licence | Not used |
| StarNet v2 | "All rights reserved" | Not used |
| PixInsight | proprietary | Published descriptions only |

This repository: **GPL-3.0-or-later** (the same as gegl-presence), which
takes both Hennig's and Siril's code.

## 5. Patents (quick check)

- Google Patents answered HTTP 503 to direct fetches in this session (as
  it did for gegl-presence); the check was done through web searches for
  "star reduction", "star size", "star mask erosion" and "astronomical
  background gradient / light pollution removal" patents. **No patent on
  star reduction or on sample-based sky background extraction for
  astronomical images was found.** Hits were unrelated (star trackers,
  plasma etch mask erosion) or general image background removal:
  US 7,085,413 B2 (image background detection by luminance thresholds),
  US 6,577,762 B1 (background surface thresholding, filed 1999, so
  expired by term), EP 1782621 A1 (background removal for spectra with
  scattering/fluorescence); their status and claims were **not read**
  (503).
- Prior art older than 20 years covers the core methods, so any patent on
  them would have lapsed by term: grey-scale morphology (Matheron,
  Serra, 1960s to 1980s), sigma-clipped sky estimation and polynomial sky
  fits (standard in astronomy for decades), thin plate splines (Duchon
  1977), Hennig's plug-in itself (2006 to 2008), Lupton's asinh (2004).
- RC Astro's pages show no patent notices (search). Not a full
  freedom-to-operate search.
