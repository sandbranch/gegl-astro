/*
 * Shared by the astro operations: the star mask, Gaussian blurs on
 * float planes, robust statistics on a fixed lattice
 *
 * astro-common.h
 * Copyright 2026 by David
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * THE STAR MASK (PLAN.md op 2; docs/research.md 2.1). Classical, no
 * neural net, written for this repository:
 *
 *  1. the luminance L of the image's RGB space (linear light);
 *     transparent and non-finite pixels take the median of L and are
 *     never stars;
 *  2. a difference of Gaussians D = G(s) * L - G(4 s) * L at the star
 *     size s: a smoothing matched to small stars minus the local sky
 *     around them, so gradients, sky glow and nebulae wider than a few
 *     s fall out and stars on a nebula stay separate;
 *  3. the noise of D: its median and robust sigma (1.4826 MAD) on a
 *     fixed lattice of at most 512 x 512 points of the whole layer. Near
 *     the edges (mirrored) D is noisier, by a factor that is exact for
 *     white noise and computed from the two kernels; D is divided by it
 *     first, so the threshold holds up to the edges;
 *  4. pixels above median + threshold x sigma, in 8-connected regions;
 *     a region is a star when its full width at half maximum is at most
 *     max-size: the area, inside the region, where the light above the
 *     local sky is at least
 *     half its peak, taken as a disc. Nebula knots and galaxy cores are
 *     left out; bright stars with a wide footprint are kept (their half
 *     maximum is still small). Objects wider than max-size are read a
 *     little smaller than they are (the blurred sky takes some of their
 *     light): by the formulas, a Gaussian of FWHM 1.1 x max-size reads
 *     about 10 % smaller;
 *  5. the stars' pixels grown by a disc of radius grow, then softened
 *     by a Gaussian of sigma feather: the mask, 0 to 1.
 *
 * The local sky: L blurred with a sigma of max-size, L clipped at that
 * plus 3 robust sigmas of the difference (so star light hardly counts),
 * blurred again.
 *
 * This follows the published outline of PixInsight's StarMask
 * (multiscale detection, threshold, growth, smoothness; deepskycolors.com)
 * without its code, and the classical detection of SExtractor (Bertin
 * and Arnouts 1996: a matched filter, a threshold in noise sigmas,
 * connected pixels).
 *
 * The whole layer is needed (a region's size is only known when all of
 * it is seen), so the operations that use the mask ask for and cache
 * their whole input. Everything is relative to the layer's corner, so
 * the result does not depend on where the layer sits. Deterministic.
 */

#ifndef ASTRO_COMMON_H
#define ASTRO_COMMON_H

#include <gegl.h>
#include <gegl-plugin.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define ASTRO_MAD_TO_SIGMA 1.4826
#define ASTRO_STATS_SIDE   512
#define ASTRO_DOG_RATIO    4.0
#define ASTRO_STRIP        128

/* The star mask settings, the same in every operation that makes a
 * mask (used inside the GEGL_PROPERTIES part of an operation). */
#define ASTRO_STAR_MASK_PROPERTIES                                            \
  property_double (threshold, _("Detection threshold"), 4.0)                \
    description (_("How far above the noise a star must rise, in noise "    \
                   "sigmas. Lower finds fainter stars and more noise."))    \
    value_range (1.0, 100.0)                                                 \
    ui_range (2.0, 20.0)                                                     \
    ui_digits (1)                                                            \
  property_double (scale, _("Star size"), 1.5)                              \
    description (_("Size of the smallest stars, in pixels (the sigma of "   \
                   "the detection filter)."))                                \
    value_range (0.5, 10.0)                                                  \
    ui_digits (2)                                                            \
  property_double (max_size, _("Largest star (FWHM)"), 15.0)               \
    description (_("Structures wider than this at half their peak, in "     \
                   "pixels, are not stars (nebula knots, galaxy cores)."))  \
    value_range (2.0, 100.0)                                                 \
    ui_range (4.0, 50.0)                                                     \
    ui_digits (0)                                                            \
  property_double (grow, _("Grow"), 1.0)                                    \
    description (_("Widens the mask around each star, in pixels."))         \
    value_range (0.0, 20.0)                                                  \
    ui_digits (1)                                                            \
  property_double (feather, _("Feather"), 1.5)                              \
    description (_("Softens the mask's edge, in pixels (sigma)."))          \
    value_range (0.0, 20.0)                                                  \
    ui_digits (1)

typedef struct
{
  gdouble threshold, scale, max_size, grow, feather;
} AstroStarMaskParams;

#define ASTRO_STAR_MASK_PARAMS(o) \
  { (o)->threshold, (o)->scale, (o)->max_size, (o)->grow, (o)->feather }

/* ---- small helpers ------------------------------------------------------ */

static inline void
astro_luminance_coefficients (const Babl *format,
                              gdouble     lum[3])
{
  lum[0] = 0.2126;
  lum[1] = 0.7152;
  lum[2] = 0.0722;
  babl_space_get_rgb_luminance (babl_format_get_space (format),
                                &lum[0], &lum[1], &lum[2]);
}

static int
astro_cmp_float (const void *a,
                 const void *b)
{
  gfloat x = *(const gfloat *) a, y = *(const gfloat *) b;

  return (x > y) - (x < y);
}

/* the k-th smallest of v[0..n-1] (reordered; quickselect, deterministic) */
static inline gfloat
astro_select (gfloat *v,
              gint    n,
              gint    k)
{
  gint lo = 0, hi = n - 1;

  while (lo < hi)
    {
      gfloat pivot = v[(lo + hi) / 2];
      gint   i = lo, j = hi;

      while (i <= j)
        {
          while (v[i] < pivot)
            i++;
          while (v[j] > pivot)
            j--;
          if (i <= j)
            {
              gfloat t = v[i];

              v[i] = v[j];
              v[j] = t;
              i++;
              j--;
            }
        }
      if (k <= j)
        hi = j;
      else if (k >= i)
        lo = i;
      else
        break;
    }
  return v[k];
}

/* median and robust sigma of v[0..n-1] (reordered) */
static void
astro_robust_stats (gfloat  *v,
                    gint     n,
                    gdouble *median,
                    gdouble *sigma)
{
  gfloat med;
  gint   i;

  if (n < 1)
    {
      *median = 0.0;
      *sigma  = 0.0;
      return;
    }
  qsort (v, n, sizeof (gfloat), astro_cmp_float);
  med = n % 2 ? v[n / 2] : 0.5f * (v[n / 2 - 1] + v[n / 2]);
  for (i = 0; i < n; i++)
    v[i] = fabsf (v[i] - med);
  qsort (v, n, sizeof (gfloat), astro_cmp_float);
  *median = med;
  *sigma  = ASTRO_MAD_TO_SIGMA * (n % 2 ? v[n / 2] : 0.5f * (v[n / 2 - 1] + v[n / 2]));
}

/* the lattice step for at most side x side points of a w x h plane */
static inline gint
astro_lattice_step (gint w,
                    gint h,
                    gint side)
{
  return MAX (1, (MAX (w, h) + side - 1) / side);
}

/* ---- Gaussian blur of a float plane, in place ----------------------------- */

typedef struct
{
  gfloat  *plane;
  gint     w, h;
  gfloat  *kernel;   /* 2 r + 1 taps */
  gint     r;
} AstroBlur;

static gfloat *
astro_gauss_kernel (gdouble sigma,
                    gint   *radius)
{
  gint    r = MAX (1, (gint) ceil (4.0 * sigma)), i;
  gfloat *k = g_new (gfloat, 2 * r + 1);
  gdouble sum = 0.0;

  for (i = -r; i <= r; i++)
    sum += k[i + r] = (gfloat) exp (-0.5 * i * i / (sigma * sigma));
  for (i = 0; i <= 2 * r; i++)
    k[i] /= sum;
  *radius = r;
  return k;
}

/* the index of x mirrored into 0 to n - 1 (edges repeated once, as
 * ... 2 1 0 | 0 1 2 ...) */
static inline gint
astro_mirror (gint x,
              gint n)
{
  if (n == 1)
    return 0;
  while (x < 0 || x >= n)
    x = x < 0 ? -x - 1 : 2 * n - x - 1;
  return x;
}

static void
astro_blur_rows (gsize    offset,
                 gsize    size,
                 gpointer data)
{
  AstroBlur *b = data;
  gint       w = b->w, r = b->r;
  gfloat    *tmp = g_new (gfloat, w + 2 * r);
  gsize      y;

  for (y = offset; y < offset + size; y++)
    {
      gfloat *row = b->plane + y * w;
      gint    x, j;

      for (x = -r; x < w + r; x++)
        tmp[x + r] = row[astro_mirror (x, w)];
      for (x = 0; x < w; x++)
        {
          const gfloat *t = tmp + x;
          gfloat        s = 0.0f;

          for (j = 0; j <= 2 * r; j++)
            s += b->kernel[j] * t[j];
          row[x] = s;
        }
    }
  g_free (tmp);
}

#define ASTRO_COLS 64

/* columns in strips of ASTRO_COLS, in place with a ring of 2 r + 1 rows;
 * the last r + 1 rows are kept first, for the mirrored rows past the end
 * (read when the rows they mirror are already written) */
static void
astro_blur_cols (gsize    offset,
                 gsize    size,
                 gpointer data)
{
  AstroBlur *b = data;
  gint       w = b->w, h = b->h, r = b->r, n = 2 * r + 1;
  gint       t0 = MAX (0, h - r - 1);
  gfloat    *ring = g_new (gfloat, (gsize) n * ASTRO_COLS);
  gfloat    *tail = g_new (gfloat, (gsize) (h - t0) * ASTRO_COLS);
  gsize      s;

  for (s = offset; s < offset + size; s++)
    {
      gint x0 = (gint) s * ASTRO_COLS, x1 = MIN (w, x0 + ASTRO_COLS), cw = x1 - x0;
      gint y, j, x;

      for (y = t0; y < h; y++)
        memcpy (tail + (gsize) (y - t0) * ASTRO_COLS, b->plane + (gsize) y * w + x0,
                cw * sizeof (gfloat));
      for (j = -r; j <= r; j++)
        {
          gint src = astro_mirror (j, h);

          memcpy (ring + (gsize) (((j % n) + n) % n) * ASTRO_COLS,
                  src >= t0 ? tail + (gsize) (src - t0) * ASTRO_COLS
                            : b->plane + (gsize) src * w + x0,
                  cw * sizeof (gfloat));
        }
      for (y = 0; y < h; y++)
        {
          gfloat *out = b->plane + (gsize) y * w + x0;
          gfloat  acc[ASTRO_COLS];
          gint    src;

          memset (acc, 0, sizeof acc);
          for (j = -r; j <= r; j++)
            {
              const gfloat *row = ring + (gsize) ((((y + j) % n) + n) % n) * ASTRO_COLS;
              gfloat        k = b->kernel[j + r];

              for (x = 0; x < cw; x++)
                acc[x] += k * row[x];
            }
          memcpy (out, acc, cw * sizeof (gfloat));
          /* row y - r leaves the ring, row y + r + 1 comes in: from the
           * plane while it is not written yet, else from the tail */
          j   = y + r + 1;
          src = astro_mirror (j, h);
          memcpy (ring + (gsize) (((j % n) + n) % n) * ASTRO_COLS,
                  src >= t0 ? tail + (gsize) (src - t0) * ASTRO_COLS
                            : b->plane + (gsize) src * w + x0,
                  cw * sizeof (gfloat));
        }
    }
  g_free (tail);
  g_free (ring);
}

/* blurs the w x h plane in place with a Gaussian of sigma (edges
 * mirrored); threads by GEGL; the same result for any number of threads */
static void
astro_gauss_blur (gfloat *plane,
                  gint    w,
                  gint    h,
                  gdouble sigma)
{
  AstroBlur b;

  if (sigma < 0.05 || w < 1 || h < 1)
    return;
  b.plane  = plane;
  b.w      = w;
  b.h      = h;
  b.kernel = astro_gauss_kernel (sigma, &b.r);
  gegl_parallel_distribute_range (h, 16, astro_blur_rows, &b);
  gegl_parallel_distribute_range ((w + ASTRO_COLS - 1) / ASTRO_COLS, 1, astro_blur_cols, &b);
  g_free (b.kernel);
}

/* ---- the noise of a difference of Gaussians near the edges --------------- */

/* For white noise, the variance of G(s) * L - G(t) * L with mirrored
 * edges at position x of a line of n: A = sum of the squared effective
 * weights of the s kernel, B of the t kernel, C of their products (the
 * 2D variance is A_x A_y - 2 C_x C_y + B_x B_y, as both kernels are
 * separable). Near an edge a mirrored pixel counts twice, so the
 * variance goes up there (up to about twice). */
static void
astro_dog_noise_1d (gint          n,
                    const gfloat *ks,
                    gint          rs,
                    const gfloat *kt,
                    gint          rt,
                    gdouble      *A,
                    gdouble      *B,
                    gdouble      *C)
{
  gdouble *ea = g_new0 (gdouble, n), *eb = g_new0 (gdouble, n);
  gint     x, j;

  for (x = 0; x < n; x++)
    {
      gint lo = n, hi = -1, m;

      for (j = -rt; j <= rt; j++)
        {
          m = astro_mirror (x + j, n);
          eb[m] += kt[j + rt];
          if (j >= -rs && j <= rs)
            ea[m] += ks[j + rs];
          lo = MIN (lo, m);
          hi = MAX (hi, m);
        }
      A[x] = B[x] = C[x] = 0.0;
      for (m = lo; m <= hi; m++)
        {
          A[x] += ea[m] * ea[m];
          B[x] += eb[m] * eb[m];
          C[x] += ea[m] * eb[m];
          ea[m] = eb[m] = 0.0;
        }
    }
  g_free (ea);
  g_free (eb);
}

/* divides the difference of Gaussians by its noise's standard deviation
 * relative to the middle of the plane (white noise), so one threshold
 * holds up to the edges */
static void
astro_dog_normalise (gfloat *dog,
                     gint    w,
                     gint    h,
                     gdouble s,
                     gdouble t)
{
  gint     rs, rt, x, y;
  gfloat  *ks = astro_gauss_kernel (s, &rs), *kt = astro_gauss_kernel (t, &rt);
  gdouble *ax = g_new (gdouble, w), *bx = g_new (gdouble, w), *cx = g_new (gdouble, w);
  gdouble *ay = g_new (gdouble, h), *by = g_new (gdouble, h), *cy = g_new (gdouble, h);
  gdouble  ai = 0, bi = 0, ci = 0, mid;

  astro_dog_noise_1d (w, ks, rs, kt, rt, ax, bx, cx);
  astro_dog_noise_1d (h, ks, rs, kt, rt, ay, by, cy);
  /* the middle of an endless line */
  for (x = -rt; x <= rt; x++)
    {
      gdouble a = (x >= -rs && x <= rs) ? ks[x + rs] : 0.0, b = kt[x + rt];

      ai += a * a;
      bi += b * b;
      ci += a * b;
    }
  mid = ai * ai - 2 * ci * ci + bi * bi;
  for (y = 0; y < h; y++)
    for (x = 0; x < w; x++)
      {
        gdouble v = ax[x] * ay[y] - 2 * cx[x] * cy[y] + bx[x] * by[y];

        dog[(gsize) y * w + x] *= (gfloat) sqrt (mid / MAX (v, 1e-30));
      }
  g_free (ks);
  g_free (kt);
  g_free (ax);
  g_free (bx);
  g_free (cx);
  g_free (ay);
  g_free (by);
  g_free (cy);
}

/* ---- the star mask ---------------------------------------------------------- */

/* the luminance of the whole box into a new plane; transparent and
 * non-finite pixels get NAN */
static gfloat *
astro_read_luminance (GeglBuffer          *input,
                      const GeglRectangle *box,
                      const Babl          *format)
{
  gint    w = box->width, h = box->height, y0;
  gfloat *lum = g_new (gfloat, (gsize) w * h);
  gfloat *strip = g_new (gfloat, (gsize) w * ASTRO_STRIP * 4);
  gdouble k[3];

  astro_luminance_coefficients (format, k);
  for (y0 = 0; y0 < h; y0 += ASTRO_STRIP)
    {
      gint   sh = MIN (ASTRO_STRIP, h - y0);
      gsize  i, n = (gsize) w * sh;
      gfloat *out = lum + (gsize) y0 * w;

      gegl_buffer_get (input, GEGL_RECTANGLE (box->x, box->y + y0, w, sh), 1.0,
                       format, strip, GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);
      for (i = 0; i < n; i++)
        {
          const gfloat *p = strip + i * 4;

          if (p[3] > 0.5f && isfinite (p[0]) && isfinite (p[1]) && isfinite (p[2]))
            out[i] = (gfloat) (k[0] * p[0] + k[1] * p[1] + k[2] * p[2]);
          else
            out[i] = NAN;
        }
    }
  g_free (strip);
  return lum;
}

/* median and robust sigma of the finite values of a plane, on a fixed
 * lattice; FALSE when there are too few */
static gboolean
astro_plane_stats (const gfloat *plane,
                   const guint8 *skip,
                   gint          w,
                   gint          h,
                   gdouble      *median,
                   gdouble      *sigma)
{
  gint    step = astro_lattice_step (w, h, ASTRO_STATS_SIDE);
  gint    nx = (w + step - 1) / step, ny = (h + step - 1) / step;
  gfloat *v = g_new (gfloat, (gsize) nx * ny);
  gint    x, y, n = 0;

  for (y = step / 2; y < h; y += step)
    for (x = step / 2; x < w; x += step)
      {
        gsize i = (gsize) y * w + x;

        if (isfinite (plane[i]) && ! (skip && skip[i]))
          v[n++] = plane[i];
      }
  astro_robust_stats (v, n, median, sigma);
  g_free (v);
  return n >= 16;
}

typedef struct
{
  gint    *stack;
  gsize    n, size;
} AstroStack;

static inline void
astro_push (AstroStack *s,
            gint        i)
{
  if (s->n == s->size)
    {
      s->size = MAX (1024, s->size * 2);
      s->stack = g_renew (gint, s->stack, s->size);
    }
  s->stack[s->n++] = i;
}

/* Makes the star mask of the whole box: w x h values from 0 to 1, in a
 * new plane (row by row from the box's corner). n_stars (may be NULL)
 * gets the number of stars found; sky (may be NULL) a new plane with the
 * local sky's luminance (the luminance blurred at max-size), which the
 * caller frees. */
static gfloat *
astro_star_mask (GeglBuffer                *input,
                 const GeglRectangle       *box,
                 const Babl                *format,
                 const AstroStarMaskParams *p,
                 gint                      *n_stars,
                 gfloat                   **sky)
{
  gint        w = box->width, h = box->height;
  gsize       n = (gsize) w * h, i;
  gfloat     *lum = astro_read_luminance (input, box, format);
  gfloat     *dog, *tmp;
  guint8     *skip = g_new0 (guint8, n);
  gint32     *label;
  gdouble     med, sigma, thr;
  AstroStack  st = { NULL, 0, 0 };
  gint        stars = 0;

  if (sky)
    *sky = NULL;
  /* 1. transparent pixels: the median luminance, never stars */
  for (i = 0; i < n; i++)
    skip[i] = ! isfinite (lum[i]);
  if (! astro_plane_stats (lum, skip, w, h, &med, &sigma))
    {
      /* (nearly) nothing opaque: no stars */
      g_free (skip);
      memset (lum, 0, n * sizeof (gfloat));
      if (sky)
        *sky = g_memdup2 (lum, n * sizeof (gfloat));
      if (n_stars)
        *n_stars = 0;
      return lum;
    }
  for (i = 0; i < n; i++)
    if (skip[i])
      lum[i] = (gfloat) med;

  /* 2. the difference of Gaussians */
  dog = g_memdup2 (lum, n * sizeof (gfloat));
  astro_gauss_blur (dog, w, h, p->scale);
  tmp = g_memdup2 (lum, n * sizeof (gfloat));
  astro_gauss_blur (tmp, w, h, ASTRO_DOG_RATIO * p->scale);
  for (i = 0; i < n; i++)
    dog[i] -= tmp[i];
  astro_dog_normalise (dog, w, h, p->scale, ASTRO_DOG_RATIO * p->scale);

  /* the sky at max-size, without the stars: the luminance blurred, the
   * luminance clipped at that plus 3 noise sigmas, blurred again; and
   * the luminance above it (for the size) */
  memcpy (tmp, lum, n * sizeof (gfloat));
  astro_gauss_blur (tmp, w, h, p->max_size);
  {
    gint    step = astro_lattice_step (w, h, ASTRO_STATS_SIDE), x, y, m = 0;
    gfloat *v = g_new (gfloat, (gsize) ((w + step - 1) / step) * ((h + step - 1) / step));
    gdouble lmed, lsig;

    for (y = step / 2; y < h; y += step)
      for (x = step / 2; x < w; x += step)
        v[m++] = lum[(gsize) y * w + x] - tmp[(gsize) y * w + x];
    astro_robust_stats (v, m, &lmed, &lsig);
    g_free (v);
    for (i = 0; i < n; i++)
      tmp[i] = MIN (lum[i], tmp[i] + (gfloat) (lmed + 3.0 * lsig));
  }
  astro_gauss_blur (tmp, w, h, p->max_size);
  for (i = 0; i < n; i++)
    lum[i] -= tmp[i];

  /* 3. the noise of the difference of Gaussians */
  astro_plane_stats (dog, skip, w, h, &med, &sigma);
  thr = med + p->threshold * MAX (sigma, 1e-12);

  /* 4. connected regions above the threshold (labels in dog's memory:
   * 0 below, -1 above and not yet seen, then 1 for stars, 2 for others) */
  label = (gint32 *) dog;
  for (i = 0; i < n; i++)
    label[i] = (! skip[i] && dog[i] > thr) ? -1 : 0;
  g_free (skip);
  for (i = 0; i < n; i++)
    if (label[i] == -1)
      {
        gdouble peak = -G_MAXDOUBLE, fwhm;
        gsize   k, half = 0;

        st.n = 0;
        label[i] = -2;
        astro_push (&st, (gint) i);
        for (k = 0; k < st.n; k++)
          {
            gint j = st.stack[k], x = j % w, y = j / w, dx, dy;

            peak = MAX (peak, lum[j]);
            for (dy = -1; dy <= 1; dy++)
              for (dx = -1; dx <= 1; dx++)
                {
                  gint xx = x + dx, yy = y + dy, jj;

                  if (xx < 0 || yy < 0 || xx >= w || yy >= h)
                    continue;
                  jj = yy * w + xx;
                  if (label[jj] == -1)
                    {
                      label[jj] = -2;
                      astro_push (&st, jj);
                    }
                }
          }
        /* the size: the full width at half maximum of the light above
         * the sky, from the area at or above half the peak */
        for (k = 0; k < st.n; k++)
          half += lum[st.stack[k]] >= 0.5 * peak;
        fwhm = 2.0 * sqrt (half / G_PI);
        {
          gint32 kind = fwhm <= p->max_size ? 1 : 2;

          for (k = 0; k < st.n; k++)
            label[st.stack[k]] = kind;
          stars += kind == 1;
        }
      }
  g_free (st.stack);

  /* 5. the stars' pixels, grown by a disc (into lum's memory) */
  memset (lum, 0, n * sizeof (gfloat));
  {
    gint g = (gint) floor (p->grow), x, y, dx, dy;
    gdouble g2 = p->grow * p->grow + 1e-9;

    for (y = 0; y < h; y++)
      for (x = 0; x < w; x++)
        if (label[(gsize) y * w + x] == 1)
          {
            gboolean edge = FALSE;

            lum[(gsize) y * w + x] = 1.0f;
            if (g < 1)
              continue;
            /* only the region's edge pixels stamp the disc */
            for (dy = -1; dy <= 1 && ! edge; dy++)
              for (dx = -1; dx <= 1; dx++)
                {
                  gint xx = x + dx, yy = y + dy;

                  if (xx >= 0 && yy >= 0 && xx < w && yy < h &&
                      label[(gsize) yy * w + xx] != 1)
                    edge = TRUE;
                }
            if (! edge)
              continue;
            for (dy = -g; dy <= g; dy++)
              for (dx = -g; dx <= g; dx++)
                {
                  gint xx = x + dx, yy = y + dy;

                  if (dx * dx + dy * dy <= g2 && xx >= 0 && yy >= 0 && xx < w && yy < h)
                    lum[(gsize) yy * w + xx] = 1.0f;
                }
          }
  }
  g_free (dog);

  /* feathered */
  astro_gauss_blur (lum, w, h, p->feather);
  for (i = 0; i < n; i++)
    lum[i] = CLAMP (lum[i], 0.0f, 1.0f);
  if (sky)
    *sky = tmp;
  else
    g_free (tmp);
  if (n_stars)
    *n_stars = stars;
  return lum;
}

/* ---- whole-image plumbing ------------------------------------------------- */

static inline GeglRectangle
astro_whole_input (GeglOperation       *operation,
                   const GeglRectangle *roi)
{
  const GeglRectangle *in = gegl_operation_source_get_bounding_box (operation, "input");

  return in && ! gegl_rectangle_is_infinite_plane (in) ? *in : *roi;
}

#endif
