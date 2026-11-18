/*
 * Background: removes the sky gradient (light pollution, sky glow,
 * vignetting) of an astronomical image, as a GEGL operation
 *
 * background.c
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
 * FIRST VERSION (docs/research.md, PLAN.md): a polynomial surface only.
 *
 * The image is cut into a grid of square cells. Per cell and channel the
 * sky level is a sigma-clipped median of a regular subsample of the cell
 * (every step-th pixel, so that a cell gives at most 64 x 64 values):
 * stars are a small part of a cell and fall outside the clip. This is the
 * sampling of Georg Hennig's Background Gradient plug-in for GIMP 2 (boxes
 * with sigma clipping) without its tensor polynomial. Cells brighter than
 * the median of all cells by more than "tolerance" robust sigmas
 * (1.4826 x MAD) are dropped, as Siril does with its samples (nebulae,
 * galaxies, bright stars); the default 5 (2 before 2026-11-18) keeps a
 * strong light dome's own cells (PLAN.md, the tolerance study). A polynomial of total degree 0 to 4 in
 * coordinates scaled to -1 to 1 is fitted to the cells by least squares,
 * and cells more than 3 robust sigmas of the residuals off the surface
 * are dropped and the surface fitted again, up to 4 times.
 *
 * The result is the image minus the surface (or divided by it) plus the
 * surface's mean level, so the sky keeps its brightness; or the surface
 * itself, to inspect it. The whole layer is needed for the fit: the
 * operation asks for and caches its whole input and output. Deterministic:
 * the same pixels give the same result, bit for bit, wherever the layer
 * starts.
 *
 * "Leave out stars": the star mask of astro-common.h at its default
 * settings first; pixels under the mask (> 0) are not sampled, for dense
 * star fields where the stars fill so much of a cell that its clipped
 * median rises.
 */

#include <glib/gi18n-lib.h>
#include "astro-common.h"

#ifdef GEGL_PROPERTIES

enum_start (astro_background_correction)
  enum_value (ASTRO_BACKGROUND_SUBTRACT, "subtract",
              N_("Subtract (light pollution, sky glow)"))
  enum_value (ASTRO_BACKGROUND_DIVIDE,   "divide",
              N_("Divide (vignetting)"))
enum_end (AstroBackgroundCorrection)

enum_start (astro_background_show)
  enum_value (ASTRO_BACKGROUND_CORRECTED, "corrected", N_("Corrected image"))
  enum_value (ASTRO_BACKGROUND_MODEL,     "model",     N_("Background model"))
enum_end (AstroBackgroundShow)

property_enum (correction, _("Correction"),
               AstroBackgroundCorrection, astro_background_correction,
               ASTRO_BACKGROUND_SUBTRACT)
  description (_("Subtract the background (sky glow and light pollution "
                 "add light) or divide by it (vignetting takes a share "
                 "of the light)."))

property_int (degree, _("Degree"), 2)
  description (_("Degree of the polynomial surface: 0 a constant level, "
                 "1 a plane (a linear gradient), 2 to 4 curved gradients. "
                 "Higher degrees follow the sky more closely but may bend "
                 "into large nebulae."))
  value_range (0, 4)

property_int (samples, _("Samples per side"), 16)
  description (_("Number of sample cells along the longer side of the "
                 "image; the cells are square."))
  value_range (4, 64)

property_double (tolerance, _("Tolerance"), 5.0)
  description (_("Cells brighter than the typical cell by more than this "
                 "many robust standard deviations are left out of the fit "
                 "(nebulae, galaxies, bright stars). Lower leaves out "
                 "more: try 2 when a large galaxy or nebula bends the "
                 "background; keep it high for strong light pollution, "
                 "whose brightest part would be left out otherwise."))
  value_range (0.2, 10.0)
  ui_digits (1)

property_boolean (exclude_stars, _("Leave out stars"), FALSE)
  description (_("Leaves the stars out of each cell's sample, with the "
                 "star mask of Star Mask at its default settings: for "
                 "dense star fields (the Milky Way), where stars raise "
                 "the sky level of a cell. Slower."))

property_boolean (keep_level, _("Keep the sky level"), TRUE)
  description (_("Adds the mean level of the background back after "
                 "subtracting (or multiplies by it after dividing), so "
                 "the sky keeps its brightness; off: the sky goes to "
                 "black (or to 1)."))

/* (not "output": that is the name of the output pad, and GIMP's dialog
 * leaves such a property out) */
property_enum (show, _("Show"),
               AstroBackgroundShow, astro_background_show,
               ASTRO_BACKGROUND_CORRECTED)
  description (_("The corrected image, or the fitted background itself."))

#else

#define GEGL_OP_FILTER
#define GEGL_OP_NAME     background
#define GEGL_OP_C_SOURCE background.c

#include "gegl-op.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAX_TERMS      15   /* total degree 4: (4 + 1) (4 + 2) / 2 */
#define MAX_SIDE       64   /* at most 64 x 64 values per cell and channel */
#define CLIP_SIGMAS    3.0
#define CLIP_ROUNDS    3
#define FIT_ROUNDS     4
#define RESIDUAL_SIGMA 3.0
#define MAD_TO_SIGMA   1.4826

typedef struct
{
  gdouble u, v;       /* centre, -1 to 1 */
  gdouble level[3];   /* sky level per channel */
  gdouble lum;        /* mean of the three */
  gboolean use;
} Cell;

static int
cmp_double (const void *a,
            const void *b)
{
  gdouble x = *(const gdouble *) a, y = *(const gdouble *) b;

  return (x > y) - (x < y);
}

/* median of v[0..n-1], reordered (quickselect of astro-common.h; the
 * same values as a sort would give) */
static gfloat
median_float (gfloat *v,
              gint    n)
{
  gfloat hi = astro_select (v, n, n / 2), lo;
  gint   i;

  if (n % 2)
    return hi;
  /* the largest of the lower half */
  lo = v[0];
  for (i = 1; i < n / 2; i++)
    lo = MAX (lo, v[i]);
  return 0.5f * (lo + hi);
}

static gdouble
median_double (gdouble *v,
               gint     n)
{
  qsort (v, n, sizeof (gdouble), cmp_double);
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

/* sigma-clipped median of v[0..n-1] (reordered); tmp holds n values */
static gdouble
clipped_median (gfloat *v,
                gint    n,
                gfloat *tmp)
{
  gfloat med = median_float (v, n);
  gint   round, i;

  for (round = 0; round < CLIP_ROUNDS && n > 4; round++)
    {
      gfloat sigma, lo, hi;
      gint   m = 0;

      for (i = 0; i < n; i++)
        tmp[i] = fabsf (v[i] - med);
      sigma = MAD_TO_SIGMA * median_float (tmp, n);
      if (sigma <= 0.0f)
        break;
      lo = med - CLIP_SIGMAS * sigma;
      hi = med + CLIP_SIGMAS * sigma;
      for (i = 0; i < n; i++)
        if (v[i] >= lo && v[i] <= hi)
          v[m++] = v[i];
      if (m == n || m < 4)
        break;
      n   = m;
      med = median_float (v, n);
    }
  return med;
}

static gint
n_terms (gint degree)
{
  return (degree + 1) * (degree + 2) / 2;
}

/* the terms u^i v^j, i + j <= degree, in the order of total degree */
static void
terms (gdouble  u,
       gdouble  v,
       gint     degree,
       gdouble *t)
{
  gdouble up[5], vp[5];
  gint    d, i, k = 0;

  up[0] = vp[0] = 1.0;
  for (i = 1; i <= 4; i++)
    {
      up[i] = up[i - 1] * u;
      vp[i] = vp[i - 1] * v;
    }
  for (d = 0; d <= degree; d++)
    for (i = d; i >= 0; i--)
      t[k++] = up[i] * vp[d - i];
}

/* solves the n x n system a x = b in place (Gaussian elimination with
 * partial pivoting); FALSE when it is singular */
static gboolean
solve (gdouble *a,
       gdouble *b,
       gint     n)
{
  gint col, row, k;

  for (col = 0; col < n; col++)
    {
      gint    piv = col;
      gdouble best = fabs (a[col * n + col]);

      for (row = col + 1; row < n; row++)
        if (fabs (a[row * n + col]) > best)
          {
            best = fabs (a[row * n + col]);
            piv  = row;
          }
      if (best < 1e-12)
        return FALSE;
      if (piv != col)
        {
          for (k = 0; k < n; k++)
            {
              gdouble t = a[col * n + k];

              a[col * n + k] = a[piv * n + k];
              a[piv * n + k] = t;
            }
          {
            gdouble t = b[col];

            b[col] = b[piv];
            b[piv] = t;
          }
        }
      for (row = col + 1; row < n; row++)
        {
          gdouble f = a[row * n + col] / a[col * n + col];

          for (k = col; k < n; k++)
            a[row * n + k] -= f * a[col * n + k];
          b[row] -= f * b[col];
        }
    }
  for (row = n - 1; row >= 0; row--)
    {
      gdouble s = b[row];

      for (k = row + 1; k < n; k++)
        s -= a[row * n + k] * b[k];
      b[row] = s / a[row * n + row];
    }
  return TRUE;
}

/* least squares fit of the used cells, per channel, into coef[3][MAX_TERMS];
 * FALSE when the system is singular */
static gboolean
fit (const Cell *cells,
     gint        n_cells,
     gint        degree,
     gdouble     coef[3][MAX_TERMS])
{
  gint    nt = n_terms (degree);
  gdouble ata[MAX_TERMS * MAX_TERMS];
  gdouble atb[3][MAX_TERMS];
  gint    i, j, k, c;

  memset (ata, 0, sizeof ata);
  memset (atb, 0, sizeof atb);
  for (i = 0; i < n_cells; i++)
    {
      gdouble t[MAX_TERMS];

      if (! cells[i].use)
        continue;
      terms (cells[i].u, cells[i].v, degree, t);
      for (j = 0; j < nt; j++)
        {
          for (k = 0; k < nt; k++)
            ata[j * nt + k] += t[j] * t[k];
          for (c = 0; c < 3; c++)
            atb[c][j] += t[j] * cells[i].level[c];
        }
    }
  for (c = 0; c < 3; c++)
    {
      gdouble a[MAX_TERMS * MAX_TERMS];

      memcpy (a, ata, sizeof a);
      if (! solve (a, atb[c], nt))
        return FALSE;
      memcpy (coef[c], atb[c], nt * sizeof (gdouble));
    }
  return TRUE;
}

static gdouble
eval (const gdouble *coef,
      gint           degree,
      gdouble        u,
      gdouble        v)
{
  gdouble t[MAX_TERMS], s = 0.0;
  gint    k, nt = n_terms (degree);

  terms (u, v, degree, t);
  for (k = 0; k < nt; k++)
    s += coef[k] * t[k];
  return s;
}

/* measures the cells of the grid (a cell with too little opaque sky
 * is not used) */
static Cell *
measure (GeglBuffer          *input,
         const GeglRectangle *box,
         const Babl          *format,
         gint                 samples,
         const gfloat        *stars,
         gint                *n_cells_out)
{
  gint    side = MAX (box->width, box->height);
  gdouble size = (gdouble) side / samples;
  gint    nx = MAX (1, (gint) ceil (box->width / size - 1e-9));
  gint    ny = MAX (1, (gint) ceil (box->height / size - 1e-9));
  gint    csize = (gint) ceil (size);
  Cell   *cells = g_new0 (Cell, nx * ny);
  gfloat *pix = g_new (gfloat, (gsize) csize * 4);
  gfloat *vals[3], *tmp;
  gint    cx, cy, c;

  for (c = 0; c < 3; c++)
    vals[c] = g_new (gfloat, MAX_SIDE * MAX_SIDE);
  tmp = g_new (gfloat, MAX_SIDE * MAX_SIDE);

  for (cy = 0; cy < ny; cy++)
    for (cx = 0; cx < nx; cx++)
      {
        Cell         *cell = &cells[cy * nx + cx];
        gint          x0 = (gint) floor (cx * size), x1 = (gint) floor ((cx + 1) * size);
        gint          y0 = (gint) floor (cy * size), y1 = (gint) floor ((cy + 1) * size);
        GeglRectangle r;
        gint          step, x, y, n = 0, total = 0, opaque = 0;

        x1 = MIN (x1, box->width);
        y1 = MIN (y1, box->height);
        if (cx == nx - 1)
          x1 = box->width;
        if (cy == ny - 1)
          y1 = box->height;
        cell->use = FALSE;
        if (x1 <= x0 || y1 <= y0)
          continue;
        gegl_rectangle_set (&r, box->x + x0, box->y + y0, x1 - x0, y1 - y0);
        step = MAX (1, (MAX (r.width, r.height) + MAX_SIDE - 1) / MAX_SIDE);
        /* only the lattice's rows are read */
        for (y = step / 2; y < r.height; y += step)
          {
            gegl_buffer_get (input, GEGL_RECTANGLE (r.x, r.y + y, r.width, 1), 1.0, format,
                             pix, GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);
            for (x = step / 2; x < r.width; x += step)
              {
                const gfloat *p = pix + (gsize) x * 4;

                total++;
                if (! (p[3] > 0.5f) || ! isfinite (p[0]) || ! isfinite (p[1]) ||
                    ! isfinite (p[2]))
                  continue;
                opaque++;
                if (stars && stars[(gsize) (y0 + y) * box->width + x0 + x] > 0.0f)
                  continue;
                for (c = 0; c < 3; c++)
                  vals[c][n] = p[c];
                n++;
              }
          }
        /* a cell counts when at least half of it is opaque and 16
         * values are left */
        if (n < 16 || 2 * opaque < total)
          continue;
        for (c = 0; c < 3; c++)
          cell->level[c] = clipped_median (vals[c], n, tmp);
        cell->lum = (cell->level[0] + cell->level[1] + cell->level[2]) / 3.0;
        cell->u   = (x0 + x1) / (gdouble) box->width - 1.0;
        cell->v   = (y0 + y1) / (gdouble) box->height - 1.0;
        cell->use = TRUE;
      }

  for (c = 0; c < 3; c++)
    g_free (vals[c]);
  g_free (tmp);
  g_free (pix);
  *n_cells_out = nx * ny;
  return cells;
}

/* drops cells brighter than the median cell by more than tolerance
 * robust sigmas */
static void
reject_bright (Cell    *cells,
               gint     n_cells,
               gdouble  tolerance)
{
  gdouble *v = g_new (gdouble, n_cells);
  gdouble  med, sigma;
  gint     i, n = 0;

  for (i = 0; i < n_cells; i++)
    if (cells[i].use)
      v[n++] = cells[i].lum;
  if (n >= 3)
    {
      med = median_double (v, n);
      for (i = 0, n = 0; i < n_cells; i++)
        if (cells[i].use)
          v[n++] = fabs (cells[i].lum - med);
      sigma = MAD_TO_SIGMA * median_double (v, n);
      if (sigma > 0.0)
        for (i = 0; i < n_cells; i++)
          if (cells[i].use && cells[i].lum > med + tolerance * sigma)
            cells[i].use = FALSE;
    }
  g_free (v);
}

static gint
count_used (const Cell *cells,
            gint        n_cells)
{
  gint i, n = 0;

  for (i = 0; i < n_cells; i++)
    n += cells[i].use;
  return n;
}

/* the surface: fits, drops outlying cells, fits again. Lowers the degree
 * when there are too few cells. Returns the degree used, -1 when there
 * is no sky at all. */
static gint
surface (Cell    *cells,
         gint     n_cells,
         gint     degree,
         gdouble  coef[3][MAX_TERMS])
{
  gdouble *res = g_new (gdouble, n_cells);
  gint     round, i;

  while (degree >= 0 && count_used (cells, n_cells) < n_terms (degree) + 2)
    degree--;
  while (degree >= 0 && ! fit (cells, n_cells, degree, coef))
    degree--;
  if (degree < 0)
    {
      g_free (res);
      return -1;
    }
  for (round = 0; round < FIT_ROUNDS; round++)
    {
      gdouble sigma;
      gint    n = 0, dropped = 0;

      for (i = 0; i < n_cells; i++)
        if (cells[i].use)
          {
            gdouble r = 0.0;
            gint    c;

            for (c = 0; c < 3; c++)
              r += cells[i].level[c] - eval (coef[c], degree, cells[i].u, cells[i].v);
            res[n++] = fabs (r / 3.0);
          }
      sigma = MAD_TO_SIGMA * median_double (res, n);
      if (sigma <= 0.0)
        break;
      for (i = 0; i < n_cells; i++)
        if (cells[i].use)
          {
            gdouble r = 0.0;
            gint    c;

            for (c = 0; c < 3; c++)
              r += cells[i].level[c] - eval (coef[c], degree, cells[i].u, cells[i].v);
            if (fabs (r / 3.0) > RESIDUAL_SIGMA * sigma &&
                n - dropped > n_terms (degree) + 2)
              {
                cells[i].use = FALSE;
                dropped++;
              }
          }
      if (! dropped || ! fit (cells, n_cells, degree, coef))
        break;
    }
  g_free (res);
  return degree;
}

static void
prepare (GeglOperation *operation)
{
  const Babl *space  = gegl_operation_get_source_space (operation, "input");
  const Babl *format = babl_format_with_space ("RGBA float", space);

  gegl_operation_set_format (operation, "input", format);
  gegl_operation_set_format (operation, "output", format);
}

static GeglRectangle
whole_input (GeglOperation       *operation,
             const GeglRectangle *roi)
{
  const GeglRectangle *in = gegl_operation_source_get_bounding_box (operation, "input");

  return in && ! gegl_rectangle_is_infinite_plane (in) ? *in : *roi;
}

static GeglRectangle
get_required_for_output (GeglOperation       *operation,
                         const gchar         *input_pad,
                         const GeglRectangle *roi)
{
  (void) input_pad;
  return whole_input (operation, roi);
}

static GeglRectangle
get_cached_region (GeglOperation       *operation,
                   const GeglRectangle *roi)
{
  return whole_input (operation, roi);
}

/* a change anywhere changes the fit, so the result everywhere */
static GeglRectangle
get_invalidated_by_change (GeglOperation       *operation,
                           const gchar         *input_pad,
                           const GeglRectangle *roi)
{
  (void) input_pad;
  return whole_input (operation, roi);
}

typedef struct
{
  GeglBuffer          *input, *output;
  const Babl          *format;
  const GeglRectangle *box;
  gdouble            (*coef)[MAX_TERMS];
  const gdouble       *mean;
  gint                 degree;
  gboolean             divide, model, keep_level;
} Apply;

/* the correction over one area (threads by GEGL; each pixel's result
 * does not depend on the split) */
static void
apply_area (const GeglRectangle *area,
            gpointer             data)
{
  const Apply        *ap = data;
  const GeglRectangle box = *ap->box;
  gint                degree = ap->degree, c, i;
  GeglBufferIterator *it;

  it = gegl_buffer_iterator_new (ap->output, area, 0, ap->format,
                                 GEGL_ACCESS_WRITE, GEGL_ABYSS_NONE, 2);
  gegl_buffer_iterator_add (it, ap->input, area, 0, ap->format,
                            GEGL_ACCESS_READ, GEGL_ABYSS_NONE);
  while (gegl_buffer_iterator_next (it))
    {
      gfloat              *out = it->items[0].data;
      const gfloat        *in  = it->items[1].data;
      const GeglRectangle *r   = &it->items[0].roi;
      gint                 x, y;

      for (y = 0; y < r->height; y++)
        {
          gdouble v = 2.0 * (r->y + y - box.y + 0.5) / box.height - 1.0;
          /* the surface on this row as a polynomial in u */
          gdouble rowc[3][5];

          for (c = 0; c < 3; c++)
            {
              gint d, k = 0;
              gdouble vp[5];

              vp[0] = 1.0;
              for (d = 1; d <= 4; d++)
                vp[d] = vp[d - 1] * v;
              memset (rowc[c], 0, sizeof rowc[c]);
              for (d = 0; d <= degree; d++)
                for (i = d; i >= 0; i--)
                  rowc[c][i] += ap->coef[c][k++] * vp[d - i];
            }
          for (x = 0; x < r->width; x++, in += 4, out += 4)
            {
              gdouble u = 2.0 * (r->x + x - box.x + 0.5) / box.width - 1.0;

              for (c = 0; c < 3; c++)
                {
                  gdouble b = rowc[c][degree];
                  gint    d;

                  for (d = degree - 1; d >= 0; d--)
                    b = b * u + rowc[c][d];
                  if (ap->model)
                    out[c] = b;
                  else if (ap->divide)
                    out[c] = b > 1e-6 ? in[c] / b * (ap->keep_level ? ap->mean[c] : 1.0) : in[c];
                  else
                    out[c] = in[c] - b + (ap->keep_level ? ap->mean[c] : 0.0);
                }
              out[3] = in[3];
            }
        }
    }
}

static gboolean
process (GeglOperation       *operation,
         GeglBuffer          *input,
         GeglBuffer          *output,
         const GeglRectangle *result,
         gint                 level)
{
  GeglProperties      *o      = GEGL_PROPERTIES (operation);
  const Babl          *format = gegl_operation_get_format (operation, "output");
  const GeglRectangle *boxp   = gegl_operation_source_get_bounding_box (operation, "input");
  GeglRectangle        box;
  gdouble              coef[3][MAX_TERMS];
  gdouble              mean[3] = { 0.0, 0.0, 0.0 };
  gboolean             divide = o->correction == ASTRO_BACKGROUND_DIVIDE;
  gboolean             model  = o->show == ASTRO_BACKGROUND_MODEL;
  Cell                *cells;
  gfloat              *stars = NULL;
  gint                 n_cells, degree, c, i, j;

  (void) level;
  if (! boxp || boxp->width < 1 || boxp->height < 1)
    return TRUE;
  box   = *boxp;
  if (o->exclude_stars)
    {
      /* the star mask at its default settings (astro-common.h) */
      const AstroStarMaskParams params = ASTRO_STAR_MASK_DEFAULT_PARAMS;

      stars = astro_star_mask (input, &box, format, &params, NULL, NULL);
    }
  cells = measure (input, &box, format, o->samples, stars, &n_cells);
  g_free (stars);
  reject_bright (cells, n_cells, o->tolerance);
  degree = surface (cells, n_cells, o->degree, coef);
  g_free (cells);

  if (degree < 0)
    {
      /* nothing to fit (a transparent layer): leave the pixels */
      gegl_buffer_copy (input, result, GEGL_ABYSS_NONE, output, result);
      return TRUE;
    }

  /* the mean level of the surface over the layer, on a 64 x 64 lattice */
  for (j = 0; j < 64; j++)
    for (i = 0; i < 64; i++)
      for (c = 0; c < 3; c++)
        mean[c] += eval (coef[c], degree, (i + 0.5) / 32.0 - 1.0, (j + 0.5) / 32.0 - 1.0);
  for (c = 0; c < 3; c++)
    mean[c] /= 64.0 * 64.0;

  {
    Apply ap = { input, output, format, &box, coef, mean, degree, divide, model,
                 o->keep_level };

    gegl_parallel_distribute_area (result, 4096, GEGL_SPLIT_STRATEGY_HORIZONTAL,
                                   apply_area, &ap);
  }
  return TRUE;
}

static void
gegl_op_class_init (GeglOpClass *klass)
{
  GeglOperationClass       *operation_class = GEGL_OPERATION_CLASS (klass);
  GeglOperationFilterClass *filter_class    = GEGL_OPERATION_FILTER_CLASS (klass);

  operation_class->prepare                   = prepare;
  operation_class->get_required_for_output   = get_required_for_output;
  operation_class->get_cached_region         = get_cached_region;
  operation_class->get_invalidated_by_change = get_invalidated_by_change;
  operation_class->threaded                  = FALSE;
  filter_class->process                      = process;

  gegl_operation_class_set_keys (operation_class,
    "name",            "astro:background",
    "title",           _("Remove Sky Background"),
    "categories",      "enhance",
    "description",     _("Fits a smooth surface to the sky between the "
                         "stars and subtracts it (light pollution, sky "
                         "glow) or divides by it (vignetting)."),
    "gimp:menu-path",  "<Image>/Filters/Astronomy",
    "gimp:menu-label", _("Remove Sky Background..."),
    NULL);
}

#endif
