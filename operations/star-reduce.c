/*
 * Reduce Stars: makes the stars of an astronomical image smaller and
 * fainter under a star mask, as a GEGL operation
 *
 * star-reduce.c
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
 * PLAN.md op 3, docs/research.md 2.2. The star mask m of astro-common.h
 * (the same settings as Star Mask), then a rank filter over a disc of
 * the given radius, repeated: the minimum (erosion), the value at a
 * rank between the minimum and the maximum ("selection": 0 is erosion,
 * 0.5 the median, 1 dilation; the published description of PixInsight's
 * Morphological Selection, none of its code). Classical grey-scale
 * morphology (Matheron, Serra). (An opening, erosion then dilation, was
 * tried and left out: it flattens a star's peak, so the star gets
 * fainter but wider at half maximum, 1.30 x on the test field.) Then out = in + amount x m x (rank - in), where rank is not
 * taken below the local sky (the luminance blurred at the mask's
 * max-size): the noise of the sky around a star is not pulled down into
 * a dark ring, and only darker values are taken. Outside the mask
 * nothing changes at all. "Keep star colour": the rank filter works on
 * the luminance, and each pixel's channels are scaled by the new over
 * the old luminance, so stars keep their colour; in this mode a pixel
 * only gets darker (a scale above 1 at faint pixels would blow up their
 * colour noise).
 *
 * The layer is done in strips of rows, in parallel; each strip reads
 * its rows and a margin of the disc's reach times the iterations from
 * the input and filters only where the mask needs it (the mask grown by
 * the reach of the later passes), so the result does not depend on the
 * strips. The whole layer is needed for the mask: the operation asks
 * for and caches its whole input and output.
 */

#include <glib/gi18n-lib.h>
#include "astro-common.h"

#ifdef GEGL_PROPERTIES

enum_start (astro_star_reduce_method)
  enum_value (ASTRO_REDUCE_EROSION,   "erosion",   N_("Erosion (minimum)"))
  enum_value (ASTRO_REDUCE_SELECTION, "selection", N_("Selection (rank)"))
enum_end (AstroStarReduceMethod)

property_double (amount, _("Amount"), 70.0)
  description (_("How far the stars move towards their reduced form, "
                 "in percent."))
  value_range (0.0, 100.0)
  ui_digits (0)

property_double (radius, _("Radius"), 1.5)
  description (_("Radius of the round neighbourhood, in pixels; a star "
                 "shrinks by about this much per iteration."))
  value_range (0.5, 8.0)
  ui_digits (1)

property_enum (method, _("Method"),
               AstroStarReduceMethod, astro_star_reduce_method,
               ASTRO_REDUCE_SELECTION)
  description (_("Erosion takes the darkest neighbour; selection a "
                 "neighbour at the given rank (gentler)."))

property_double (selection, _("Selection"), 0.25)
  description (_("The rank for the selection method: 0 the darkest "
                 "neighbour (erosion), 0.5 the median. Above 0.5 little "
                 "changes, as only darker values are taken."))
  value_range (0.0, 1.0)
  ui_digits (2)
  ui_meta ("visible", "method {selection}")

property_int (iterations, _("Iterations"), 1)
  description (_("How many times the filter is applied."))
  value_range (1, 4)

property_boolean (luminance_only, _("Keep star colour"), TRUE)
  description (_("Works on the brightness only and keeps each pixel's "
                 "colour; off: each channel on its own."))

ASTRO_STAR_MASK_PROPERTIES

#else

#define GEGL_OP_FILTER
#define GEGL_OP_NAME     star_reduce
#define GEGL_OP_C_SOURCE star-reduce.c

#include "gegl-op.h"

static void
prepare (GeglOperation *operation)
{
  const Babl *space  = gegl_operation_get_source_space (operation, "input");
  const Babl *format = babl_format_with_space ("RGBA float", space);

  gegl_operation_set_format (operation, "input", format);
  gegl_operation_set_format (operation, "output", format);
}

static GeglRectangle
get_required_for_output (GeglOperation       *operation,
                         const gchar         *input_pad,
                         const GeglRectangle *roi)
{
  (void) input_pad;
  return astro_whole_input (operation, roi);
}

static GeglRectangle
get_cached_region (GeglOperation       *operation,
                   const GeglRectangle *roi)
{
  return astro_whole_input (operation, roi);
}

static GeglRectangle
get_invalidated_by_change (GeglOperation       *operation,
                           const gchar         *input_pad,
                           const GeglRectangle *roi)
{
  (void) input_pad;
  return astro_whole_input (operation, roi);
}

typedef struct
{
  gint  n;
  gint *dx, *dy;
} Disc;

static Disc
disc_new (gdouble radius)
{
  gint    r = (gint) floor (radius), x, y;
  gdouble r2 = radius * radius + 1e-9;
  Disc    d;

  d.n  = 0;
  d.dx = g_new (gint, (2 * r + 1) * (2 * r + 1));
  d.dy = g_new (gint, (2 * r + 1) * (2 * r + 1));
  for (y = -r; y <= r; y++)
    for (x = -r; x <= r; x++)
      if (x * x + y * y <= r2)
        {
          d.dx[d.n] = x;
          d.dy[d.n] = y;
          d.n++;
        }
  return d;
}

#define STRIP_ROWS 128

typedef struct
{
  GeglBuffer          *input, *output;
  const Babl          *format;
  GeglRectangle        box, result;
  const gfloat        *mask, *sky;
  Disc                 disc;
  gdouble              lumk[3], amount, rank;
  gint                 iterations, step;
  gboolean             luminance_only;
} Reduce;

/* the strip's rows within the result, written from rows of pix (which
 * start at box row e0) */
static void
write_rows (const Reduce *rd,
            const gfloat *pix,
            gint          e0,
            gint          y0,
            gint          y1)
{
  GeglRectangle rows, part;
  gint          w = rd->box.width;

  gegl_rectangle_set (&rows, rd->box.x, rd->box.y + y0, w, y1 - y0);
  if (! gegl_rectangle_intersect (&part, &rows, &rd->result))
    return;
  gegl_buffer_set (rd->output, &part, 0, rd->format,
                   pix + ((gsize) (part.y - rd->box.y - e0) * w + (part.x - rd->box.x)) * 4,
                   (gint) (w * 4 * sizeof (gfloat)));
}

static void
reduce_strip (const Reduce *rd,
              gint          y0,
              gint          y1)
{
  gint     w = rd->box.width, h = rd->box.height;
  gint     reach = rd->step * rd->iterations;
  gint     e0 = MAX (0, y0 - reach), e1 = MIN (h, y1 + reach), eh = e1 - e0;
  gint     n_planes = rd->luminance_only ? 1 : 3, c, it, x, y, i, k;
  gint     grow = rd->step * (rd->iterations - 1);
  gsize    n = (gsize) w * eh, j;
  gboolean any = FALSE;
  gfloat  *pix, *cur[3], *nxt[3], *tmp;
  guint8  *act, *row;

  for (y = y0; y < y1 && ! any; y++)
    for (x = 0; x < w; x++)
      if (rd->mask[(gsize) y * w + x] > 0.0f)
        {
          any = TRUE;
          break;
        }
  pix = g_new (gfloat, n * 4);
  gegl_buffer_get (rd->input, GEGL_RECTANGLE (rd->box.x, rd->box.y + e0, w, eh), 1.0,
                   rd->format, pix, GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);
  if (! any)
    {
      write_rows (rd, pix, e0, y0, y1);
      g_free (pix);
      return;
    }

  /* where the filter is needed: the strip's mask, grown by the reach of
   * the later passes (a square, separably) */
  act = g_new0 (guint8, n);
  row = g_new (guint8, MAX (w, eh));
  for (y = y0; y < y1; y++)
    for (x = 0; x < w; x++)
      act[(gsize) (y - e0) * w + x] = rd->mask[(gsize) y * w + x] > 0.0f;
  if (grow > 0)
    {
      for (y = 0; y < eh; y++)
        {
          guint8 *a = act + (gsize) y * w;
          gint    last = -1 - grow;

          memcpy (row, a, w);
          /* distance to the nearest set pixel on the left, then right */
          for (x = 0; x < w; x++)
            {
              if (row[x])
                last = x;
              a[x] = x - last <= grow;
            }
          last = w + grow;
          for (x = w - 1; x >= 0; x--)
            {
              if (row[x])
                last = x;
              a[x] |= last - x <= grow;
            }
        }
      for (x = 0; x < w; x++)
        {
          gint last = -1 - grow;

          for (y = 0; y < eh; y++)
            row[y] = act[(gsize) y * w + x];
          for (y = 0; y < eh; y++)
            {
              if (row[y])
                last = y;
              act[(gsize) y * w + x] = y - last <= grow;
            }
          last = eh + grow;
          for (y = eh - 1; y >= 0; y--)
            {
              if (row[y])
                last = y;
              act[(gsize) y * w + x] |= last - y <= grow;
            }
        }
    }
  g_free (row);

  for (c = 0; c < n_planes; c++)
    {
      cur[c] = g_new (gfloat, n);
      nxt[c] = g_new (gfloat, n);
    }
  for (j = 0; j < n; j++)
    {
      const gfloat *p = pix + j * 4;

      if (rd->luminance_only)
        cur[0][j] = rd->lumk[0] * p[0] + rd->lumk[1] * p[1] + rd->lumk[2] * p[2];
      else
        for (c = 0; c < 3; c++)
          cur[c][j] = p[c];
    }

  /* the passes, at the needed pixels only (edges repeated) */
  tmp = g_new (gfloat, rd->disc.n);
  k   = (gint) floor (rd->rank * (rd->disc.n - 1) + 0.5);
  for (it = 0; it < rd->iterations; it++)
    for (c = 0; c < n_planes; c++)
      {
        gfloat *s = cur[c], *d = nxt[c];

        memcpy (d, s, n * sizeof (gfloat));
        for (y = 0; y < eh; y++)
          for (x = 0; x < w; x++)
            {
              if (! act[(gsize) y * w + x])
                continue;
              for (i = 0; i < rd->disc.n; i++)
                {
                  gint xx = CLAMP (x + rd->disc.dx[i], 0, w - 1);
                  gint yy = CLAMP (y + rd->disc.dy[i], 0, eh - 1);

                  tmp[i] = s[(gsize) yy * w + xx];
                }
              d[(gsize) y * w + x] = astro_select (tmp, rd->disc.n, k);
            }
        cur[c] = d;
        nxt[c] = s;
      }
  g_free (tmp);

  /* the strip's own rows, blended by the mask; never below the local sky
   * (so the noise around a star is not pulled down into a dark ring) */
  for (y = y0; y < y1; y++)
    for (x = 0; x < w; x++)
      {
        gsize   p = (gsize) y * w + x, q = (gsize) (y - e0) * w + x;
        gfloat *v = pix + q * 4;
        gdouble a, l;

        if (! (rd->mask[p] > 0.0f))
          continue;
        a = rd->amount * rd->mask[p];
        l = rd->lumk[0] * v[0] + rd->lumk[1] * v[1] + rd->lumk[2] * v[2];
        if (rd->luminance_only)
          {
            gdouble target = MAX (cur[0][q], rd->sky[p]);

            /* (dark or negative: the scale is not defined, left) */
            if (l > 1e-9 && target < l)
              {
                gdouble f = 1.0 + a * (target / l - 1.0);

                for (c = 0; c < 3; c++)
                  v[c] *= f;
              }
          }
        else
          for (c = 0; c < 3; c++)
            {
              /* the sky's share of this channel */
              gdouble floor_c = l > 1e-9 ? v[c] * MIN (1.0, rd->sky[p] / l) : v[c];
              gdouble target  = MAX (cur[c][q], floor_c);

              if (target < v[c])
                v[c] += a * (target - v[c]);
            }
      }
  write_rows (rd, pix, e0, y0, y1);

  for (c = 0; c < n_planes; c++)
    {
      g_free (cur[c]);
      g_free (nxt[c]);
    }
  g_free (act);
  g_free (pix);
}

static void
reduce_strips (gsize    offset,
               gsize    size,
               gpointer data)
{
  const Reduce *rd = data;
  gsize         s;

  for (s = offset; s < offset + size; s++)
    reduce_strip (rd, (gint) s * STRIP_ROWS,
                  MIN (rd->box.height, (gint) (s + 1) * STRIP_ROWS));
}

static gboolean
process (GeglOperation       *operation,
         GeglBuffer          *input,
         GeglBuffer          *output,
         const GeglRectangle *result,
         gint                 level)
{
  GeglProperties           *o      = GEGL_PROPERTIES (operation);
  const Babl               *format = gegl_operation_get_format (operation, "output");
  const GeglRectangle      *boxp   = gegl_operation_source_get_bounding_box (operation, "input");
  const AstroStarMaskParams params = ASTRO_STAR_MASK_PARAMS (o);
  gfloat                   *mask, *sky;
  Reduce                    rd;

  (void) level;
  if (! boxp || boxp->width < 1 || boxp->height < 1 || o->amount <= 0.0)
    {
      gegl_buffer_copy (input, result, GEGL_ABYSS_NONE, output, result);
      return TRUE;
    }
  mask = astro_star_mask (input, boxp, format, &params, NULL, &sky);

  rd.input          = input;
  rd.output         = output;
  rd.format         = format;
  rd.box            = *boxp;
  rd.result         = *result;
  rd.mask           = mask;
  rd.sky            = sky;
  rd.disc           = disc_new (o->radius);
  rd.amount         = o->amount / 100.0;
  rd.rank           = o->method == ASTRO_REDUCE_EROSION ? 0.0 : o->selection;
  rd.iterations     = o->iterations;
  rd.step           = (gint) floor (o->radius);
  rd.luminance_only = o->luminance_only;
  astro_luminance_coefficients (format, rd.lumk);
  gegl_parallel_distribute_range ((boxp->height + STRIP_ROWS - 1) / STRIP_ROWS, 1,
                                  reduce_strips, &rd);
  g_free (rd.disc.dx);
  g_free (rd.disc.dy);
  g_free (mask);
  g_free (sky);
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
    "name",            "astro:star-reduce",
    "title",           _("Reduce Stars"),
    "categories",      "enhance",
    "description",     _("Makes the stars smaller and fainter under a star "
                         "mask, so nebulae and galaxies stand out; the sky "
                         "and the nebulae outside the mask stay as they "
                         "are."),
    "gimp:menu-path",  "<Image>/Filters/Astronomy",
    "gimp:menu-label", _("Reduce Stars..."),
    NULL);
}

#endif
