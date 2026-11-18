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
 * Only the mask's regions are filtered (each with a margin of the
 * radius times the iterations, read from the input, so regions close to
 * each other do not see each other's results). The whole layer is
 * needed for the mask: the operation asks for and caches its whole
 * input and output.
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

/* one pass of the rank filter over a w x h plane: rank 0 to 1 of the
 * values in the disc (edges repeated) */
static void
rank_pass (const gfloat *src,
           gfloat       *dst,
           gint          w,
           gint          h,
           const Disc   *d,
           gdouble       rank,
           gfloat       *tmp)
{
  gint k = (gint) floor (rank * (d->n - 1) + 0.5), x, y, i;

  for (y = 0; y < h; y++)
    for (x = 0; x < w; x++)
      {
        for (i = 0; i < d->n; i++)
          {
            gint xx = CLAMP (x + d->dx[i], 0, w - 1), yy = CLAMP (y + d->dy[i], 0, h - 1);

            tmp[i] = src[yy * w + xx];
          }
        dst[y * w + x] = astro_select (tmp, d->n, k);
      }
}

/* the reduced plane, after all passes */
static void
reduce_plane (gfloat     *plane,
              gfloat     *spare,
              gint        w,
              gint        h,
              const Disc *d,
              gint        method,
              gdouble     selection,
              gint        iterations)
{
  gfloat *tmp = g_new (gfloat, d->n);
  gint    it;

  for (it = 0; it < iterations; it++)
    {
      rank_pass (plane, spare, w, h, d,
                 method == ASTRO_REDUCE_EROSION ? 0.0 : selection, tmp);
      memcpy (plane, spare, (gsize) w * h * sizeof (gfloat));
    }
  g_free (tmp);
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
  gdouble                   amount = o->amount / 100.0;
  gdouble                   lumk[3];
  GeglRectangle             box;
  gfloat                   *mask, *sky;
  guint8                   *seen;
  gint                      w, h, margin, x, y;
  Disc                      disc;
  AstroStack                st = { NULL, 0, 0 };

  (void) level;
  gegl_buffer_copy (input, result, GEGL_ABYSS_NONE, output, result);
  if (! boxp || boxp->width < 1 || boxp->height < 1 || amount <= 0.0)
    return TRUE;
  box  = *boxp;
  w    = box.width;
  h    = box.height;
  mask = astro_star_mask (input, &box, format, &params, NULL, &sky);
  astro_luminance_coefficients (format, lumk);
  disc   = disc_new (o->radius);
  margin = (gint) ceil (o->radius) * o->iterations + 1;
  seen   = g_new0 (guint8, (gsize) w * h);

  /* each connected region of the mask (m > 0), in the order of its
   * first pixel */
  for (y = 0; y < h; y++)
    for (x = 0; x < w; x++)
      {
        gsize         i = (gsize) y * w + x, k;
        gint          x0 = x, x1 = x, y0 = y, y1 = y, rw, rh, c, j;
        GeglRectangle r;
        gfloat       *pix, *plane[3], *spare;
        gint          n_planes;

        if (seen[i] || ! (mask[i] > 0.0f))
          continue;
        seen[i] = 1;
        st.n = 0;
        astro_push (&st, (gint) i);
        for (k = 0; k < st.n; k++)
          {
            gint p = st.stack[k], px = p % w, py = p / w, dx, dy;

            x0 = MIN (x0, px);
            x1 = MAX (x1, px);
            y0 = MIN (y0, py);
            y1 = MAX (y1, py);
            for (dy = -1; dy <= 1; dy++)
              for (dx = -1; dx <= 1; dx++)
                {
                  gint xx = px + dx, yy = py + dy;
                  gsize q;

                  if (xx < 0 || yy < 0 || xx >= w || yy >= h)
                    continue;
                  q = (gsize) yy * w + xx;
                  if (! seen[q] && mask[q] > 0.0f)
                    {
                      seen[q] = 1;
                      astro_push (&st, (gint) q);
                    }
                }
          }

        /* the region with its margin, from the input */
        x0 = MAX (0, x0 - margin);
        y0 = MAX (0, y0 - margin);
        x1 = MIN (w - 1, x1 + margin);
        y1 = MIN (h - 1, y1 + margin);
        rw = x1 - x0 + 1;
        rh = y1 - y0 + 1;
        gegl_rectangle_set (&r, box.x + x0, box.y + y0, rw, rh);
        pix = g_new (gfloat, (gsize) rw * rh * 4);
        gegl_buffer_get (input, &r, 1.0, format, pix, GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);
        n_planes = o->luminance_only ? 1 : 3;
        for (c = 0; c < n_planes; c++)
          plane[c] = g_new (gfloat, (gsize) rw * rh);
        spare = g_new (gfloat, (gsize) rw * rh);
        for (j = 0; j < rw * rh; j++)
          {
            const gfloat *p = pix + (gsize) j * 4;

            if (o->luminance_only)
              plane[0][j] = lumk[0] * p[0] + lumk[1] * p[1] + lumk[2] * p[2];
            else
              for (c = 0; c < 3; c++)
                plane[c][j] = p[c];
          }
        for (c = 0; c < n_planes; c++)
          reduce_plane (plane[c], spare, rw, rh, &disc, o->method, o->selection,
                        o->iterations);

        /* the region's own pixels, blended by the mask; never below the
         * local sky (so the noise around a star is not pulled down into a
         * dark ring) */
        for (k = 0; k < st.n; k++)
          {
            gint    p = st.stack[k], px = p % w - x0, py = p / w - y0;
            gsize   q = (gsize) py * rw + px;
            gfloat *v = pix + q * 4;
            gdouble a = amount * mask[p];
            gdouble l = lumk[0] * v[0] + lumk[1] * v[1] + lumk[2] * v[2];

            if (o->luminance_only)
              {
                gdouble target = MAX (plane[0][q], sky[p]);

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
                  gdouble floor_c = l > 1e-9 ? v[c] * MIN (1.0, sky[p] / l) : v[c];
                  gdouble target  = MAX (plane[c][q], floor_c);

                  if (target < v[c])
                    v[c] += a * (target - v[c]);
                }
          }
        /* written back only where the region is (other pixels of the
         * rectangle keep what the output has) */
        {
          gfloat *cur = g_new (gfloat, (gsize) rw * rh * 4);

          gegl_buffer_get (output, &r, 1.0, format, cur, GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);
          for (k = 0; k < st.n; k++)
            {
              gint  p = st.stack[k];
              gsize q = (gsize) (p / w - y0) * rw + (p % w - x0);

              memcpy (cur + q * 4, pix + q * 4, 3 * sizeof (gfloat));
            }
          gegl_buffer_set (output, &r, 0, format, cur, GEGL_AUTO_ROWSTRIDE);
          g_free (cur);
        }
        for (c = 0; c < n_planes; c++)
          g_free (plane[c]);
        g_free (spare);
        g_free (pix);
      }

  g_free (st.stack);
  g_free (seen);
  g_free (disc.dx);
  g_free (disc.dy);
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
