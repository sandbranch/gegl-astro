/*
 * Star Colour: brings out the colour of the stars under a star mask, as
 * a GEGL operation
 *
 * star-color.c
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
 * PLAN.md op 4, docs/research.md 2.4. The star mask m of astro-common.h
 * (the same settings as Star Mask). Each pixel's chroma in linear RGB,
 * c - Y (Y its luminance in the image's RGB space), is scaled by
 * 1 + saturation x m x b, and Y stays exactly as it was; the scale is
 * limited so that no channel goes below 0. b is Siril's background
 * factor (src/filters/saturation.c, read; written anew here): 0 for
 * pixels darker than (median + sigma) x background, with median and
 * sigma the median and robust sigma of the luminance on a lattice of
 * the whole layer, rising to 1 over one more sigma (Siril switches
 * hard), so the sky's noise between and around the stars does not get
 * coloured.
 *
 * "Fill clipped cores" (our own idea, to be looked at on real data):
 * the core of a bright star clips to white in all three channels, while
 * the ring just outside it keeps the star's colour. Each clipped region
 * under the mask (all channels at 0.95 or above) takes the mean
 * chromaticity (c / Y) of the ring of 2 pixels around it, scaled so that
 * its largest channel keeps its value, blended by the mask.
 *
 * The whole layer is needed (the mask, the statistics): the operation
 * asks for and caches its whole input and output.
 */

#include <glib/gi18n-lib.h>
#include "astro-common.h"

#ifdef GEGL_PROPERTIES

property_double (saturation, _("Saturation"), 50.0)
  description (_("How much the stars' colour is raised, in percent of "
                 "their chroma."))
  value_range (0.0, 300.0)
  ui_range (0.0, 200.0)
  ui_digits (0)

property_double (background, _("Background factor"), 1.0)
  description (_("Pixels darker than (median + sigma) times this stay as "
                 "they are, so the sky's noise is not coloured (as "
                 "Siril's saturation). 0 colours everything under the "
                 "mask."))
  value_range (0.0, 10.0)
  ui_range (0.0, 3.0)
  ui_digits (2)

property_boolean (fill_cores, _("Fill clipped cores"), TRUE)
  description (_("The white, clipped cores of bright stars take the "
                 "colour of the ring around them."))

ASTRO_STAR_MASK_PROPERTIES

#else

#define GEGL_OP_FILTER
#define GEGL_OP_NAME     star_color
#define GEGL_OP_C_SOURCE star-color.c

#include "gegl-op.h"

#define CLIP_LEVEL 0.95f
#define RING       2

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

/* the median and robust sigma of the luminance of opaque pixels, on a
 * fixed lattice of the whole box */
static void
luminance_stats (GeglBuffer          *input,
                 const GeglRectangle *box,
                 const Babl          *format,
                 const gdouble        k[3],
                 gdouble             *median,
                 gdouble             *sigma)
{
  gint    step = astro_lattice_step (box->width, box->height, ASTRO_STATS_SIDE);
  gint    nx = (box->width + step - 1) / step, ny = (box->height + step - 1) / step;
  gfloat *row = g_new (gfloat, (gsize) box->width * 4);
  gfloat *v = g_new (gfloat, (gsize) nx * ny);
  gint    x, y, n = 0;

  for (y = step / 2; y < box->height; y += step)
    {
      gegl_buffer_get (input, GEGL_RECTANGLE (box->x, box->y + y, box->width, 1), 1.0,
                       format, row, GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);
      for (x = step / 2; x < box->width; x += step)
        {
          const gfloat *p = row + (gsize) x * 4;

          if (p[3] > 0.5f && isfinite (p[0]) && isfinite (p[1]) && isfinite (p[2]))
            v[n++] = k[0] * p[0] + k[1] * p[1] + k[2] * p[2];
        }
    }
  astro_robust_stats (v, n, median, sigma);
  g_free (v);
  g_free (row);
}

/* the chroma of p scaled by f about its luminance y, limited so that no
 * channel goes below 0 */
static inline void
boost (gfloat  *p,
       gdouble  y,
       gdouble  f)
{
  gint c;

  if (! (y > 0.0))
    return;
  for (c = 0; c < 3; c++)
    if (p[c] < y && f > 1.0)
      /* p[c] + (f - 1) (p[c] - y) >= 0 */
      f = MIN (f, 1.0 + p[c] / (y - p[c]));
  for (c = 0; c < 3; c++)
    p[c] = y + (p[c] - y) * f;
}

#define STRIP_ROWS 128

typedef struct
{
  gsize   index;     /* the pixel, in the box */
  gfloat  rgb[3];    /* its colour (the ring's, scaled), and */
  gfloat  m;         /* its mask */
} Core;

typedef struct
{
  GeglBuffer    *input, *output;
  const Babl    *format;
  GeglRectangle  box, result;
  const gfloat  *mask;
  guint8        *clipped;   /* NULL, or set by the first pass */
  const Core    *cores;
  gsize          n_cores;
  gdouble        k[3], amount, lo, hi;
  gboolean       write;     /* the second pass */
} Pass;

static int
cmp_core (const void *a,
          const void *b)
{
  gsize x = ((const Core *) a)->index, y = ((const Core *) b)->index;

  return (x > y) - (x < y);
}

/* first pass: marks the clipped pixels under the mask; second pass:
 * the boost and the cores' colours, written */
static void
pass_strips (gsize    offset,
             gsize    size,
             gpointer data)
{
  const Pass *ps = data;
  gint        w = ps->box.width, h = ps->box.height;
  gfloat     *pix = g_new (gfloat, (gsize) w * STRIP_ROWS * 4);
  gsize       s;

  for (s = offset; s < offset + size; s++)
    {
      gint  y0 = (gint) s * STRIP_ROWS, y1 = MIN (h, y0 + STRIP_ROWS);
      gsize i0 = (gsize) y0 * w, i1 = (gsize) y1 * w, i, lo, hi;

      gegl_buffer_get (ps->input, GEGL_RECTANGLE (ps->box.x, ps->box.y + y0, w, y1 - y0),
                       1.0, ps->format, pix, GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);
      if (! ps->write)
        {
          for (i = i0; i < i1; i++)
            {
              const gfloat *p = pix + (i - i0) * 4;

              ps->clipped[i] = ps->mask[i] > 0.0f && p[0] >= CLIP_LEVEL &&
                               p[1] >= CLIP_LEVEL && p[2] >= CLIP_LEVEL;
            }
          continue;
        }
      if (ps->amount > 0.0)
        for (i = i0; i < i1; i++)
          {
            gfloat *p = pix + (i - i0) * 4;
            gdouble y, t;

            if (! (ps->mask[i] > 0.0f))
              continue;
            y = ps->k[0] * p[0] + ps->k[1] * p[1] + ps->k[2] * p[2];
            t = ps->hi > ps->lo ? CLAMP ((y - ps->lo) / (ps->hi - ps->lo), 0.0, 1.0)
                                : (y >= ps->lo);
            if (t > 0.0)
              boost (p, y, 1.0 + ps->amount * ps->mask[i] * t);
          }
      /* the cores in this strip (sorted by index: the first at or after
       * i0, by bisection) */
      lo = 0;
      hi = ps->n_cores;
      while (lo < hi)
        {
          gsize mid = (lo + hi) / 2;

          if (ps->cores[mid].index < i0)
            lo = mid + 1;
          else
            hi = mid;
        }
      for (i = lo; i < ps->n_cores && ps->cores[i].index < i1; i++)
        {
          const Core *c = &ps->cores[i];
          gfloat     *p = pix + (c->index - i0) * 4;
          gint        ch;

          for (ch = 0; ch < 3; ch++)
            p[ch] = (1.0f - c->m) * p[ch] + c->m * c->rgb[ch];
        }
      {
        GeglRectangle rows, part;

        gegl_rectangle_set (&rows, ps->box.x, ps->box.y + y0, w, y1 - y0);
        if (gegl_rectangle_intersect (&part, &rows, &ps->result))
          gegl_buffer_set (ps->output, &part, 0, ps->format,
                           pix + ((gsize) (part.y - rows.y) * w + (part.x - rows.x)) * 4,
                           (gint) (w * 4 * sizeof (gfloat)));
      }
    }
  g_free (pix);
}

/* the clipped regions: each takes the mean chromaticity (c / Y) of the
 * unclipped pixels under the mask within RING of it, scaled so that its
 * largest channel keeps the pixel's largest value; a list sorted by
 * pixel */
static Core *
find_cores (const Pass *ps,
            gsize      *n_out)
{
  gint       w = ps->box.width, h = ps->box.height;
  gsize      n = (gsize) w * h, i, n_cores = 0, size = 0;
  guint8    *seen = g_new0 (guint8, n);
  AstroStack st = { NULL, 0, 0 };
  Core      *cores = NULL;

  for (i = 0; i < n; i++)
    if (ps->clipped[i] && ! seen[i])
      {
        gsize         q;
        gint          x0 = w, x1 = -1, y0 = h, y1 = -1, x, y, rw, rh, ring = 0, c;
        gdouble       ratio[3] = { 0, 0, 0 }, top;
        gfloat       *pix;
        GeglRectangle r;

        st.n = 0;
        seen[i] = 1;
        astro_push (&st, (gint) i);
        for (q = 0; q < st.n; q++)
          {
            gint p = st.stack[q], px = p % w, py = p / w, dx, dy;

            x0 = MIN (x0, px);
            x1 = MAX (x1, px);
            y0 = MIN (y0, py);
            y1 = MAX (y1, py);
            for (dy = -1; dy <= 1; dy++)
              for (dx = -1; dx <= 1; dx++)
                {
                  gint  xx = px + dx, yy = py + dy;
                  gsize j;

                  if (xx < 0 || yy < 0 || xx >= w || yy >= h)
                    continue;
                  j = (gsize) yy * w + xx;
                  if (ps->clipped[j] && ! seen[j])
                    {
                      seen[j] = 1;
                      astro_push (&st, (gint) j);
                    }
                }
          }
        x0 = MAX (0, x0 - RING);
        y0 = MAX (0, y0 - RING);
        x1 = MIN (w - 1, x1 + RING);
        y1 = MIN (h - 1, y1 + RING);
        rw = x1 - x0 + 1;
        rh = y1 - y0 + 1;
        gegl_rectangle_set (&r, ps->box.x + x0, ps->box.y + y0, rw, rh);
        pix = g_new (gfloat, (gsize) rw * rh * 4);
        gegl_buffer_get (ps->input, &r, 1.0, ps->format, pix, GEGL_AUTO_ROWSTRIDE,
                         GEGL_ABYSS_NONE);
        for (y = y0; y <= y1; y++)
          for (x = x0; x <= x1; x++)
            {
              gsize         j = (gsize) y * w + x;
              const gfloat *p = pix + ((gsize) (y - y0) * rw + (x - x0)) * 4;
              gboolean      near = FALSE;
              gint          dx, dy;
              gdouble       yl;

              if (ps->clipped[j] || ! (ps->mask[j] > 0.0f))
                continue;
              /* within RING of this region (its pixels are seen == 1) */
              for (dy = -RING; dy <= RING && ! near; dy++)
                for (dx = -RING; dx <= RING; dx++)
                  {
                    gint xx = x + dx, yy = y + dy;

                    if (xx >= 0 && yy >= 0 && xx < w && yy < h &&
                        seen[(gsize) yy * w + xx] == 1)
                      {
                        near = TRUE;
                        break;
                      }
                  }
              yl = ps->k[0] * p[0] + ps->k[1] * p[1] + ps->k[2] * p[2];
              if (! near || ! (yl > 0.0))
                continue;
              for (c = 0; c < 3; c++)
                ratio[c] += p[c] / yl;
              ring++;
            }
        top = MAX (ratio[0], MAX (ratio[1], ratio[2]));
        if (ring > 0 && top > 0.0)
          for (q = 0; q < st.n; q++)
            {
              gsize         j = st.stack[q];
              gint          px = (gint) (j % w) - x0, py = (gint) (j / w) - y0;
              const gfloat *p = pix + ((gsize) py * rw + px) * 4;
              gdouble       v = MAX (p[0], MAX (p[1], p[2]));
              Core         *core;

              if (n_cores == size)
                {
                  size  = MAX (256, size * 2);
                  cores = g_renew (Core, cores, size);
                }
              core = &cores[n_cores++];
              core->index = j;
              core->m     = ps->mask[j];
              for (c = 0; c < 3; c++)
                core->rgb[c] = (gfloat) (v * ratio[c] / top);
            }
        g_free (pix);
        /* this region is done: 2, so the next region's ring test does
         * not see it */
        for (q = 0; q < st.n; q++)
          seen[st.stack[q]] = 2;
      }
  g_free (st.stack);
  g_free (seen);
  if (n_cores > 1)
    qsort (cores, n_cores, sizeof (Core), cmp_core);
  *n_out = n_cores;
  return cores;
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
  gdouble                   median, sigma;
  gint                      n_strips;
  Core                     *cores = NULL;
  Pass                      ps;

  (void) level;
  if (! boxp || boxp->width < 1 || boxp->height < 1 ||
      (o->saturation <= 0.0 && ! o->fill_cores))
    {
      gegl_buffer_copy (input, result, GEGL_ABYSS_NONE, output, result);
      return TRUE;
    }
  memset (&ps, 0, sizeof ps);
  ps.input  = input;
  ps.output = output;
  ps.format = format;
  ps.box    = *boxp;
  ps.result = *result;
  ps.mask   = astro_star_mask (input, boxp, format, &params, NULL, NULL);
  astro_luminance_coefficients (format, ps.k);
  luminance_stats (input, boxp, format, ps.k, &median, &sigma);
  ps.amount = o->saturation / 100.0;
  ps.lo     = (median + sigma) * o->background;
  ps.hi     = ps.lo + sigma;
  if (o->background <= 0.0)
    ps.lo = ps.hi = -G_MAXDOUBLE;
  n_strips = (boxp->height + STRIP_ROWS - 1) / STRIP_ROWS;

  if (o->fill_cores)
    {
      ps.clipped = g_new0 (guint8, (gsize) boxp->width * boxp->height);
      gegl_parallel_distribute_range (n_strips, 1, pass_strips, &ps);
      cores = find_cores (&ps, &ps.n_cores);
      ps.cores = cores;
      g_free (ps.clipped);
      ps.clipped = NULL;
    }
  ps.write = TRUE;
  gegl_parallel_distribute_range (n_strips, 1, pass_strips, &ps);
  g_free (cores);
  g_free ((gpointer) ps.mask);
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
    "name",            "astro:star-color",
    "title",           _("Star Colour"),
    "categories",      "color",
    "description",     _("Brings out the colour of the stars under a star "
                         "mask and gives clipped white cores the colour of "
                         "their halo; the sky and nebulae stay as they "
                         "are."),
    "gimp:menu-path",  "<Image>/Filters/Astronomy",
    "gimp:menu-label", _("Star Colour..."),
    NULL);
}

#endif
