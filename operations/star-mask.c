/*
 * Star Mask: a mask of the stars of an astronomical image, as a GEGL
 * operation
 *
 * star-mask.c
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
 * The mask (white = star) to look at, to copy into a layer mask, or to
 * guide other tools; the same mask as inside Reduce Stars and Star
 * Colour with the same settings (astro-common.h has the method). Shown
 * as grey, or over the image in red. Alpha is kept. The whole layer is
 * needed: the operation asks for and caches its whole input and output.
 */

#include <glib/gi18n-lib.h>
#include "astro-common.h"

#ifdef GEGL_PROPERTIES

enum_start (astro_star_mask_show)
  enum_value (ASTRO_STAR_MASK_MASK,    "mask",    N_("Mask"))
  enum_value (ASTRO_STAR_MASK_OVERLAY, "overlay", N_("Over the image, in red"))
enum_end (AstroStarMaskShow)

ASTRO_STAR_MASK_PROPERTIES

property_enum (show, _("Show"),
               AstroStarMaskShow, astro_star_mask_show, ASTRO_STAR_MASK_MASK)
  description (_("The mask in grey (white = star), or the image with the "
                 "mask in red."))

#else

#define GEGL_OP_FILTER
#define GEGL_OP_NAME     star_mask
#define GEGL_OP_C_SOURCE star-mask.c

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
  gboolean                  overlay = o->show == ASTRO_STAR_MASK_OVERLAY;
  GeglRectangle             box;
  GeglBufferIterator       *it;
  gfloat                   *mask;

  (void) level;
  if (! boxp || boxp->width < 1 || boxp->height < 1)
    return TRUE;
  box  = *boxp;
  mask = astro_star_mask (input, &box, format, &params, NULL, NULL);

  it = gegl_buffer_iterator_new (output, result, 0, format,
                                 GEGL_ACCESS_WRITE, GEGL_ABYSS_NONE, 2);
  gegl_buffer_iterator_add (it, input, result, 0, format,
                            GEGL_ACCESS_READ, GEGL_ABYSS_NONE);
  while (gegl_buffer_iterator_next (it))
    {
      gfloat              *out = it->items[0].data;
      const gfloat        *in  = it->items[1].data;
      const GeglRectangle *r   = &it->items[0].roi;
      gint                 x, y;

      for (y = 0; y < r->height; y++)
        {
          const gfloat *m = mask + (gsize) (r->y + y - box.y) * box.width + (r->x - box.x);

          for (x = 0; x < r->width; x++, in += 4, out += 4)
            {
              if (overlay)
                {
                  gfloat a = 0.6f * m[x];

                  out[0] = in[0] * (1.0f - a) + a;
                  out[1] = in[1] * (1.0f - a);
                  out[2] = in[2] * (1.0f - a);
                }
              else
                out[0] = out[1] = out[2] = m[x];
              out[3] = in[3];
            }
        }
    }
  g_free (mask);
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
    "name",            "astro:star-mask",
    "title",           _("Star Mask"),
    "categories",      "enhance",
    "description",     _("A mask of the stars (white = star): small, "
                         "concentrated points of light above the noise; "
                         "nebulae and the sky are left out."),
    "gimp:menu-path",  "<Image>/Filters/Astronomy",
    "gimp:menu-label", _("Star Mask..."),
    NULL);
}

#endif
