/*
 * Stretch: the arcsinh stretch and the Generalised Hyperbolic Stretch
 * (GHS) of linear astronomical data, as a GEGL operation
 *
 * stretch.c
 * Copyright 2026 by David
 *
 * The GHS coefficients and the transform (ghs_setup and ghs_eval below)
 * are ported from Siril's src/filters/ght.c (GHTsetup, GHT), master
 * 499fd209 of 2026-11-09:
 *
 *   This file is part of Siril, an astronomy image processor.
 *   Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 *   Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 *   Reference site is https://siril.org
 *
 * Siril's ght.c is under the GNU General Public License version 3 or
 * (at your option) any later version, as is this file. The GHS equations
 * are by Dave Payne and Mike Cranfield (ghsastro.co.uk). Changes in the
 * port: double instead of float, only the forward ("normal") and the
 * inverse transform of the five local intensity cases (Siril's GHS asinh
 * and linear modes are left out), the protection points kept in order
 * (LP <= SP <= HP) instead of refused, and one fix: in the inverse for
 * b = -1, the transformed symmetry point is - q0 q (Siril: + q0 q, so
 * its inverse took the wrong branch below SP; found by the check
 * against the published equations, tests/stretch-reference.py).
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
 * Both curves map 0 to 0 and 1 to 1 and rise monotonically. The black
 * point is taken off first: x' = max (0, (x - bp) / (1 - bp)), and the
 * values are clipped to 0 to 1 (as Siril does). "Preserve colour": the
 * curve is applied to the luminance Y of the image's RGB space, and all
 * three channels are multiplied by curve (Y) / Y, so each pixel keeps its
 * hue and its ratios of R, G and B (Lupton et al. 2004); a channel that
 * would go above 1 scales the pixel down instead (Siril's "rescale").
 * "Per channel": the curve on each channel alone.
 *
 * A point operation, except with "automatic black and symmetry point":
 * then the median and the robust sigma (1.4826 x MAD) of the luminance
 * (or of each channel) on a fixed lattice of at most 256 x 256 points of
 * the whole layer set bp = max (0, median - 3 sigma) and the symmetry
 * point to where the median lands; the operation then asks for its whole
 * input.
 */

#include <glib/gi18n-lib.h>

#ifdef GEGL_PROPERTIES

enum_start (astro_stretch_curve)
  enum_value (ASTRO_STRETCH_GHS,     "ghs",     N_("Generalised hyperbolic (GHS)"))
  enum_value (ASTRO_STRETCH_ARCSINH, "arcsinh", N_("Arcsinh"))
enum_end (AstroStretchCurve)

enum_start (astro_stretch_colour)
  enum_value (ASTRO_STRETCH_PRESERVE,    "preserve",    N_("Preserve colour"))
  enum_value (ASTRO_STRETCH_PER_CHANNEL, "per-channel", N_("Per channel"))
enum_end (AstroStretchColour)

property_enum (curve, _("Curve"),
               AstroStretchCurve, astro_stretch_curve, ASTRO_STRETCH_GHS)
  description (_("The generalised hyperbolic stretch (five parameters) "
                 "or the arcsinh stretch (one)."))

property_double (black_point, _("Black point"), 0.0)
  description (_("Taken off first; the rest is stretched to fill 0 to 1."))
  value_range (0.0, 0.99)
  ui_range (0.0, 0.2)
  ui_digits (4)

property_double (stretch, _("Stretch (ln (D + 1))"), 0.0)
  description (_("GHS: the strength of the stretch; 0 leaves the image "
                 "as it is."))
  value_range (0.0, 10.0)
  ui_digits (3)
  ui_meta ("visible", "curve {ghs}")

property_double (local_intensity, _("Local intensity (b)"), 0.0)
  description (_("GHS: how much the stretch is focused around the "
                 "symmetry point. -1 logarithmic, 0 exponential, 1 "
                 "harmonic; higher values narrow the range of brightness "
                 "that is stretched."))
  value_range (-5.0, 15.0)
  ui_digits (2)
  ui_meta ("visible", "curve {ghs}")

property_double (symmetry_point, _("Symmetry point (SP)"), 0.0)
  description (_("GHS: the brightness where the stretch adds the most "
                 "contrast, usually just above the sky background."))
  value_range (0.0, 1.0)
  ui_digits (4)
  ui_meta ("visible", "curve {ghs}")

property_double (protect_shadows, _("Shadow protection (LP)"), 0.0)
  description (_("GHS: below this level the curve is a straight line, "
                 "so dark noise is not stretched. At most the symmetry "
                 "point."))
  value_range (0.0, 1.0)
  ui_digits (4)
  ui_meta ("visible", "curve {ghs}")

property_double (protect_highlights, _("Highlight protection (HP)"), 1.0)
  description (_("GHS: above this level the curve is a straight line, "
                 "so bright stars keep their shape. At least the "
                 "symmetry point."))
  value_range (0.0, 1.0)
  ui_digits (4)
  ui_meta ("visible", "curve {ghs}")

property_boolean (inverse, _("Inverse"), FALSE)
  description (_("GHS: the inverse transform, which undoes the same "
                 "stretch."))
  ui_meta ("visible", "curve {ghs}")

property_double (beta, _("Arcsinh strength"), 10.0)
  description (_("Arcsinh: x becomes asinh (beta x) / asinh (beta); 0 "
                 "leaves the image as it is."))
  value_range (0.0, 1000.0)
  ui_range (0.0, 200.0)
  ui_gamma (2.0)
  ui_digits (1)
  ui_meta ("visible", "curve {arcsinh}")

property_enum (colour, _("Colour"),
               AstroStretchColour, astro_stretch_colour, ASTRO_STRETCH_PRESERVE)
  description (_("Preserve colour: one factor from the luminance for "
                 "all three channels. Per channel: each channel "
                 "stretched alone, which pales colours and can shift "
                 "them."))

/* (not "auto": that is a keyword of C) */
property_boolean (auto_points, _("Automatic black and symmetry point"), FALSE)
  description (_("Black point three robust sigmas below the median of "
                 "the image, and the symmetry point at the median; the "
                 "two settings are then not used."))

#else

#define GEGL_OP_FILTER
#define GEGL_OP_NAME     stretch
#define GEGL_OP_C_SOURCE stretch.c

#include "gegl-op.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define AUTO_SIDE    256
#define AUTO_SIGMAS  3.0
#define MAD_TO_SIGMA 1.4826

/* ---- GHS, ported from Siril's ght.c ----------------------------------- */

typedef struct
{
  gdouble B, D, LP, SP, HP;
  gboolean inverse;
  gdouble qlp, q0, qwp, q1, q, b1, a2, b2, c2, d2, e2, a3, b3, c3, d3, e3, a4, b4;
  gdouble LPT, SPT, HPT;
} Ghs;

/* GHTsetup (STRETCH_PAYNE_NORMAL and STRETCH_PAYNE_INVERSE); D here is
 * already e^stretch - 1 */
static void
ghs_setup (Ghs     *c,
           gdouble  B,
           gdouble  D,
           gdouble  LP,
           gdouble  SP,
           gdouble  HP,
           gboolean inverse)
{
  memset (c, 0, sizeof *c);
  c->B = B;
  c->D = D;
  c->LP = LP;
  c->SP = SP;
  c->HP = HP;
  c->inverse = inverse;
  if (D == 0.0)
    return;

  if (B == -1.0)
    {
      c->qlp = -log1p (D * (SP - LP));
      c->q0  = c->qlp - D * LP / (1.0 + D * (SP - LP));
      c->qwp = log1p (D * (HP - SP));
      c->q1  = c->qwp + D * (1.0 - HP) / (1.0 + D * (HP - SP));
      c->q   = 1.0 / (c->q1 - c->q0);
      if (! inverse)
        {
          c->b1 = (1.0 + D * (SP - LP)) / (D * c->q);
          c->a2 = -c->q0 * c->q;
          c->b2 = -c->q;
          c->c2 = 1.0 + D * SP;
          c->d2 = -D;
          c->a3 = -c->q0 * c->q;
          c->b3 = c->q;
          c->c3 = 1.0 - D * SP;
          c->d3 = D;
          c->a4 = (c->qwp - c->q0 - D * HP / (1.0 + D * (HP - SP))) * c->q;
          c->b4 = c->q * D / (1.0 + D * (HP - SP));
        }
      else
        {
          c->LPT = (c->qlp - c->q0) * c->q;
          /* (Siril has + q0 q, which is negative, so the lower branch
           * was never taken: the forward curve gives - q0 q at SP) */
          c->SPT = -c->q0 * c->q;
          c->HPT = (c->qwp - c->q0) * c->q;
          c->b1 = (1.0 + D * (SP - LP)) / (D * c->q);
          c->a2 = (1.0 + D * SP) / D;
          c->b2 = -1.0 / D;
          c->c2 = -c->q0;
          c->d2 = -1.0 / c->q;
          c->a3 = -(1.0 - D * SP) / D;
          c->b3 = 1.0 / D;
          c->c3 = c->q0;
          c->d3 = 1.0 / c->q;
          c->a4 = HP + (c->q0 - c->qwp) * (1.0 + D * (HP - SP)) / D;
          c->b4 = (1.0 + D * (HP - SP)) / (c->q * D);
        }
    }
  else if (B < 0.0)
    {
      B = -B;
      c->qlp = (1.0 - pow (1.0 + D * B * (SP - LP), (B - 1.0) / B)) / (B - 1.0);
      c->q0  = c->qlp - D * LP * pow (1.0 + D * B * (SP - LP), -1.0 / B);
      c->qwp = (pow (1.0 + D * B * (HP - SP), (B - 1.0) / B) - 1.0) / (B - 1.0);
      c->q1  = c->qwp + D * (1.0 - HP) * pow (1.0 + D * B * (HP - SP), -1.0 / B);
      c->q   = 1.0 / (c->q1 - c->q0);
      if (! inverse)
        {
          c->b1 = D * pow (1.0 + D * B * (SP - LP), -1.0 / B) * c->q;
          c->a2 = (1.0 / (B - 1.0) - c->q0) * c->q;
          c->b2 = -c->q / (B - 1.0);
          c->c2 = 1.0 + D * B * SP;
          c->d2 = -D * B;
          c->e2 = (B - 1.0) / B;
          c->a3 = (-1.0 / (B - 1.0) - c->q0) * c->q;
          c->b3 = c->q / (B - 1.0);
          c->c3 = 1.0 - D * B * SP;
          c->d3 = D * B;
          c->e3 = (B - 1.0) / B;
          c->a4 = (c->qwp - c->q0 - D * HP * pow (1.0 + D * B * (HP - SP), -1.0 / B)) * c->q;
          c->b4 = D * pow (1.0 + D * B * (HP - SP), -1.0 / B) * c->q;
        }
      else
        {
          c->LPT = (c->qlp - c->q0) * c->q;
          c->SPT = -c->q0 * c->q;
          c->HPT = (c->qwp - c->q0) * c->q;
          c->b1 = pow (1.0 + D * B * (SP - LP), 1.0 / B) / (c->q * D);
          c->a2 = (1.0 + D * B * SP) / (D * B);
          c->b2 = -1.0 / (D * B);
          c->c2 = -c->q0 * (B - 1.0) + 1.0;
          c->d2 = (1.0 - B) / c->q;
          c->e2 = B / (B - 1.0);
          c->a3 = (D * B * SP - 1.0) / (D * B);
          c->b3 = 1.0 / (D * B);
          c->c3 = 1.0 + c->q0 * (B - 1.0);
          c->d3 = (B - 1.0) / c->q;
          c->e3 = B / (B - 1.0);
          c->a4 = (c->q0 - c->qwp) / (D * pow (1.0 + D * B * (HP - SP), -1.0 / B)) + HP;
          c->b4 = 1.0 / (D * pow (1.0 + D * B * (HP - SP), -1.0 / B) * c->q);
        }
    }
  else if (B == 0.0)
    {
      c->qlp = exp (-D * (SP - LP));
      c->q0  = c->qlp - D * LP * exp (-D * (SP - LP));
      c->qwp = 2.0 - exp (-D * (HP - SP));
      c->q1  = c->qwp + D * (1.0 - HP) * exp (-D * (HP - SP));
      c->q   = 1.0 / (c->q1 - c->q0);
      if (! inverse)
        {
          c->b1 = D * exp (-D * (SP - LP)) * c->q;
          c->a2 = -c->q0 * c->q;
          c->b2 = c->q;
          c->c2 = -D * SP;
          c->d2 = D;
          c->a3 = (2.0 - c->q0) * c->q;
          c->b3 = -c->q;
          c->c3 = D * SP;
          c->d3 = -D;
          c->a4 = (c->qwp - c->q0 - D * HP * exp (-D * (HP - SP))) * c->q;
          c->b4 = D * exp (-D * (HP - SP)) * c->q;
        }
      else
        {
          c->LPT = (c->qlp - c->q0) * c->q;
          c->SPT = (1.0 - c->q0) * c->q;
          c->HPT = (c->qwp - c->q0) * c->q;
          c->b1 = 1.0 / (D * exp (-D * (SP - LP)) * c->q);
          c->a2 = SP;
          c->b2 = 1.0 / D;
          c->c2 = c->q0;
          c->d2 = 1.0 / c->q;
          c->a3 = SP;
          c->b3 = -1.0 / D;
          c->c3 = 2.0 - c->q0;
          c->d3 = -1.0 / c->q;
          c->a4 = (c->q0 - c->qwp) / (D * exp (-D * (HP - SP))) + HP;
          c->b4 = 1.0 / (D * exp (-D * (HP - SP)) * c->q);
        }
    }
  else /* B > 0 */
    {
      c->qlp = pow (1.0 + D * B * (SP - LP), -1.0 / B);
      c->q0  = c->qlp - D * LP * pow (1.0 + D * B * (SP - LP), -(1.0 + B) / B);
      c->qwp = 2.0 - pow (1.0 + D * B * (HP - SP), -1.0 / B);
      c->q1  = c->qwp + D * (1.0 - HP) * pow (1.0 + D * B * (HP - SP), -(1.0 + B) / B);
      c->q   = 1.0 / (c->q1 - c->q0);
      if (! inverse)
        {
          c->b1 = D * pow (1.0 + D * B * (SP - LP), -(1.0 + B) / B) * c->q;
          c->a2 = -c->q0 * c->q;
          c->b2 = c->q;
          c->c2 = 1.0 + D * B * SP;
          c->d2 = -D * B;
          c->e2 = -1.0 / B;
          c->a3 = (2.0 - c->q0) * c->q;
          c->b3 = -c->q;
          c->c3 = 1.0 - D * B * SP;
          c->d3 = D * B;
          c->e3 = -1.0 / B;
          c->a4 = (c->qwp - c->q0 - D * HP * pow (1.0 + D * B * (HP - SP), -(B + 1.0) / B)) * c->q;
          c->b4 = D * pow (1.0 + D * B * (HP - SP), -(B + 1.0) / B) * c->q;
        }
      else
        {
          c->LPT = (c->qlp - c->q0) * c->q;
          c->SPT = (1.0 - c->q0) * c->q;
          c->HPT = (c->qwp - c->q0) * c->q;
          c->b1 = 1.0 / (D * pow (1.0 + D * B * (SP - LP), -(1.0 + B) / B) * c->q);
          c->a2 = 1.0 / (D * B) + SP;
          c->b2 = -1.0 / (D * B);
          c->c2 = c->q0;
          c->d2 = 1.0 / c->q;
          c->e2 = -B;
          c->a3 = -1.0 / (D * B) + SP;
          c->b3 = 1.0 / (D * B);
          c->c3 = 2.0 - c->q0;
          c->d3 = -1.0 / c->q;
          c->e3 = -B;
          c->a4 = (c->q0 - c->qwp) / (D * pow (1.0 + D * B * (HP - SP), -(B + 1.0) / B)) + HP;
          c->b4 = 1.0 / (D * pow (1.0 + D * B * (HP - SP), -(B + 1.0) / B) * c->q);
        }
    }
}

/* GHT for x in 0 to 1 (black point already taken off) */
static inline gdouble
ghs_eval (const Ghs *c,
          gdouble    x)
{
  if (c->D == 0.0)
    return x;
  if (! c->inverse)
    {
      if (x < c->LP)
        return c->b1 * x;
      if (x >= c->HP)
        return c->a4 + c->b4 * x;
      if (x < c->SP)
        {
          if (c->B == -1.0)
            return c->a2 + c->b2 * log (c->c2 + c->d2 * x);
          if (c->B == 0.0)
            return c->a2 + c->b2 * exp (c->c2 + c->d2 * x);
          return c->a2 + c->b2 * pow (c->c2 + c->d2 * x, c->e2);
        }
      if (c->B == -1.0)
        return c->a3 + c->b3 * log (c->c3 + c->d3 * x);
      if (c->B == 0.0)
        return c->a3 + c->b3 * exp (c->c3 + c->d3 * x);
      return c->a3 + c->b3 * pow (c->c3 + c->d3 * x, c->e3);
    }
  if (x < c->LPT)
    return c->b1 * x;
  if (x >= c->HPT)
    return c->a4 + c->b4 * x;
  if (x < c->SPT)
    {
      if (c->B == -1.0)
        return c->a2 + c->b2 * exp (c->c2 + c->d2 * x);
      if (c->B == 0.0)
        return c->a2 + c->b2 * log (c->c2 + c->d2 * x);
      return c->a2 + c->b2 * pow (c->c2 + c->d2 * x, c->e2);
    }
  if (c->B == -1.0)
    return c->a3 + c->b3 * exp (c->c3 + c->d3 * x);
  if (c->B == 0.0)
    return c->a3 + c->b3 * log (c->c3 + c->d3 * x);
  return c->a3 + c->b3 * pow (c->c3 + c->d3 * x, c->e3);
}

/* ---- the operation ---------------------------------------------------- */

typedef struct
{
  gint    curve;
  gdouble bp[3];      /* black point per channel (all the same unless auto per channel) */
  Ghs     ghs[3];     /* per channel (differ only with auto per channel) */
  gdouble beta, asinh_beta;
} Curve;

static inline gdouble
curve_eval (const Curve *k,
            gint         c,
            gdouble      x)
{
  if (k->curve == ASTRO_STRETCH_ARCSINH)
    return k->beta == 0.0 ? x : asinh (k->beta * x) / k->asinh_beta;
  return ghs_eval (&k->ghs[c], x);
}

static inline gdouble
black (gdouble x,
       gdouble bp)
{
  x = (x - bp) / (1.0 - bp);
  return x < 0.0 ? 0.0 : x > 1.0 ? 1.0 : x;
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
whole_or_roi (GeglOperation       *operation,
              const GeglRectangle *roi)
{
  GeglProperties      *o  = GEGL_PROPERTIES (operation);
  const GeglRectangle *in = gegl_operation_source_get_bounding_box (operation, "input");

  if (o->auto_points && in && ! gegl_rectangle_is_infinite_plane (in))
    return *in;
  return *roi;
}

static GeglRectangle
get_required_for_output (GeglOperation       *operation,
                         const gchar         *input_pad,
                         const GeglRectangle *roi)
{
  (void) input_pad;
  return whole_or_roi (operation, roi);
}

static GeglRectangle
get_invalidated_by_change (GeglOperation       *operation,
                           const gchar         *input_pad,
                           const GeglRectangle *roi)
{
  (void) input_pad;
  return whole_or_roi (operation, roi);
}

static int
cmp_float (const void *a,
           const void *b)
{
  gfloat x = *(const gfloat *) a, y = *(const gfloat *) b;

  return (x > y) - (x < y);
}

/* median and robust sigma of v[0..n-1] (reordered) */
static void
robust_stats (gfloat  *v,
              gint     n,
              gdouble *median,
              gdouble *sigma)
{
  gfloat med;
  gint   i;

  qsort (v, n, sizeof (gfloat), cmp_float);
  med = n % 2 ? v[n / 2] : 0.5f * (v[n / 2 - 1] + v[n / 2]);
  for (i = 0; i < n; i++)
    v[i] = fabsf (v[i] - med);
  qsort (v, n, sizeof (gfloat), cmp_float);
  *median = med;
  *sigma  = MAD_TO_SIGMA * (n % 2 ? v[n / 2] : 0.5f * (v[n / 2 - 1] + v[n / 2]));
}

/* the automatic black point and symmetry point, per channel or from the
 * luminance, on a fixed lattice of the whole input (opaque pixels only) */
static gboolean
auto_points (GeglBuffer          *input,
             const GeglRectangle *box,
             const Babl          *format,
             gboolean             per_channel,
             const gdouble        lum[3],
             gdouble              bp[3],
             gdouble              sp[3])
{
  gint    step = MAX (1, (MAX (box->width, box->height) + AUTO_SIDE - 1) / AUTO_SIDE);
  gint    nx = (box->width - step / 2 + step - 1) / step;
  gint    ny = (box->height - step / 2 + step - 1) / step;
  gfloat *row = g_new (gfloat, (gsize) box->width * 4);
  gfloat *vals[3];
  gint    n = 0, x, y, c;

  for (c = 0; c < 3; c++)
    vals[c] = g_new (gfloat, (gsize) MAX (1, nx * ny));
  for (y = step / 2; y < box->height; y += step)
    {
      gegl_buffer_get (input, GEGL_RECTANGLE (box->x, box->y + y, box->width, 1),
                       1.0, format, row, GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);
      for (x = step / 2; x < box->width; x += step)
        {
          const gfloat *p = row + (gsize) x * 4;

          if (! (p[3] > 0.5f) || ! isfinite (p[0]) || ! isfinite (p[1]) || ! isfinite (p[2]))
            continue;
          if (per_channel)
            for (c = 0; c < 3; c++)
              vals[c][n] = p[c];
          else
            vals[0][n] = lum[0] * p[0] + lum[1] * p[1] + lum[2] * p[2];
          n++;
        }
    }
  g_free (row);
  if (n >= 16)
    for (c = 0; c < (per_channel ? 3 : 1); c++)
      {
        gdouble med, sigma;

        robust_stats (vals[c], n, &med, &sigma);
        bp[c] = CLAMP (med - AUTO_SIGMAS * sigma, 0.0, 0.99);
        sp[c] = CLAMP ((med - bp[c]) / (1.0 - bp[c]), 0.0, 1.0);
      }
  if (n >= 16 && ! per_channel)
    for (c = 1; c < 3; c++)
      {
        bp[c] = bp[0];
        sp[c] = sp[0];
      }
  for (c = 0; c < 3; c++)
    g_free (vals[c]);
  return n >= 16;
}

static gboolean
process (GeglOperation       *operation,
         GeglBuffer          *input,
         GeglBuffer          *output,
         const GeglRectangle *result,
         gint                 level)
{
  GeglProperties     *o       = GEGL_PROPERTIES (operation);
  const Babl         *format  = gegl_operation_get_format (operation, "output");
  gboolean            per_ch  = o->colour == ASTRO_STRETCH_PER_CHANNEL;
  gdouble             lum[3]  = { 0.2126, 0.7152, 0.0722 };
  gdouble             bp[3], sp[3];
  Curve               k;
  GeglBufferIterator *it;
  gint                c;

  (void) level;
  babl_space_get_rgb_luminance (babl_format_get_space (format), &lum[0], &lum[1], &lum[2]);
  for (c = 0; c < 3; c++)
    {
      bp[c] = o->black_point;
      sp[c] = o->symmetry_point;
    }
  if (o->auto_points)
    {
      const GeglRectangle *box = gegl_operation_source_get_bounding_box (operation, "input");

      if (box && box->width > 0 && box->height > 0)
        auto_points (input, box, format, per_ch, lum, bp, sp);
    }

  k.curve      = o->curve;
  k.beta       = o->beta;
  k.asinh_beta = asinh (o->beta);
  for (c = 0; c < 3; c++)
    {
      gdouble SP = sp[c];
      gdouble LP = MIN (o->protect_shadows, SP);
      gdouble HP = MAX (o->protect_highlights, SP);

      k.bp[c] = bp[c];
      ghs_setup (&k.ghs[c], o->local_intensity, expm1 (o->stretch), LP, SP, HP, o->inverse);
    }

  it = gegl_buffer_iterator_new (output, result, 0, format,
                                 GEGL_ACCESS_WRITE, GEGL_ABYSS_NONE, 2);
  gegl_buffer_iterator_add (it, input, result, 0, format,
                            GEGL_ACCESS_READ, GEGL_ABYSS_NONE);
  while (gegl_buffer_iterator_next (it))
    {
      gfloat       *out = it->items[0].data;
      const gfloat *in  = it->items[1].data;
      glong         i, n = it->length;

      for (i = 0; i < n; i++, in += 4, out += 4)
        {
          gdouble x[3];

          for (c = 0; c < 3; c++)
            x[c] = black (isfinite (in[c]) ? in[c] : 0.0, k.bp[c]);
          if (per_ch)
            {
              for (c = 0; c < 3; c++)
                out[c] = curve_eval (&k, c, x[c]);
            }
          else
            {
              gdouble Y = lum[0] * x[0] + lum[1] * x[1] + lum[2] * x[2];
              gdouble f, m;

              if (Y <= 0.0)
                {
                  /* black stays black (the curve maps 0 to 0) */
                  for (c = 0; c < 3; c++)
                    out[c] = 0.0f;
                }
              else
                {
                  f = curve_eval (&k, 0, MIN (Y, 1.0)) / Y;
                  m = MAX (x[0], MAX (x[1], x[2])) * f;
                  if (m > 1.0)
                    f /= m;
                  for (c = 0; c < 3; c++)
                    out[c] = x[c] * f;
                }
            }
          out[3] = in[3];
        }
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
  operation_class->get_invalidated_by_change = get_invalidated_by_change;
  filter_class->process                      = process;

  gegl_operation_class_set_keys (operation_class,
    "name",            "astro:stretch",
    "title",           _("Astro Stretch"),
    "categories",      "color",
    "description",     _("Brings out the faint parts of linear astronomical "
                         "data: the generalised hyperbolic stretch (GHS, "
                         "after Siril) or the arcsinh stretch."),
    "gimp:menu-path",  "<Image>/Filters/Astronomy",
    "gimp:menu-label", _("Stretch..."),
    NULL);
}

#endif
