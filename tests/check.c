/*
 * Numerical checks of the operations on synthetic star fields with known
 * stars (Gaussian PSFs), a known sky gradient and known noise: how well
 * the operations recover what was put in. Prints PASS or FAIL per check
 * with the numbers, and "N passed, M failed".
 *
 *   check <build folder>
 *
 * check.c
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
 */

#include <gegl.h>
#include <glib/gstdio.h>
#include <locale.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define W 600
#define H 400
#define FMT "RGBA float"

static const gchar *const modules[] = { "astro-background", "astro-stretch" };

/* the operations, for the checks that apply to each */
static const gchar *const ops[] = { "astro:background", "astro:stretch" };

/* the build folder, where checks may leave files for tests/check.sh */
static const gchar *build_dir;

static gint n_passed, n_failed, n_messages;

static void
report (const gchar *name,
        gboolean     ok,
        const gchar *fmt,
        ...)
{
  va_list ap;

  printf ("%s  %s: ", ok ? "PASS" : "FAIL", name);
  va_start (ap, fmt);
  vprintf (fmt, ap);
  va_end (ap);
  printf ("\n");
  if (ok)
    n_passed++;
  else
    n_failed++;
}

static void
log_handler (const gchar    *domain,
             GLogLevelFlags  level,
             const gchar    *message,
             gpointer        data)
{
  (void) data;
  if (level & (G_LOG_LEVEL_WARNING | G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_ERROR))
    {
      n_messages++;
      fprintf (stderr, "message from %s: %s\n", domain ? domain : "?", message);
    }
}

/* Links the modules of the build folder into a folder of their own and
 * points GEGL_PATH there and at GEGL's own operations only, so installed
 * copies are not loaded. */
static gchar *
modules_setup (const gchar *build)
{
  gchar *dir = g_dir_make_tmp ("astro-modules-XXXXXX", NULL);
  gchar *path;
  guint  i;

  if (! dir)
    return NULL;
  for (i = 0; i < G_N_ELEMENTS (modules); i++)
    {
      gchar *name   = g_strconcat (modules[i], ".so", NULL);
      gchar *from   = g_build_filename (build, name, NULL);
      gchar *target = g_canonicalize_filename (from, NULL);
      gchar *link   = g_build_filename (dir, name, NULL);

      if (symlink (target, link) != 0)
        {
          fprintf (stderr, "cannot link %s\n", target);
          return NULL;
        }
      g_free (name);
      g_free (from);
      g_free (target);
      g_free (link);
    }
  path = g_strconcat (dir, G_SEARCHPATH_SEPARATOR_S, GEGL_PLUGINSDIR, NULL);
  g_setenv ("GEGL_PATH", path, TRUE);
  g_free (path);
  return dir;
}

static void
modules_cleanup (gchar *dir)
{
  guint i;

  for (i = 0; i < G_N_ELEMENTS (modules); i++)
    {
      gchar *name = g_strconcat (modules[i], ".so", NULL);
      gchar *link = g_build_filename (dir, name, NULL);

      g_unlink (link);
      g_free (name);
      g_free (link);
    }
  g_rmdir (dir);
  g_free (dir);
}

/* ---- the synthetic field ---------------------------------------------- */

typedef struct
{
  gfloat *sky;     /* the gradient alone, RGB */
  gfloat *stars;   /* the stars alone, RGB */
  gfloat *noise;   /* the noise alone, RGB */
  gfloat *image;   /* RGBA: sky + stars + noise (or times, for vignetting) */
} Field;

static gdouble
gauss (GRand *rnd)
{
  /* Box-Muller */
  gdouble a = g_rand_double_range (rnd, 1e-12, 1.0), b = g_rand_double (rnd);

  return sqrt (-2.0 * log (a)) * cos (2.0 * G_PI * b);
}

/* the gradient: a quadratic in u, v (-1 to 1), a little different per
 * channel, like light pollution that is redder low in the sky */
static void
sky_quadratic (gdouble u, gdouble v, gfloat *rgb)
{
  rgb[0] = 0.060 + 0.030 * u + 0.025 * v + 0.012 * u * u - 0.008 * u * v + 0.010 * v * v;
  rgb[1] = 0.050 + 0.025 * u + 0.020 * v + 0.010 * u * u - 0.006 * u * v + 0.008 * v * v;
  rgb[2] = 0.045 + 0.020 * u + 0.012 * v + 0.008 * u * u - 0.004 * u * v + 0.006 * v * v;
}

/* vignetting: 1 in the centre, falls off with r^2 */
static void
sky_vignette (gdouble u, gdouble v, gfloat *rgb)
{
  gdouble f = 1.0 - 0.35 * (u * u + v * v) / 2.0;
  gint    c;

  for (c = 0; c < 3; c++)
    rgb[c] = 0.08 * f;
}

/* n Gaussian stars with sigma 0.8 to 2.5 px and peaks spread like real
 * star counts (many faint, few bright), colours from red to blue; a
 * nebula (a wide Gaussian) when nebula > 0; noise of sigma noise_sigma.
 * Multiplicative sky (vignetting) multiplies stars and nebula too. */
static Field
make_field (guint32  seed,
            gint     n_stars,
            void   (*sky) (gdouble, gdouble, gfloat *),
            gboolean multiply,
            gdouble  nebula,
            gdouble  noise_sigma)
{
  Field  f;
  GRand *rnd = g_rand_new_with_seed (seed);
  gint   i, x, y, c;

  f.sky   = g_new0 (gfloat, W * H * 3);
  f.stars = g_new0 (gfloat, W * H * 3);
  f.noise = g_new0 (gfloat, W * H * 3);
  f.image = g_new0 (gfloat, W * H * 4);

  for (i = 0; i < n_stars; i++)
    {
      gdouble sx = g_rand_double_range (rnd, 0, W), sy = g_rand_double_range (rnd, 0, H);
      gdouble sigma = g_rand_double_range (rnd, 0.8, 2.5);
      gdouble peak = 0.02 * pow (g_rand_double_range (rnd, 1e-3, 1.0), -0.7);
      gdouble t = g_rand_double (rnd);
      gdouble col[3] = { 0.8 + 0.4 * t, 1.0, 1.2 - 0.4 * t };
      gint    rad = (gint) ceil (4 * sigma);

      peak = MIN (peak, 0.9);
      for (y = MAX (0, (gint) sy - rad); y <= MIN (H - 1, (gint) sy + rad); y++)
        for (x = MAX (0, (gint) sx - rad); x <= MIN (W - 1, (gint) sx + rad); x++)
          {
            gdouble dx = x + 0.5 - sx, dy = y + 0.5 - sy;
            gdouble g = peak * exp (-(dx * dx + dy * dy) / (2 * sigma * sigma));

            for (c = 0; c < 3; c++)
              f.stars[(y * W + x) * 3 + c] += g * col[c];
          }
    }
  if (nebula > 0)
    for (y = 0; y < H; y++)
      for (x = 0; x < W; x++)
        {
          gdouble dx = x - 0.62 * W, dy = y - 0.45 * H;
          gdouble g = nebula * exp (-(dx * dx + dy * dy) / (2 * 55.0 * 55.0));

          f.stars[(y * W + x) * 3 + 0] += g;
          f.stars[(y * W + x) * 3 + 1] += 0.4 * g;
          f.stars[(y * W + x) * 3 + 2] += 0.6 * g;
        }
  for (y = 0; y < H; y++)
    for (x = 0; x < W; x++)
      {
        gdouble u = 2.0 * (x + 0.5) / W - 1.0, v = 2.0 * (y + 0.5) / H - 1.0;
        gfloat  s[3];

        sky (u, v, s);
        for (c = 0; c < 3; c++)
          {
            gint   k = (y * W + x) * 3 + c;
            gfloat n = noise_sigma * gauss (rnd);

            f.sky[k]   = s[c];
            f.noise[k] = n;
            if (multiply)
              f.image[(y * W + x) * 4 + c] = s[c] / 0.08f * (0.08f + f.stars[k]) + n;
            else
              f.image[(y * W + x) * 4 + c] = s[c] + f.stars[k] + n;
          }
        f.image[(y * W + x) * 4 + 3] = 1.0f;
      }
  g_rand_free (rnd);
  return f;
}

static void
free_field (Field *f)
{
  g_free (f->sky);
  g_free (f->stars);
  g_free (f->noise);
  g_free (f->image);
}

/* runs op on the pixels (w x h RGBA float at x0, y0) with the
 * properties as name, value pairs (NULL ended), into a new array */
static gfloat *
run_valist (const gfloat *pixels,
            gint          w,
            gint          h,
            gint          x0,
            gint          y0,
            const gchar  *op,
            va_list       ap)
{
  GeglRectangle r = { x0, y0, w, h };
  const Babl   *fmt = babl_format (FMT);
  GeglBuffer   *in = gegl_buffer_new (&r, fmt);
  GeglBuffer   *out = NULL;
  GeglNode     *g = gegl_node_new ();
  GeglNode     *src, *node, *sink;
  gfloat       *res = g_new (gfloat, (gsize) w * h * 4);

  gegl_buffer_set (in, &r, 0, fmt, pixels, GEGL_AUTO_ROWSTRIDE);
  src  = gegl_node_new_child (g, "operation", "gegl:buffer-source", "buffer", in, NULL);
  node = gegl_node_new_child (g, "operation", op, NULL);
  gegl_node_set_valist (node, va_arg (ap, const gchar *), ap);
  sink = gegl_node_new_child (g, "operation", "gegl:buffer-sink", "buffer", &out, NULL);
  gegl_node_link_many (src, node, sink, NULL);
  gegl_node_process (sink);
  gegl_buffer_get (out, &r, 1.0, fmt, res, GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);
  g_object_unref (out);
  g_object_unref (g);
  g_object_unref (in);
  return res;
}

/* W x H at x0, y0 */
static gfloat *
run (const gfloat *pixels,
     gint          x0,
     gint          y0,
     const gchar  *op,
     ...)
{
  va_list ap;
  gfloat *res;

  va_start (ap, op);
  res = run_valist (pixels, W, H, x0, y0, op, ap);
  va_end (ap);
  return res;
}

/* any size, at 0, 0 */
static gfloat *
run_size (const gfloat *pixels,
          gint          w,
          gint          h,
          const gchar  *op,
          ...)
{
  va_list ap;
  gfloat *res;

  va_start (ap, op);
  res = run_valist (pixels, w, h, 0, 0, op, ap);
  va_end (ap);
  return res;
}

/* amplitude (max - min) of the gradient, over all channels */
static gdouble
amplitude (const gfloat *sky)
{
  gdouble lo = G_MAXDOUBLE, hi = -G_MAXDOUBLE;
  gint    i;

  for (i = 0; i < W * H * 3; i++)
    {
      lo = MIN (lo, sky[i]);
      hi = MAX (hi, sky[i]);
    }
  return hi - lo;
}

/* ---- the checks ------------------------------------------------------- */

static void
check_registration (void)
{
  guint k;

  for (k = 0; k < G_N_ELEMENTS (ops); k++)
    {
      const gchar *op = ops[k];
      const gchar *base = strchr (op, ':') + 1;
      gboolean     ok = gegl_has_operation (op);
      const gchar *path = ok ? gegl_operation_get_key (op, "gimp:menu-path") : NULL;
      gchar       *name = g_strdup_printf ("%s_registered", base);

      report (name, ok && g_strcmp0 (path, "<Image>/Filters/Astronomy") == 0,
              "operation %s, menu %s", ok ? "found" : "missing", path ? path : "none");
      g_free (name);

      /* no property named like a pad: GIMP's dialog leaves such a one out */
      if (ok)
        {
          guint        n, i;
          GParamSpec **props = gegl_operation_list_properties (op, &n);
          const gchar *clash = NULL;

          for (i = 0; i < n; i++)
            if (! strcmp (props[i]->name, "input") || ! strcmp (props[i]->name, "output") ||
                ! strcmp (props[i]->name, "aux"))
              clash = props[i]->name;
          name = g_strdup_printf ("%s_no_property_named_like_a_pad", base);
          report (name, clash == NULL, "%u properties%s%s",
                  n, clash ? ", clash: " : "", clash ? clash : "");
          g_free (name);
          g_free (props);
        }
    }
}

/* a sky without gradient stays as it is */
static void
check_flat (void)
{
  gfloat *p = g_new (gfloat, W * H * 4);
  gint    i, k;

  for (i = 0; i < W * H; i++)
    {
      p[i * 4 + 0] = 0.1f;
      p[i * 4 + 1] = 0.2f;
      p[i * 4 + 2] = 0.3f;
      p[i * 4 + 3] = 1.0f;
    }
  for (k = 0; k < 2; k++)
    {
      gfloat *o = run (p, 0, 0, "astro:background", "correction", k, NULL);
      gdouble d = 0;

      for (i = 0; i < W * H * 4; i++)
        d = MAX (d, fabs (o[i] - p[i]));
      report (k ? "background_flat_divide_unchanged" : "background_flat_subtract_unchanged",
              d < 1e-5, "max difference %.2g", d);
      g_free (o);
    }
  g_free (p);
}

/* the model against the true sky: max error over the image, as a share
 * of the gradient's amplitude */
static gdouble
model_error (const Field *f,
             const gfloat *model)
{
  gdouble e = 0;
  gint    i, c;

  for (i = 0; i < W * H; i++)
    for (c = 0; c < 3; c++)
      e = MAX (e, fabs (model[i * 4 + c] - f->sky[i * 3 + c]));
  return e / amplitude (f->sky);
}

static void
check_gradient_recovery (void)
{
  Field   f = make_field (1, 600, sky_quadratic, FALSE, 0.0, 0.003);
  gfloat *model = run (f.image, 0, 0, "astro:background", "show", 1, NULL);
  gfloat *corr  = run (f.image, 0, 0, "astro:background", NULL);
  gdouble e = model_error (&f, model), lo = G_MAXDOUBLE, hi = -G_MAXDOUBLE, flux_in = 0, flux_out = 0;
  gdouble level[3] = { 0, 0, 0 };
  gint    i, c;

  report ("background_quadratic_model", e < 0.02,
          "600 stars, noise 0.003: model within %.2f %% of the gradient's amplitude "
          "(%.3f)", 100 * e, amplitude (f.sky));

  /* what is left of the sky after the correction: corrected minus stars
   * and noise, which should be flat (the mean level) */
  for (i = 0; i < W * H; i++)
    for (c = 0; c < 3; c++)
      level[c] += f.sky[i * 3 + c] / (W * H);
  for (i = 0; i < W * H; i++)
    for (c = 0; c < 3; c++)
      {
        gdouble left = corr[i * 4 + c] - f.stars[i * 3 + c] - f.noise[i * 3 + c] - level[c];

        lo = MIN (lo, left);
        hi = MAX (hi, left);
        flux_in  += f.stars[i * 3 + c];
        flux_out += corr[i * 4 + c] - f.noise[i * 3 + c] - level[c];
      }
  report ("background_quadratic_flat_after", (hi - lo) / amplitude (f.sky) < 0.04,
          "sky left after the correction spans %.2f %% of the gradient "
          "(from %.2f %% to %.2f %% around the mean level)",
          100 * (hi - lo) / amplitude (f.sky), 100 * lo / amplitude (f.sky),
          100 * hi / amplitude (f.sky));
  report ("background_star_flux_kept", fabs (flux_out / flux_in - 1) < 0.01,
          "star flux after / before %.4f", flux_out / flux_in);
  g_free (model);
  g_free (corr);
  free_field (&f);
}

/* a plane with degree 1: exact up to noise */
static void
check_degrees (void)
{
  gint d;

  for (d = 1; d <= 4; d++)
    {
      Field   f = make_field (2, 400, sky_quadratic, FALSE, 0.0, 0.003);
      gfloat *model = run (f.image, 0, 0, "astro:background", "show", 1, "degree", d, NULL);
      gdouble e = model_error (&f, model);
      /* degree 1 cannot follow a quadratic; 2 and up can */
      gdouble limit = d == 1 ? 0.25 : 0.03;
      gchar  *name = g_strdup_printf ("background_degree_%d", d);

      report (name, e < limit, "model within %.2f %% of the amplitude (limit %.0f %%)",
              100 * e, 100 * limit);
      g_free (name);
      g_free (model);
      free_field (&f);
    }
}

/* a nebula of a fifth of the image: its cells are left out */
static void
check_nebula (void)
{
  Field   f = make_field (3, 500, sky_quadratic, FALSE, 0.06, 0.003);
  gfloat *model = run (f.image, 0, 0, "astro:background", "show", 1, NULL);
  gdouble e = model_error (&f, model);

  report ("background_nebula_left_out", e < 0.06,
          "nebula of peak 0.06 (about the gradient's amplitude): model within %.2f %%", 100 * e);
  g_free (model);
  free_field (&f);
}

static void
check_vignetting (void)
{
  Field   f = make_field (4, 500, sky_vignette, TRUE, 0.0, 0.002);
  gfloat *corr = run (f.image, 0, 0, "astro:background", "correction", 1, NULL);
  gdouble lo = G_MAXDOUBLE, hi = -G_MAXDOUBLE, mean = 0, level = 0;
  gint    i, n = 0;

  for (i = 0; i < W * H; i++)
    level += f.sky[i * 3 + 1] / (W * H);
  /* the sky after the division, where there are no stars: flat (the
   * noise was divided by the sky and multiplied by the level too) */
  for (i = 0; i < W * H; i++)
    if (f.stars[i * 3 + 1] < 1e-4)
      {
        gdouble s = corr[i * 4 + 1] - f.noise[i * 3 + 1] * level / f.sky[i * 3 + 1];

        lo = MIN (lo, s);
        hi = MAX (hi, s);
        mean += s;
        n++;
      }
  mean /= n;
  report ("background_vignetting_divided", (hi - lo) / mean < 0.03,
          "sky after dividing out 35 %% vignetting: spans %.2f %% of its level",
          100 * (hi - lo) / mean);
  g_free (corr);
  free_field (&f);
}

/* the same pixels at another place give the same result, bit for bit;
 * and two runs are the same */
static void
check_determinism (void)
{
  Field   f = make_field (5, 400, sky_quadratic, FALSE, 0.02, 0.003);
  gfloat *a = run (f.image, 0, 0, "astro:background", NULL);
  gfloat *b = run (f.image, 0, 0, "astro:background", NULL);
  gfloat *c = run (f.image, 137, -61, "astro:background", NULL);

  report ("background_deterministic", memcmp (a, b, sizeof (gfloat) * W * H * 4) == 0,
          "two runs %s", memcmp (a, b, sizeof (gfloat) * W * H * 4) ? "differ" : "identical");
  report ("background_offset_invariant", memcmp (a, c, sizeof (gfloat) * W * H * 4) == 0,
          "buffer at 137, -61 %s", memcmp (a, c, sizeof (gfloat) * W * H * 4) ? "differs" : "identical");
  g_free (a);
  g_free (b);
  g_free (c);
  free_field (&f);
}

/* transparent pixels do not count, and alpha is kept */
static void
check_alpha (void)
{
  Field   f = make_field (6, 400, sky_quadratic, FALSE, 0.0, 0.003);
  gfloat *model, *corr;
  gdouble e, da = 0;
  gint    x, y, i;

  /* a transparent block full of junk */
  for (y = 50; y < 200; y++)
    for (x = 40; x < 260; x++)
      {
        gfloat *p = f.image + (y * W + x) * 4;

        p[0] = p[1] = p[2] = 5.0f;
        p[3] = 0.0f;
      }
  for (i = 0; i < W * H; i++)
    if (i % 7 == 0 && f.image[i * 4 + 3] > 0)
      f.image[i * 4 + 3] = 0.6f;
  model = run (f.image, 0, 0, "astro:background", "show", 1, NULL);
  corr  = run (f.image, 0, 0, "astro:background", NULL);
  e     = model_error (&f, model);
  for (i = 0; i < W * H; i++)
    da = MAX (da, fabs (corr[i * 4 + 3] - f.image[i * 4 + 3]));
  report ("background_transparent_ignored", e < 0.03,
          "transparent block of junk: model within %.2f %%", 100 * e);
  report ("background_alpha_kept", da == 0, "max alpha difference %.2g", da);
  g_free (model);
  g_free (corr);
  free_field (&f);
}

/* ---- astro:stretch ------------------------------------------------------ */

typedef struct
{
  gboolean arcsinh;
  gdouble  b, stretch, lp, sp, hp, bp, beta;
  gboolean inverse;
} StretchCase;

static const StretchCase stretch_cases[] =
  {
    /* GHS: b, stretch, LP, SP, HP, BP; inverse */
    { FALSE,  0.0,  0.0, 0.0,  0.0,  1.0,  0.0,  0, FALSE },   /* identity */
    { FALSE, -1.0,  3.0, 0.0,  0.1,  1.0,  0.0,  0, FALSE },
    { FALSE, -5.0,  2.0, 0.02, 0.2,  0.9,  0.0,  0, FALSE },
    { FALSE, -2.5,  4.0, 0.0,  0.05, 0.8,  0.0,  0, FALSE },
    { FALSE, -0.5,  5.0, 0.01, 0.1,  1.0,  0.0,  0, FALSE },
    { FALSE,  0.0,  5.0, 0.0,  0.0,  1.0,  0.0,  0, FALSE },
    { FALSE,  0.0,  3.0, 0.05, 0.25, 0.7,  0.0,  0, FALSE },
    { FALSE,  0.5,  6.0, 0.0,  0.02, 1.0,  0.0,  0, FALSE },
    { FALSE,  1.0,  5.0, 0.01, 0.08, 0.95, 0.0,  0, FALSE },
    { FALSE,  3.0,  8.0, 0.0,  0.1,  0.9,  0.0,  0, FALSE },
    { FALSE, 15.0, 10.0, 0.0,  0.05, 1.0,  0.0,  0, FALSE },
    { FALSE,  2.0,  4.0, 0.3,  0.3,  0.3,  0.0,  0, FALSE },   /* LP = SP = HP */
    { FALSE,  6.0,  7.0, 0.0,  1.0,  1.0,  0.0,  0, FALSE },   /* SP at 1 */
    { FALSE,  1.0,  3.0, 0.4,  0.2,  0.1,  0.0,  0, FALSE },   /* LP, HP out of order */
    { FALSE,  1.0,  4.0, 0.0,  0.05, 1.0,  0.02, 0, FALSE },   /* black point */
    { FALSE, -1.0,  3.0, 0.0,  0.1,  1.0,  0.0,  0, TRUE },
    { FALSE, -2.5,  4.0, 0.0,  0.05, 0.8,  0.0,  0, TRUE },
    { FALSE,  0.0,  3.0, 0.05, 0.25, 0.7,  0.0,  0, TRUE },
    { FALSE,  1.0,  5.0, 0.01, 0.08, 0.95, 0.0,  0, TRUE },
    { FALSE,  3.0,  8.0, 0.0,  0.1,  0.9,  0.0,  0, TRUE },
    /* arcsinh: beta, BP */
    { TRUE, 0, 0, 0, 0, 0, 0.0,  0.0,    FALSE },
    { TRUE, 0, 0, 0, 0, 0, 0.0,  10.0,   FALSE },
    { TRUE, 0, 0, 0, 0, 0, 0.0,  100.0,  FALSE },
    { TRUE, 0, 0, 0, 0, 0, 0.03, 1000.0, FALSE },
  };

#define RAMP_N 5120

/* runs a stretch case per channel on a grey ramp; returns the outputs
 * (channel 0) and fills x with the inputs */
static gfloat *
stretch_ramp (const StretchCase *k,
              gfloat            *x,
              gboolean           inverse)
{
  gfloat *p = g_new (gfloat, RAMP_N * 4), *o, *y = g_new (gfloat, RAMP_N);
  gdouble sp = MIN (k->sp, 0.99);
  gint    i;

  for (i = 0; i < RAMP_N; i++)
    {
      /* 4096 even steps from 0 to 1, then 1024 around the symmetry point */
      gfloat v = i < 4096 ? i / 4095.0f : (gfloat) (sp + (i - 4096 - 512) / 512.0 * 0.01);

      v = CLAMP (v, 0.0f, 1.0f);
      x[i] = v;
      p[i * 4 + 0] = p[i * 4 + 1] = p[i * 4 + 2] = v;
      p[i * 4 + 3] = 1.0f;
    }
  o = run_size (p, RAMP_N, 1, "astro:stretch", "curve", k->arcsinh ? 1 : 0,
                "colour", 1, "black-point", k->bp, "stretch", k->stretch,
                "local-intensity", k->b, "symmetry-point", k->sp,
                "protect-shadows", k->lp, "protect-highlights", k->hp,
                "inverse", inverse, "beta", k->beta, NULL);
  for (i = 0; i < RAMP_N; i++)
    y[i] = o[i * 4];
  g_free (o);
  g_free (p);
  return y;
}

/* every case on a ramp: written for tests/stretch-reference.py (the
 * published equations, written anew in Python) and checked here for 0 and
 * 1 kept and a rising curve */
static void
check_stretch_curves (void)
{
  gchar  *path = g_build_filename (build_dir, "stretch-op.txt", NULL);
  FILE   *f = fopen (path, "w");
  gfloat *x = g_new (gfloat, RAMP_N);
  guint   k;
  gint    i, bad_mono = 0, bad_ends = 0;
  gdouble worst_end = 0;

  for (k = 0; k < G_N_ELEMENTS (stretch_cases); k++)
    {
      const StretchCase *c = &stretch_cases[k];
      gfloat            *y = stretch_ramp (c, x, c->inverse);
      gdouble            y0, y1, e;
      gboolean           mono = TRUE;

      if (f)
        {
          if (c->arcsinh)
            fprintf (f, "case %u arcsinh %.17g %.17g\n", k, c->beta, c->bp);
          else
            fprintf (f, "case %u ghs %.17g %.17g %.17g %.17g %.17g %.17g %d\n", k, c->b,
                     c->stretch, c->lp, c->sp, c->hp, c->bp, c->inverse);
          for (i = 0; i < RAMP_N; i++)
            fprintf (f, "%.9g %.9g\n", x[i], y[i]);
        }
      /* the even ramp: rising, 0 (or the black point) to 0, 1 to 1 */
      for (i = 1; i < 4096; i++)
        if (y[i] < y[i - 1])
          mono = FALSE;
      y0 = y[0];
      y1 = y[4095];
      e = MAX (fabs (y0), fabs (y1 - 1));
      worst_end = MAX (worst_end, e);
      bad_mono += ! mono;
      bad_ends += e > 1e-6;
      g_free (y);
    }
  if (f)
    fclose (f);
  report ("stretch_written_for_reference", f != NULL, "%s", path);
  report ("stretch_monotonic", bad_mono == 0, "%u cases, %d not rising",
          (guint) G_N_ELEMENTS (stretch_cases), bad_mono);
  report ("stretch_0_and_1_kept", bad_ends == 0,
          "%u cases, largest error at 0 or 1: %.2g", (guint) G_N_ELEMENTS (stretch_cases),
          worst_end);
  g_free (x);
  g_free (path);
}

/* the inverse undoes the stretch */
static void
check_stretch_inverse (void)
{
  gfloat *x = g_new (gfloat, RAMP_N);
  gdouble worst = 0;
  guint   k;
  gint    i, n = 0;

  for (k = 1; k < G_N_ELEMENTS (stretch_cases); k++)
    {
      const StretchCase *c = &stretch_cases[k];
      gfloat            *y, *p, *back;

      if (c->arcsinh || c->inverse || c->bp > 0)
        continue;
      y = stretch_ramp (c, x, FALSE);
      p = g_new (gfloat, RAMP_N * 4);
      for (i = 0; i < RAMP_N; i++)
        {
          p[i * 4 + 0] = p[i * 4 + 1] = p[i * 4 + 2] = y[i];
          p[i * 4 + 3] = 1.0f;
        }
      back = run_size (p, RAMP_N, 1, "astro:stretch", "colour", 1, "stretch", c->stretch,
                       "local-intensity", c->b, "symmetry-point", c->sp,
                       "protect-shadows", c->lp, "protect-highlights", c->hp,
                       "inverse", TRUE, NULL);
      /* only where the stretch did not flatten the curve below float
       * resolution (the slope of the stretch above 1/64) */
      for (i = 1; i < 4095; i++)
        if ((y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1]) > 1.0 / 64)
          worst = MAX (worst, fabs (back[i * 4] - x[i]));
      n++;
      g_free (back);
      g_free (p);
      g_free (y);
    }
  report ("stretch_inverse_undoes", worst < 2e-5,
          "%d GHS cases forward then inverse: largest error %.2g", n, worst);
  g_free (x);
}

/* preserve colour: R:G:B kept, the luminance follows the curve, and a
 * pixel that would leave 0 to 1 is scaled down whole */
static void
check_stretch_colour (void)
{
  const gint n = 4000;
  gfloat    *p = g_new (gfloat, n * 4), *grey = g_new (gfloat, n * 4), *o, *g;
  GRand     *rnd = g_rand_new_with_seed (11);
  gdouble    lum[3], ratio = 0, ylum = 0, over = 0;
  gint       i, c, scaled = 0;

  babl_space_get_rgb_luminance (babl_space ("sRGB"), &lum[0], &lum[1], &lum[2]);
  for (i = 0; i < n; i++)
    {
      gdouble Y = 0;

      for (c = 0; c < 3; c++)
        {
          /* most pixels dark, some bright and coloured */
          p[i * 4 + c] = i < n / 2 ? g_rand_double_range (rnd, 0.0, 0.1)
                                   : g_rand_double_range (rnd, 0.0, 0.9);
          Y += lum[c] * p[i * 4 + c];
        }
      p[i * 4 + 3] = 1.0f;
      grey[i * 4 + 0] = grey[i * 4 + 1] = grey[i * 4 + 2] = Y;
      grey[i * 4 + 3] = 1.0f;
    }
  o = run_size (p, n, 1, "astro:stretch", "stretch", 5.0, "local-intensity", 2.0,
                "symmetry-point", 0.05, NULL);
  g = run_size (grey, n, 1, "astro:stretch", "stretch", 5.0, "local-intensity", 2.0,
                "symmetry-point", 0.05, "colour", 1, NULL);
  for (i = 0; i < n; i++)
    {
      gdouble Yo = 0, m = 0, f = 0;

      for (c = 0; c < 3; c++)
        {
          Yo += lum[c] * o[i * 4 + c];
          m = MAX (m, o[i * 4 + c]);
          if (p[i * 4 + c] > 0.02f)
            f = o[i * 4 + c] / p[i * 4 + c];
        }
      /* the ratios: each channel times the same factor */
      for (c = 0; c < 3; c++)
        if (p[i * 4 + c] > 0.02f)
          ratio = MAX (ratio, fabs (o[i * 4 + c] / p[i * 4 + c] / f - 1));
      over = MAX (over, m - 1);
      if (g[i * 4] * MAX (p[i * 4], MAX (p[i * 4 + 1], p[i * 4 + 2])) /
          grey[i * 4] > 1.0 + 1e-6)
        scaled++;
      else
        ylum = MAX (ylum, fabs (Yo - g[i * 4]));
    }
  report ("stretch_preserve_ratios", ratio < 1e-5,
          "%d pixels: channel ratios kept to %.2g (relative)", n, ratio);
  report ("stretch_preserve_luminance", ylum < 2e-6,
          "luminance follows the curve to %.2g (%d pixels scaled down to stay in range)",
          ylum, scaled);
  report ("stretch_preserve_in_range", over <= 1e-6 && scaled > 0,
          "largest channel %.7f above 1, %d pixels scaled", over, scaled);
  g_free (o);
  g_free (g);
  g_free (p);
  g_free (grey);
  g_rand_free (rnd);
}

static int
cmp_gfloat (const void *a,
            const void *b)
{
  gfloat x = *(const gfloat *) a, y = *(const gfloat *) b;

  return (x > y) - (x < y);
}

static gdouble
median_of (const gfloat *img,
           gint          c)
{
  gfloat *v = g_new (gfloat, W * H);
  gdouble m;
  gint    i;

  for (i = 0; i < W * H; i++)
    v[i] = img[i * 4 + c];
  qsort (v, W * H, sizeof (gfloat), cmp_gfloat);
  m = v[W * H / 2];
  g_free (v);
  return m;
}

/* the automatic black and symmetry point on a sky of known level and
 * noise, per channel: like the same settings by hand; whole-image, so the same bits
 * wherever the layer sits */
static void
check_stretch_auto (void)
{
  Field   f = make_field (21, 300, sky_vignette, FALSE, 0.0, 0.004);
  gfloat *a, *b, *m;
  gdouble ma, mm, zeros = 0;
  gint    i;

  /* a flat sky of 0.08 instead of the vignetting: the known level */
  for (i = 0; i < W * H; i++)
    {
      gint c;

      for (c = 0; c < 3; c++)
        f.image[i * 4 + c] = 0.08f + f.stars[i * 3 + c] + f.noise[i * 3 + c];
    }
  a = run (f.image, 0, 0, "astro:stretch", "stretch", 6.0, "local-intensity", 3.0,
           "auto-points", TRUE, "colour", 1, NULL);
  b = run (f.image, 77, -33, "astro:stretch", "stretch", 6.0, "local-intensity", 3.0,
           "auto-points", TRUE, "colour", 1, NULL);
  /* by hand: the median and robust sigma of the whole green channel,
   * black point 3 sigma below the median, symmetry point at the median */
  {
    gfloat *v = g_new (gfloat, W * H);
    gdouble med, sig, bp;

    for (i = 0; i < W * H; i++)
      v[i] = f.image[i * 4 + 1];
    qsort (v, W * H, sizeof (gfloat), cmp_gfloat);
    med = v[W * H / 2];
    for (i = 0; i < W * H; i++)
      v[i] = fabsf (v[i] - (gfloat) med);
    qsort (v, W * H, sizeof (gfloat), cmp_gfloat);
    sig = 1.4826 * v[W * H / 2];
    g_free (v);
    bp = med - 3 * sig;
    m = run (f.image, 0, 0, "astro:stretch", "stretch", 6.0, "local-intensity", 3.0,
             "black-point", bp, "symmetry-point", (med - bp) / (1 - bp), "colour", 1, NULL);
    printf ("NOTE  stretch_auto: whole image green median %.5f, robust sigma %.5f "
            "(made with 0.08 and 0.004)\n", med, sig);
  }
  ma = median_of (a, 1);
  mm = median_of (m, 1);
  for (i = 0; i < W * H; i++)
    zeros += a[i * 4 + 1] == 0.0f;
  zeros /= W * H;
  report ("stretch_auto_like_by_hand", fabs (ma / mm - 1) < 0.02,
          "sky median after the stretch %.4f, with the points by hand %.4f", ma, mm);
  /* 3 sigma below the median: 0.13 % of the noise is under the black point */
  report ("stretch_auto_black_point", zeros > 0.0005 && zeros < 0.003,
          "%.3f %% of the pixels at 0 (Gaussian noise: 0.135 %%)", 100 * zeros);
  report ("stretch_auto_offset_invariant", memcmp (a, b, sizeof (gfloat) * W * H * 4) == 0,
          "buffer at 77, -33 %s", memcmp (a, b, sizeof (gfloat) * W * H * 4) ? "differs" : "identical");
  g_free (a);
  g_free (b);
  g_free (m);
  free_field (&f);
}

/* the default settings change nothing (stretch 0), alpha is kept */
static void
check_stretch_default (void)
{
  Field   f = make_field (22, 300, sky_quadratic, FALSE, 0.0, 0.003);
  gfloat *o;
  gdouble d = 0;
  gint    i;

  for (i = 0; i < W * H; i++)
    f.image[i * 4 + 3] = (i % 5) / 4.0f;
  o = run (f.image, 0, 0, "astro:stretch", NULL);
  for (i = 0; i < W * H * 4; i++)
    d = MAX (d, fabs (o[i] - ((i % 4) == 3 ? f.image[i] : CLAMP (f.image[i], 0.0f, 1.0f))));
  report ("stretch_default_unchanged", d < 1e-6,
          "defaults (stretch 0): max difference %.2g, alpha kept", d);
  g_free (o);
  free_field (&f);
}

int
main (int    argc,
      char **argv)
{
  gchar *dir;

  if (argc < 2)
    {
      fprintf (stderr, "usage: %s <build folder>\n", argv[0]);
      return 2;
    }
  setlocale (LC_ALL, "C");
  build_dir = argv[1];
  dir = modules_setup (argv[1]);
  if (! dir)
    return 2;
  gegl_init (&argc, &argv);
  g_log_set_default_handler (log_handler, NULL);

  check_registration ();
  check_flat ();
  check_gradient_recovery ();
  check_degrees ();
  check_nebula ();
  check_vignetting ();
  check_determinism ();
  check_alpha ();
  check_stretch_default ();
  check_stretch_curves ();
  check_stretch_inverse ();
  check_stretch_colour ();
  check_stretch_auto ();

  report ("no_warnings_or_criticals", n_messages == 0, "%d", n_messages);
  gegl_exit ();
  if (! g_getenv ("ASTRO_CHECK_KEEP_MODULE"))
    modules_cleanup (dir);
  else
    g_free (dir);
  printf ("%d passed, %d failed\n", n_passed, n_failed);
  return n_failed ? 1 : 0;
}
