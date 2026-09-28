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

static const gchar *const modules[] = { "astro-background" };

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

/* runs op on the pixels (W x H RGBA float at x0, y0) with the properties
 * as name, value pairs (NULL ended), into a new array */
static gfloat *
run (const gfloat *pixels,
     gint          x0,
     gint          y0,
     const gchar  *op,
     ...)
{
  GeglRectangle r = { x0, y0, W, H };
  const Babl   *fmt = babl_format (FMT);
  GeglBuffer   *in = gegl_buffer_new (&r, fmt);
  GeglBuffer   *out = NULL;
  GeglNode     *g = gegl_node_new ();
  GeglNode     *src, *node, *sink;
  gfloat       *res = g_new (gfloat, W * H * 4);
  va_list       ap;

  gegl_buffer_set (in, &r, 0, fmt, pixels, GEGL_AUTO_ROWSTRIDE);
  src  = gegl_node_new_child (g, "operation", "gegl:buffer-source", "buffer", in, NULL);
  node = gegl_node_new_child (g, "operation", op, NULL);
  va_start (ap, op);
  gegl_node_set_valist (node, va_arg (ap, const gchar *), ap);
  va_end (ap);
  sink = gegl_node_new_child (g, "operation", "gegl:buffer-sink", "buffer", &out, NULL);
  gegl_node_link_many (src, node, sink, NULL);
  gegl_node_process (sink);
  gegl_buffer_get (out, &r, 1.0, fmt, res, GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);
  g_object_unref (out);
  g_object_unref (g);
  g_object_unref (in);
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
  gboolean ok = gegl_has_operation ("astro:background");
  const gchar *path = ok ? gegl_operation_get_key ("astro:background", "gimp:menu-path") : NULL;

  report ("background_registered", ok && g_strcmp0 (path, "<Image>/Filters/Astronomy") == 0,
          "operation %s, menu %s", ok ? "found" : "missing", path ? path : "none");
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
  gfloat *model = run (f.image, 0, 0, "astro:background", "output", 1, NULL);
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
      gfloat *model = run (f.image, 0, 0, "astro:background", "output", 1, "degree", d, NULL);
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
  gfloat *model = run (f.image, 0, 0, "astro:background", "output", 1, NULL);
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
  model = run (f.image, 0, 0, "astro:background", "output", 1, NULL);
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

  report ("no_warnings_or_criticals", n_messages == 0, "%d", n_messages);
  gegl_exit ();
  if (! g_getenv ("ASTRO_CHECK_KEEP_MODULE"))
    modules_cleanup (dir);
  else
    g_free (dir);
  printf ("%d passed, %d failed\n", n_passed, n_failed);
  return n_failed ? 1 : 0;
}
