/*
 * How long the operations take on a 24 megapixel star field
 *
 * bench.c
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
 *   bench <build folder> [width height]
 *
 * Times each operation on a 6000 x 4000 synthetic star field (a sky
 * gradient, 30000 Gaussian stars, noise; fixed seed) in float, with
 * GEGL's threads: the best of three runs after a first one (which
 * allocates the output). GEGL's gegl:gaussian-blur (std-dev 40) is timed
 * the same way, for scale. (After gegl-presence's tests/bench.c.)
 */

#include <gegl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "modules.h"

static gint W = 6000, H = 4000;

static GeglBuffer *
make_field (void)
{
  GeglBuffer *b   = gegl_buffer_new (GEGL_RECTANGLE (0, 0, W, H), babl_format ("RGBA float"));
  gfloat     *img = g_new (gfloat, (gsize) W * H * 4);
  GRand      *rnd = g_rand_new_with_seed (7);
  gsize       i;
  gint        s, x, y, c;

  for (y = 0; y < H; y++)
    for (x = 0; x < W; x++)
      {
        gdouble u = 2.0 * x / W - 1, v = 2.0 * y / H - 1;
        gdouble sky = 0.05 + 0.02 * u + 0.015 * v + 0.01 * u * u;
        gfloat *p = img + ((gsize) y * W + x) * 4;

        p[0] = sky * 1.2;
        p[1] = sky;
        p[2] = sky * 0.9;
        p[3] = 1.0f;
      }
  for (s = 0; s < 30000; s++)
    {
      gdouble sx = g_rand_double_range (rnd, 0, W), sy = g_rand_double_range (rnd, 0, H);
      gdouble sig = g_rand_double_range (rnd, 0.8, 2.5);
      gdouble peak = MIN (0.9, 0.02 * pow (g_rand_double_range (rnd, 1e-3, 1.0), -0.7));
      gint    r = (gint) ceil (4 * sig);

      for (y = MAX (0, (gint) sy - r); y <= MIN (H - 1, (gint) sy + r); y++)
        for (x = MAX (0, (gint) sx - r); x <= MIN (W - 1, (gint) sx + r); x++)
          {
            gdouble dx = x + 0.5 - sx, dy = y + 0.5 - sy;
            gdouble g = peak * exp (-(dx * dx + dy * dy) / (2 * sig * sig));

            for (c = 0; c < 3; c++)
              img[((gsize) y * W + x) * 4 + c] += g;
          }
    }
  /* noise: a cheap uniform sum (Irwin-Hall of 4), sigma 0.003 */
  for (i = 0; i < (gsize) W * H * 4; i++)
    if (i % 4 != 3)
      img[i] += 0.003f * 1.732f * (gfloat) (g_rand_double (rnd) + g_rand_double (rnd) +
                                            g_rand_double (rnd) + g_rand_double (rnd) - 2.0);
  gegl_buffer_set (b, GEGL_RECTANGLE (0, 0, W, H), 0, babl_format ("RGBA float"), img,
                   GEGL_AUTO_ROWSTRIDE);
  g_free (img);
  g_rand_free (rnd);
  return b;
}

/* seconds for the operation over the whole image, the best of three;
 * props is a space separated list of name=value */
static gdouble
time_op (GeglBuffer  *in,
         const gchar *op,
         const gchar *props)
{
  GeglBuffer *out  = gegl_buffer_new (GEGL_RECTANGLE (0, 0, W, H), gegl_buffer_get_format (in));
  gchar     **p    = g_strsplit (props ? props : "", " ", -1);
  gdouble     best = G_MAXDOUBLE;
  gint        run, i;

  for (run = 0; run < 4; run++)
    {
      GeglNode *g    = gegl_node_new ();
      GeglNode *src  = gegl_node_new_child (g, "operation", "gegl:buffer-source",
                                            "buffer", in, NULL);
      GeglNode *node = gegl_node_new_child (g, "operation", op, NULL);
      GeglNode *sink = gegl_node_new_child (g, "operation", "gegl:write-buffer",
                                            "buffer", out, NULL);
      gint64    t0;

      for (i = 0; p[i]; i++)
        if (p[i][0] && ! node_set_from_string (node, p[i]))
          fprintf (stderr, "cannot set %s\n", p[i]);
      gegl_node_link_many (src, node, sink, NULL);
      t0 = g_get_monotonic_time ();
      gegl_node_process (sink);
      if (run > 0)
        best = MIN (best, (g_get_monotonic_time () - t0) / 1e6);
      g_object_unref (g);
    }
  g_strfreev (p);
  g_object_unref (out);
  return best;
}

int
main (int    argc,
      char **argv)
{
  static const gchar *const rows[][3] =
    {
      { "gegl:gaussian-blur", "std-dev-x=40 std-dev-y=40", "gegl:gaussian-blur 40" },
      { "astro:background", "", "Remove Sky Background, defaults" },
      { "astro:background", "samples=32 degree=4", "Remove Sky Background, 32 cells, degree 4" },
      { "astro:background", "exclude-stars=true", "Remove Sky Background, stars left out" },
      { "astro:stretch", "stretch=4 local-intensity=2 symmetry-point=0.05", "Astro Stretch, GHS" },
      { "astro:stretch", "curve=arcsinh beta=50", "Astro Stretch, arcsinh" },
      { "astro:stretch", "stretch=4 local-intensity=2 auto-points=true", "Astro Stretch, GHS, automatic points" },
      { "astro:star-mask", "", "Star Mask, defaults" },
      { "astro:star-reduce", "", "Reduce Stars, defaults" },
      { "astro:star-reduce", "method=erosion amount=100 iterations=2", "Reduce Stars, erosion, 2 iterations" },
    };
  gchar      *dir;
  GeglBuffer *img;
  gint        threads;
  guint       r;

  if (argc != 2 && argc != 4)
    {
      fprintf (stderr, "usage: %s <build folder> [width height]\n", argv[0]);
      return 2;
    }
  if (argc == 4)
    {
      W = atoi (argv[2]);
      H = atoi (argv[3]);
    }
  dir = modules_setup (argv[1]);
  if (! dir)
    return 2;
  gegl_init (NULL, NULL);
  /* room for the images in memory, as GIMP gives GEGL; swapping to disk
   * would be what is timed otherwise */
  g_object_set (gegl_config (), "tile-cache-size", (guint64) 8 << 30, NULL);
  g_object_get (gegl_config (), "threads", &threads, NULL);
  img = make_field ();
  printf ("%d x %d pixels (%.1f megapixels), float, GEGL threads: %d\n\n",
          W, H, W * H / 1e6, threads);
  for (r = 0; r < G_N_ELEMENTS (rows); r++)
    {
      if (! gegl_has_operation (rows[r][0]))
        continue;
      printf ("  %-46s %6.3f s\n", rows[r][2], time_op (img, rows[r][0], rows[r][1]));
      fflush (stdout);
    }
  g_object_unref (img);
  gegl_exit ();
  modules_cleanup (dir);
  return 0;
}
