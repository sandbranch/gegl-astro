/*
 * The checks of astro:star-mask and astro:star-reduce: star fields with
 * a known list of stars (position, size, peak), included by check.c
 * (it shares its helpers: report, run, run_size, make_field).
 *
 * check-stars.c
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

#define MAX_STARS 1024
#define N_SNR     6

/* peak over noise sigma of the star classes */
static const gdouble snr_class[N_SNR] = { 3, 5, 8, 12, 20, 50 };

typedef struct
{
  gdouble x, y, sigma, peak;
  gint    cls;        /* index into snr_class, or -1 for bright stars */
} Star;

typedef struct
{
  gint    n;
  Star    star[MAX_STARS];
  gdouble sky, noise;
  gfloat *image;      /* RGBA */
  gfloat *stars;      /* the stars alone, luminance (grey stars) */
  gfloat *nebula;     /* the nebula alone */
} StarField;

/* isolated grey stars on a jittered grid of spacing px: the classes of
 * snr_class in turn (and every seventh a bright star, peak 0.5 to 0.9),
 * sigma 0.9 to 2.2 px; a flat sky; Gaussian noise; an optional nebula
 * (a wide Gaussian, sigma 60 px, and a compact knot, sigma 8 px: FWHM
 * 18.8 px, above the default largest star of 15) */
static StarField *
star_field (guint32 seed,
            gint    spacing,
            gdouble noise,
            gdouble nebula,
            gdouble colour)
{
  StarField *f = g_new0 (StarField, 1);
  GRand     *rnd = g_rand_new_with_seed (seed);
  gint       gx, gy, i, x, y, c;

  f->sky    = 0.05;
  f->noise  = noise;
  f->image  = g_new0 (gfloat, W * H * 4);
  f->stars  = g_new0 (gfloat, W * H);
  f->nebula = g_new0 (gfloat, W * H);
  for (gy = spacing / 2; gy < H - spacing / 2 + 1 && f->n < MAX_STARS; gy += spacing)
    for (gx = spacing / 2; gx < W - spacing / 2 + 1 && f->n < MAX_STARS; gx += spacing)
      {
        Star *s = &f->star[f->n];

        s->x     = gx + g_rand_double_range (rnd, -0.2, 0.2) * spacing;
        s->y     = gy + g_rand_double_range (rnd, -0.2, 0.2) * spacing;
        s->sigma = g_rand_double_range (rnd, 0.9, 2.2);
        if (f->n % 7 == 6)
          {
            s->cls  = -1;
            s->peak = g_rand_double_range (rnd, 0.5, 0.9);
          }
        else
          {
            s->cls  = f->n % N_SNR;
            s->peak = snr_class[s->cls] * MAX (noise, 1e-4);
          }
        f->n++;
      }
  for (i = 0; i < f->n; i++)
    {
      const Star *s = &f->star[i];
      gint        r = (gint) ceil (5 * s->sigma);

      for (y = MAX (0, (gint) s->y - r); y <= MIN (H - 1, (gint) s->y + r); y++)
        for (x = MAX (0, (gint) s->x - r); x <= MIN (W - 1, (gint) s->x + r); x++)
          {
            gdouble dx = x + 0.5 - s->x, dy = y + 0.5 - s->y;

            f->stars[y * W + x] += s->peak * exp (-(dx * dx + dy * dy) / (2 * s->sigma * s->sigma));
          }
    }
  if (nebula > 0)
    for (y = 0; y < H; y++)
      for (x = 0; x < W; x++)
        {
          gdouble dx = x - 0.4 * W, dy = y - 0.5 * H;
          gdouble kx = x - 0.75 * W, ky = y - 0.3 * H;

          f->nebula[y * W + x] = nebula * exp (-(dx * dx + dy * dy) / (2 * 60.0 * 60.0)) +
                                 0.5 * nebula * exp (-(kx * kx + ky * ky) / (2 * 8.0 * 8.0));
        }
  for (i = 0; i < W * H; i++)
    {
      for (c = 0; c < 3; c++)
        {
          /* colour > 0: the stars reddish (R up, B down by that share) */
          gdouble tint = c == 0 ? 1 + colour : c == 2 ? 1 - colour : 1;

          f->image[i * 4 + c] = f->sky + f->nebula[i] + f->stars[i] * tint +
                                noise * gauss (rnd);
        }
      f->image[i * 4 + 3] = 1.0f;
    }
  g_rand_free (rnd);
  return f;
}

static void
star_field_free (StarField *f)
{
  g_free (f->image);
  g_free (f->stars);
  g_free (f->nebula);
  g_free (f);
}

static gfloat
at (const gfloat *rgba,
    gdouble       x,
    gdouble       y)
{
  gint ix = CLAMP ((gint) x, 0, W - 1), iy = CLAMP ((gint) y, 0, H - 1);

  return rgba[(iy * W + ix) * 4];
}

/* connected regions of mask >= 0.5 that hold no star's centre */
static gint
false_regions (const gfloat    *mask,
               const StarField *f)
{
  gint   *label = g_new0 (gint, W * H);
  gint   *stack = g_new (gint, W * H);
  gint    i, n_regions = 0, n_false = 0;
  gboolean *has;

  for (i = 0; i < W * H; i++)
    if (mask[i * 4] >= 0.5f && ! label[i])
      {
        gint top = 0;

        n_regions++;
        label[i] = n_regions;
        stack[top++] = i;
        while (top)
          {
            gint p = stack[--top], px = p % W, py = p / W, dx, dy;

            for (dy = -1; dy <= 1; dy++)
              for (dx = -1; dx <= 1; dx++)
                {
                  gint xx = px + dx, yy = py + dy, q;

                  if (xx < 0 || yy < 0 || xx >= W || yy >= H)
                    continue;
                  q = yy * W + xx;
                  if (mask[q * 4] >= 0.5f && ! label[q])
                    {
                      label[q] = n_regions;
                      stack[top++] = q;
                    }
                }
          }
      }
  has = g_new0 (gboolean, n_regions + 1);
  for (i = 0; f && i < f->n; i++)
    {
      gint ix = CLAMP ((gint) f->star[i].x, 0, W - 1), iy = CLAMP ((gint) f->star[i].y, 0, H - 1);

      has[label[iy * W + ix]] = TRUE;
    }
  for (i = 1; i <= n_regions; i++)
    n_false += ! has[i];

  g_free (has);
  g_free (stack);
  g_free (label);
  return n_false;
}

/* ---- astro:star-mask ---------------------------------------------------- */

static void
check_star_mask_completeness (void)
{
  StarField *f = star_field (31, 20, 0.004, 0.0, 0.0);
  gfloat    *m = run (f->image, 0, 0, "astro:star-mask", NULL);
  gint       found[N_SNR] = { 0 }, total[N_SNR] = { 0 }, i, k, nf, bright = 0, nbright = 0;
  gdouble    share_sum = 0;
  gint       share_n = 0;
  GString   *s = g_string_new (NULL);

  for (i = 0; i < f->n; i++)
    {
      const Star *st = &f->star[i];
      gboolean    hit = at (m, st->x, st->y) >= 0.5f;

      if (st->cls < 0)
        {
          nbright++;
          bright += hit;
        }
      else
        {
          total[st->cls]++;
          found[st->cls] += hit;
        }
      /* the share of the star's light under the mask (SNR 20 and up) */
      if (st->cls < 0 || snr_class[st->cls] >= 20)
        {
          gint    r = (gint) ceil (4 * st->sigma), x, y;
          gdouble in = 0, all = 0;

          for (y = MAX (0, (gint) st->y - r); y <= MIN (H - 1, (gint) st->y + r); y++)
            for (x = MAX (0, (gint) st->x - r); x <= MIN (W - 1, (gint) st->x + r); x++)
              {
                gdouble dx = x + 0.5 - st->x, dy = y + 0.5 - st->y;
                gdouble g = st->peak * exp (-(dx * dx + dy * dy) / (2 * st->sigma * st->sigma));

                all += g;
                in  += g * m[(y * W + x) * 4];
              }
          share_sum += in / all;
          share_n++;
        }
    }
  for (k = 0; k < N_SNR; k++)
    g_string_append_printf (s, "%sSNR %g: %d/%d", k ? ", " : "", snr_class[k], found[k], total[k]);
  printf ("NOTE  star_mask_completeness: threshold 4, %s, bright %d/%d\n", s->str, bright, nbright);
  report ("star_mask_complete_snr_12_up",
          found[3] + found[4] + found[5] + bright >= 0.97 * (total[3] + total[4] + total[5] + nbright),
          "%d of %d stars with peak at least 12 noise sigmas found (limit 97 %%)",
          found[3] + found[4] + found[5] + bright, total[3] + total[4] + total[5] + nbright);
  report ("star_mask_complete_snr_8", found[2] >= 0.8 * total[2],
          "%d of %d at 8 sigmas (limit 80 %%)", found[2], total[2]);
  report ("star_mask_faint_noise_left_out", found[0] <= 0.5 * total[0],
          "%d of %d at 3 sigmas marked (below the threshold of 4; limit 50 %%)",
          found[0], total[0]);
  nf = false_regions (m, f);
  report ("star_mask_false_detections", nf <= 3,
          "%d regions without a star in %d stars (limit 3)", nf, f->n);
  report ("star_mask_flux_inside", share_sum / share_n >= 0.75,
          "mean share of a star's light under the mask %.3f (%d stars at 20 sigmas and up; "
          "limit 0.75)", share_sum / share_n, share_n);
  g_string_free (s, TRUE);
  g_free (m);
  star_field_free (f);
}

/* pure noise (and a gradient): no stars; four fields, about a megapixel */
static void
check_star_mask_noise_only (void)
{
  gint t, nf = 0, n3 = 0;

  for (t = 0; t < 4; t++)
    {
      Field   f = make_field (32 + t, 0, sky_quadratic, FALSE, 0.0, 0.004);
      gfloat *m = run (f.image, 0, 0, "astro:star-mask", NULL);
      gfloat *m3 = run (f.image, 0, 0, "astro:star-mask", "threshold", 3.0, NULL);

      nf += false_regions (m, NULL);
      n3 += false_regions (m3, NULL);
      g_free (m3);
      g_free (m);
      free_field (&f);
    }
  report ("star_mask_noise_only", nf / (4 * W * H / 1e6) <= 8,
          "noise of 0.004 on a gradient, 4 fields of 0.24 megapixels: %d false regions, "
          "%.1f per megapixel (limit 8; at threshold 3: %.0f per megapixel)", nf,
          nf / (4 * W * H / 1e6), n3 / (4 * W * H / 1e6));
}

/* the blur of astro-common.h against a plain sum with mirrored edges,
 * also with kernels longer than the plane */
static void
check_blur (void)
{
  static const gint sizes[][2] = { { 37, 23 }, { 300, 5 }, { 4, 90 }, { 1, 1 } };
  static const gdouble sigmas[] = { 0.7, 1.5, 6.0, 25.0 };
  GRand  *rnd = g_rand_new_with_seed (51);
  gdouble worst = 0;
  guint   a, b;

  for (a = 0; a < G_N_ELEMENTS (sizes); a++)
    for (b = 0; b < G_N_ELEMENTS (sigmas); b++)
      {
        gint     w = sizes[a][0], h = sizes[a][1], x, y, i, j, r;
        gfloat  *p = g_new (gfloat, w * h), *k = astro_gauss_kernel (sigmas[b], &r);
        gdouble *want = g_new (gdouble, w * h);

        for (i = 0; i < w * h; i++)
          p[i] = g_rand_double (rnd);
        for (y = 0; y < h; y++)
          for (x = 0; x < w; x++)
            {
              gdouble sum = 0;

              for (j = -r; j <= r; j++)
                for (i = -r; i <= r; i++)
                  sum += (gdouble) k[j + r] * k[i + r] *
                         p[astro_mirror (y + j, h) * w + astro_mirror (x + i, w)];
              want[y * w + x] = sum;
            }
        astro_gauss_blur (p, w, h, sigmas[b]);
        for (i = 0; i < w * h; i++)
          worst = MAX (worst, fabs (p[i] - want[i]));
        g_free (p);
        g_free (k);
        g_free (want);
      }
  report ("common_blur_like_plain_sum", worst < 1e-5,
          "16 planes and sigmas: largest difference %.2g", worst);
  g_rand_free (rnd);
}

/* a nebula and a compact knot: not in the mask; the stars on the nebula
 * found */
static void
check_star_mask_nebula (void)
{
  StarField *f = star_field (33, 24, 0.003, 0.15, 0.0);
  gfloat    *m = run (f->image, 0, 0, "astro:star-mask", NULL);
  gdouble    neb_mask = 0, neb_n = 0;
  gint       i, on = 0, on_found = 0;

  for (i = 0; i < W * H; i++)
    if (f->nebula[i] > 0.03 && f->stars[i] < 0.1 * f->noise)
      {
        neb_mask += m[i * 4];
        neb_n++;
      }
  for (i = 0; i < f->n; i++)
    {
      const Star *st = &f->star[i];

      if (f->nebula[(gint) st->y * W + (gint) st->x] > 0.05 &&
          (st->cls < 0 || snr_class[st->cls] >= 12))
        {
          on++;
          on_found += at (m, st->x, st->y) >= 0.5f;
        }
    }
  report ("star_mask_nebula_left_out", neb_mask / neb_n < 0.02,
          "mean mask on the nebula away from stars %.4f (limit 0.02)", neb_mask / neb_n);
  report ("star_mask_knot_left_out", at (m, 0.75 * W, 0.3 * H) < 0.5f,
          "mask at the knot's centre (sigma 8 px, FWHM 18.8) %.3f", at (m, 0.75 * W, 0.3 * H));
  report ("star_mask_stars_on_nebula", on_found >= 0.9 * on,
          "%d of %d stars (12 sigmas and up) on the nebula found", on_found, on);
  g_free (m);
  star_field_free (f);
}

static void
check_star_mask_misc (void)
{
  StarField *f = star_field (34, 22, 0.004, 0.05, 0.0);
  gfloat    *a = run (f->image, 0, 0, "astro:star-mask", NULL);
  gfloat    *b = run (f->image, 0, 0, "astro:star-mask", NULL);
  gfloat    *c = run (f->image, -41, 93, "astro:star-mask", NULL);
  gfloat    *o = run (f->image, 0, 0, "astro:star-mask", "show", 1, NULL);
  gdouble    d = 0, junk = 0;
  gint       i, x, y;

  report ("star_mask_deterministic", memcmp (a, b, sizeof (gfloat) * W * H * 4) == 0,
          "two runs %s", memcmp (a, b, sizeof (gfloat) * W * H * 4) ? "differ" : "identical");
  report ("star_mask_offset_invariant", memcmp (a, c, sizeof (gfloat) * W * H * 4) == 0,
          "buffer at -41, 93 %s", memcmp (a, c, sizeof (gfloat) * W * H * 4) ? "differs" : "identical");
  for (i = 0; i < W * H; i++)
    if (a[i * 4] == 0.0f)
      for (x = 0; x < 4; x++)
        d = MAX (d, fabs (o[i * 4 + x] - f->image[i * 4 + x]));
  report ("star_mask_overlay", d == 0, "outside the mask the overlay is the image: %.2g", d);

  /* a transparent block of junk: no stars there */
  for (y = 100; y < 250; y++)
    for (x = 300; x < 500; x++)
      {
        gfloat *p = f->image + (y * W + x) * 4;

        p[0] = p[1] = p[2] = (x + y) % 5 ? 0.0f : 3.0f;
        p[3] = 0.0f;
      }
  g_free (b);
  b = run (f->image, 0, 0, "astro:star-mask", NULL);
  for (y = 110; y < 240; y++)
    for (x = 310; x < 490; x++)
      junk = MAX (junk, b[(y * W + x) * 4]);
  report ("star_mask_transparent_ignored", junk == 0, "largest mask in the block %.2g", junk);
  g_free (a);
  g_free (b);
  g_free (c);
  g_free (o);
  star_field_free (f);
}

/* ---- astro:star-reduce -------------------------------------------------- */

/* the FWHM of a star: the area, in a disc of 4 sigma of the true star,
 * where green minus the sky is at least half its peak, taken as a disc.
 * (Second moments over such a disc were tried first, as PLAN.md said:
 * they weight the edge by r^2, so a sky offset of a few thousandths
 * there decides the result, e.g. the darkened noise of an erosion made
 * stars look 25 % smaller. The half maximum does not see the sky.) */
static gdouble
star_fwhm (const gfloat *rgba,
           const Star   *st,
           gdouble       sky)
{
  gint    r = (gint) ceil (4 * st->sigma), x, y, half = 0;
  gdouble peak = 0;

  for (y = MAX (0, (gint) st->y - r); y <= MIN (H - 1, (gint) st->y + r); y++)
    for (x = MAX (0, (gint) st->x - r); x <= MIN (W - 1, (gint) st->x + r); x++)
      peak = MAX (peak, rgba[(y * W + x) * 4 + 1] - sky);
  for (y = MAX (0, (gint) st->y - r); y <= MIN (H - 1, (gint) st->y + r); y++)
    for (x = MAX (0, (gint) st->x - r); x <= MIN (W - 1, (gint) st->x + r); x++)
      {
        gdouble dx = x + 0.5 - st->x, dy = y + 0.5 - st->y;

        if (dx * dx + dy * dy <= r * r && rgba[(y * W + x) * 4 + 1] - sky >= 0.5 * peak)
          half++;
      }
  return 2.0 * sqrt (half / G_PI);
}

static int
cmp_double_asc (const void *a,
                const void *b)
{
  gdouble x = *(const gdouble *) a, y = *(const gdouble *) b;

  return (x > y) - (x < y);
}

/* the median FWHM ratio (after / before) of the bright stars and of the
 * stars at 50 sigmas */
static gdouble
fwhm_ratio (const StarField *f,
            const gfloat    *out,
            gdouble         *flux_ratio)
{
  gdouble r[MAX_STARS], fl[MAX_STARS], med;
  gint    i, n = 0;

  for (i = 0; i < f->n; i++)
    {
      const Star *st = &f->star[i];
      gint        rr = (gint) ceil (4 * st->sigma), x, y;
      gdouble     a = 0, b = 0;

      if (st->cls >= 0 && snr_class[st->cls] < 50)
        continue;
      r[n] = star_fwhm (out, st, f->sky) / star_fwhm (f->image, st, f->sky);
      for (y = MAX (0, (gint) st->y - rr); y <= MIN (H - 1, (gint) st->y + rr); y++)
        for (x = MAX (0, (gint) st->x - rr); x <= MIN (W - 1, (gint) st->x + rr); x++)
          {
            a += out[(y * W + x) * 4 + 1] - f->sky;
            b += f->image[(y * W + x) * 4 + 1] - f->sky;
          }
      fl[n] = a / b;
      n++;
    }
  qsort (r, n, sizeof (gdouble), cmp_double_asc);
  qsort (fl, n, sizeof (gdouble), cmp_double_asc);
  med = r[n / 2];
  if (flux_ratio)
    *flux_ratio = fl[n / 2];
  return med;
}

static void
check_star_reduce_fwhm (void)
{
  StarField *f = star_field (41, 24, 0.002, 0.0, 0.0);
  static const struct { const gchar *label; const gchar *method; gint m; gdouble amount, radius; gint it; gdouble limit; } runs[] =
    {
      { "default", "selection 0.25", 1, 70, 1.5, 1, 0.92 },
      { "erosion_100", "erosion", 0, 100, 1.5, 1, 0.85 },
      { "erosion_100_2_iterations", "erosion", 0, 100, 1.5, 2, 0.80 },
    };
  guint   k;
  gdouble prev = 1.0;

  for (k = 0; k < G_N_ELEMENTS (runs); k++)
    {
      gfloat *o = run (f->image, 0, 0, "astro:star-reduce", "method", runs[k].m,
                       "amount", runs[k].amount, "radius", runs[k].radius,
                       "iterations", runs[k].it, NULL);
      gdouble flux, ratio = fwhm_ratio (f, o, &flux);
      gchar  *name = g_strdup_printf ("star_reduce_fwhm_%s", runs[k].label);

      report (name, ratio < runs[k].limit,
              "%s, amount %.0f, radius %.1f, %d iteration%s: FWHM after / before %.3f "
              "(limit %.2f), flux after / before %.3f", runs[k].method, runs[k].amount,
              runs[k].radius, runs[k].it, runs[k].it > 1 ? "s" : "", ratio, runs[k].limit, flux);
      if (k == 2)
        report ("star_reduce_more_iterations_smaller", ratio < prev,
                "2 iterations %.3f, 1 iteration %.3f", ratio, prev);
      prev = ratio;
      g_free (name);
      g_free (o);
    }
  star_field_free (f);
}

/* outside the mask nothing changes, nebula included */
static void
check_star_reduce_nebula (void)
{
  StarField *f = star_field (42, 24, 0.003, 0.15, 0.0);
  gfloat    *o = run (f->image, 0, 0, "astro:star-reduce", "amount", 100.0, NULL);
  gfloat    *m = run (f->image, 0, 0, "astro:star-mask", NULL);
  gdouble    out_change = 0, neb_out = 0, all_change = 0, neb_all = 0, d = 0;
  gint       i, c;

  for (i = 0; i < W * H; i++)
    {
      gdouble ch = o[i * 4 + 1] - f->image[i * 4 + 1];

      neb_all    += f->nebula[i];
      all_change += f->stars[i] < 0.1 * f->noise ? ch : 0;
      if (m[i * 4] == 0.0f)
        {
          out_change += ch;
          neb_out    += f->nebula[i];
          for (c = 0; c < 4; c++)
            d = MAX (d, fabs (o[i * 4 + c] - f->image[i * 4 + c]));
        }
    }
  report ("star_reduce_outside_mask_unchanged", d == 0 && fabs (out_change) <= 0.001 * neb_out,
          "outside the mask: largest change %.2g, nebula flux changed by %.2g %%", d,
          100 * out_change / neb_out);
  printf ("NOTE  star_reduce_nebula: away from stars (under the mask's edges too) the "
          "nebula's flux changed by %.3f %%\n", 100 * all_change / neb_all);
  g_free (o);
  g_free (m);
  star_field_free (f);
}

static int
cmp_float_asc (const void *a,
               const void *b)
{
  gfloat x = *(const gfloat *) a, y = *(const gfloat *) b;

  return (x > y) - (x < y);
}

/* no dark rings: in the ring 2 to 6 star sigmas out, the median after
 * erosion at 100 % is not below the sky minus 3 noise sigmas (nearly no
 * noise, and noise of 0.003); the deepest ring is reported in noise
 * sigmas */
static void
check_star_reduce_rings (void)
{
  gdouble noises[2] = { 1e-4, 0.003 };
  gint    t;

  for (t = 0; t < 2; t++)
    {
      StarField *f = star_field (43 + t, 24, noises[t], 0.0, 0.0);
      gfloat    *o = run (f->image, 0, 0, "astro:star-reduce", "method", 0, "amount", 100.0, NULL);
      gfloat    *ring = g_new (gfloat, W * H);
      gint       i, bad = 0;
      gdouble    worst = -G_MAXDOUBLE;
      gchar     *name;

      for (i = 0; i < f->n; i++)
        {
          const Star *st = &f->star[i];
          gint        r = (gint) ceil (6 * st->sigma), x, y, cnt = 0;
          gdouble     depth;

          for (y = MAX (0, (gint) st->y - r); y <= MIN (H - 1, (gint) st->y + r); y++)
            for (x = MAX (0, (gint) st->x - r); x <= MIN (W - 1, (gint) st->x + r); x++)
              {
                gdouble dx = x + 0.5 - st->x, dy = y + 0.5 - st->y, d2 = dx * dx + dy * dy;

                if (d2 >= 4 * st->sigma * st->sigma && d2 <= r * r)
                  ring[cnt++] = o[(y * W + x) * 4 + 1];
              }
          qsort (ring, cnt, sizeof (gfloat), cmp_float_asc);
          depth = (f->sky - ring[cnt / 2]) / f->noise;
          worst = MAX (worst, depth);
          bad += depth > 3;
        }
      name = g_strdup_printf ("star_reduce_no_dark_rings_noise_%g", noises[t]);
      report (name, bad == 0,
              "erosion 100 %%: %d of %d stars with the ring's median below sky - 3 noise; "
              "deepest ring median %.2f noise sigmas %s the sky", bad, f->n, fabs (worst),
              worst > 0 ? "below" : "above");
      g_free (name);
      g_free (ring);
      g_free (o);
      star_field_free (f);
    }
}

static void
check_star_reduce_misc (void)
{
  StarField *f = star_field (45, 22, 0.003, 0.05, 0.3);
  gfloat    *a = run (f->image, 0, 0, "astro:star-reduce", NULL);
  gfloat    *b = run (f->image, 52, -17, "astro:star-reduce", NULL);
  gfloat    *z = run (f->image, 0, 0, "astro:star-reduce", "amount", 0.0, NULL);
  gdouble    d = 0, ratio = 0;
  gint       i, changed = 0;

  report ("star_reduce_offset_invariant", memcmp (a, b, sizeof (gfloat) * W * H * 4) == 0,
          "buffer at 52, -17 %s", memcmp (a, b, sizeof (gfloat) * W * H * 4) ? "differs" : "identical");
  for (i = 0; i < W * H * 4; i++)
    d = MAX (d, fabs (z[i] - f->image[i]));
  report ("star_reduce_amount_0_unchanged", d == 0, "max difference %.2g", d);
  /* keep star colour: R / G the same where the pixel changed */
  for (i = 0; i < W * H; i++)
    if (a[i * 4 + 1] != f->image[i * 4 + 1] && f->image[i * 4 + 1] > 0.02f)
      {
        ratio = MAX (ratio, fabs ((a[i * 4] / a[i * 4 + 1]) /
                                  (f->image[i * 4] / f->image[i * 4 + 1]) - 1));
        changed++;
      }
  report ("star_reduce_colour_kept", ratio < 1e-5 && changed > 0,
          "%d changed pixels: R / G kept to %.2g (relative)", changed, ratio);
  g_free (a);
  g_free (b);
  g_free (z);
  star_field_free (f);
}
