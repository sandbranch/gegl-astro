# Probe (not a test): the polynomial against the spline with the robust
# cell rejection of spline_fit_robust (ASTRO_SPLINE_* environment
# variables pick the variant), on the synthetic scenes and the real
# images of tests/real.py. Run like tests/real.sh, with GEGL_PATH at this
# branch's modules: python3 tests/probes/background-spline-placement.py tests
#
# Copyright 2026 David
# SPDX-License-Identifier: GPL-3.0-or-later
import os, sys
sys.path.insert(0, sys.argv[1])
import real, numpy as np
VARIANTS = [
    ('poly deg 2', {'degree': 2}, {}),
    ('poly deg 4', {'degree': 4}, {}),
    ('spline robust', {'model': 1}, {}),
    ('spline robust above', {'model': 1}, {'ASTRO_SPLINE_ABOVE': '1'}),
    ('spline robust s0.8', {'model': 1, 'smoothing': 0.8}, {}),
    ('spline old (local)', {'model': 1}, {'ASTRO_SPLINE_LOCAL': '1', 'ASTRO_SPLINE_ROUNDS': '0'}),
]
def setenv(env):
    for k in ('ASTRO_SPLINE_ABOVE', 'ASTRO_SPLINE_LOCAL', 'ASTRO_SPLINE_ROUNDS'):
        os.environ.pop(k, None)
    os.environ.update(env)
print('columns:', [v[0] for v in VARIANTS])
print('synthetic: max / mean error, % of the sky range (tolerance 5)')
for kind in ['quadratic', 'dome', 'band', 'nebula', 'big-nebula', 'galaxy']:
    img, sky, extra = real.synthetic_scene(kind)
    amp = sky.max() - sky.min(); row = '%-11s' % kind
    for name, props, env in VARIANTS:
        setenv(env)
        m = real.run('astro:background', dict(props, show='model'), img)
        e = np.abs(m[..., 1] - sky) / amp
        row += ' | %-6.1f/%5.2f' % (100 * e.max(), 100 * e.mean())
    print(row, flush=True)
print('real + glow: max error quadratic / dome / band, and bend')
for im in real.manifest():
    a = real.load(im['file'], half=True)
    mask = real.star_mask(a)
    amp = max(float(np.median(real.sky_boxes(a, mask))), 0.02)
    h, w = a.shape[:2]; grads = real.gradients(h, w)
    row = '%-14s' % im['file'][:14]
    for name, props, env in VARIANTS:
        setenv(env)
        m0 = real.run('astro:background', dict(props, show='model'), a)
        lm = real.lum(m0)
        bend = 255 * (real.srgb(np.percentile(lm, 99)) - real.srgb(np.percentile(lm, 1)))
        errs = [real.model_error(a, m0, g, amp, props)[0] for g in grads.values()]
        row += ' | %4.0f %4.0f %4.0f b%4.1f' % (*errs, bend)
    print(row, flush=True)
