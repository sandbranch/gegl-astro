# Probe (not a test): astro:background's polynomial and spline surfaces
# on synthetic scenes (dome, band, nebula), max and mean error of the
# model as a share of the sky's amplitude. Run with GEGL_PATH pointing at
# a folder with astro-background.so and GEGL's own operations:
#   python3 tests/probes/background-scenes.py dome|band|nebula
# Needs numpy and the GEGL introspection data.
#
# Copyright 2026 David
# SPDX-License-Identifier: GPL-3.0-or-later
import sys, gi
import numpy as np
gi.require_version('Gegl', '0.4')
from gi.repository import Gegl
Gegl.init(None)
W, H = 600, 400
y, x = np.mgrid[0:H, 0:W]
u = 2 * (x + 0.5) / W - 1; v = 2 * (y + 0.5) / H - 1
kind = sys.argv[1]
rng = np.random.default_rng(1)
if kind == 'dome':
    dx = (u + 1.3) * 0.5; dy = (v - 1.3) * 0.5 * H / W
    sky = 0.045 + 0.07 * np.exp(-np.sqrt(dx * dx + dy * dy) / 0.3)
    extra = 0
elif kind == 'band':
    sky = 0.05 + 0.01 * u + 0.06 * np.exp(-(H - y - 0.5) / (0.15 * H))
    extra = 0
elif kind == 'nebula':
    sky = 0.050 + 0.025 * u + 0.020 * v + 0.010 * u * u - 0.006 * u * v + 0.008 * v * v
    extra = 0.06 * np.exp(-((x - 0.62 * W) ** 2 + (y - 0.45 * H) ** 2) / (2 * 55.0 ** 2))
img = np.zeros((H, W, 4), np.float32)
for c in range(3):
    img[..., c] = sky + extra + rng.normal(0, 0.003, (H, W))
img[..., 3] = 1
def run(props):
    buf = Gegl.Buffer.new("RGBA float", 0, 0, W, H)
    buf.set(Gegl.Rectangle.new(0, 0, W, H), "RGBA float", img.tobytes())
    g = Gegl.Node(); src = g.create_child('gegl:buffer-source'); src.set_property('buffer', buf)
    op = g.create_child('astro:background')
    for k, val in props.items(): op.set_property(k, val)
    out = Gegl.Buffer.new("RGBA float", 0, 0, W, H)
    sink = g.create_child('gegl:write-buffer'); sink.set_property('buffer', out)
    src.link(op); op.link(sink); sink.process()
    return np.frombuffer(out.get(Gegl.Rectangle.new(0, 0, W, H), 1.0, "RGBA float", Gegl.AbyssPolicy.NONE), np.float32).reshape(H, W, 4)
amp = sky.max() - sky.min()
for props in [dict(show=1, degree=4), dict(show=1, degree=4, tolerance=10.0), dict(show=1, model=1), dict(show=1, model=1, tolerance=10.0), dict(show=1, model=1, tolerance=10.0, samples=32)]:
    m = run(props)[..., 1]
    e = np.abs(m - sky); i = np.unravel_index(e.argmax(), e.shape)
    print(props, 'max %.2f %%' % (100 * e.max() / amp), 'at', tuple(int(t) for t in i), 'mean %.3f %%' % (100 * e.mean() / amp))
