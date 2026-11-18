# Probe (not a test): astro:background on a noise-free quadratic sky;
# where the model's largest error is. GEGL_PATH as background-scenes.py.
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
sky = 0.060 + 0.030 * u + 0.025 * v + 0.012 * u * u - 0.008 * u * v + 0.010 * v * v
img = np.zeros((H, W, 4), np.float32); img[..., 0] = img[..., 1] = img[..., 2] = sky; img[..., 3] = 1
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
for props in [dict(show=1), dict(show=1, model=1), dict(show=1, model=1, tolerance=10.0)]:
    m = run(props)[..., 1]
    e = np.abs(m - sky); i = np.unravel_index(e.argmax(), e.shape)
    print(props, 'max err %.2f %%' % (100 * e.max() / amp), 'at', i, 'mean err %.3f %%' % (100 * e.mean() / amp))
