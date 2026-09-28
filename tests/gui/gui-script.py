# Runs inside GIMP on a Broadway display (tests/gui/start.sh). Opens the
# image ASTRO_GUI_IMAGE (generated, or CC0 or public domain; never
# private photos), scaled to 1200 pixels across, or without it a
# generated star field with a sky gradient, and shows it in a window with
# its layer selected. (After gegl-presence's tests/gui/gui-script.py.)
#
# Then waits: when ASTRO_GUI_OUT/save exists it saves the image as
# ASTRO_GUI_OUT/gui-<name>.xcf (ASTRO_GUI_NAME) and writes
# ASTRO_GUI_OUT/saved; when ASTRO_GUI_OUT/quit exists (or after 10
# minutes) it quits GIMP.
#
# Copyright 2026 David
# SPDX-License-Identifier: GPL-3.0-or-later
import array
import math
import os
import random
import time

import gi
gi.require_version('Gimp', '3.0')
gi.require_version('Gegl', '0.4')
from gi.repository import Gimp, Gegl, Gio  # noqa: E402

out = os.environ['ASTRO_GUI_OUT']
name = os.environ.get('ASTRO_GUI_NAME', 'gui')
photo = os.environ.get('ASTRO_GUI_IMAGE', '')


def picture(W, H):
    """a star field on a sky that brightens towards the bottom right"""
    rnd = random.Random(7)
    sky = []
    for y in range(H):
        for x in range(W):
            u, v = x / (W - 1), y / (H - 1)
            b = 0.05 + 0.12 * u * v + 0.04 * v
            sky.append([b * 1.25, b * 1.05, b * 0.8])
    for _ in range(700):
        sx, sy = rnd.uniform(0, W), rnd.uniform(0, H)
        sig = rnd.uniform(0.7, 2.2)
        peak = min(1.0, 0.05 * rnd.uniform(1e-3, 1) ** -0.6)
        t = rnd.random()
        col = (0.8 + 0.4 * t, 1.0, 1.2 - 0.4 * t)
        r = int(math.ceil(4 * sig))
        for y in range(max(0, int(sy) - r), min(H, int(sy) + r + 1)):
            for x in range(max(0, int(sx) - r), min(W, int(sx) + r + 1)):
                g = peak * math.exp(-((x + 0.5 - sx) ** 2 + (y + 0.5 - sy) ** 2) / (2 * sig * sig))
                p = sky[y * W + x]
                for c in range(3):
                    p[c] += g * col[c]
    a = array.array('f')
    for p in sky:
        a.extend([min(q, 1.0) for q in p] + [1.0])
    return a


if photo and os.path.exists(photo):
    image = Gimp.file_load(Gimp.RunMode.NONINTERACTIVE, Gio.File.new_for_path(photo))
    w, h = image.get_width(), image.get_height()
    if w > 1200:
        image.scale(1200, round(h * 1200 / w))
    layer = image.get_layers()[0]
else:
    W, H = 900, 600
    image = Gimp.Image.new(W, H, Gimp.ImageBaseType.RGB)
    layer = Gimp.Layer.new(image, 'picture', W, H, Gimp.ImageType.RGBA_IMAGE, 100,
                           Gimp.LayerMode.NORMAL)
    image.insert_layer(layer, None, 0)
    buf = layer.get_buffer()
    buf.set(Gegl.Rectangle.new(0, 0, W, H), "R'G'B'A float", picture(W, H).tobytes())
    buf.flush()
    layer.update(0, 0, W, H)
image.set_selected_layers([layer])
image.clean_all()
display = Gimp.Display.new(image)
Gimp.displays_flush()
with open(os.path.join(out, 'ready'), 'w') as f:
    f.write('ready\n')

t0 = time.time()
while not os.path.exists(os.path.join(out, 'quit')) and time.time() - t0 < 600:
    if os.path.exists(os.path.join(out, 'save')):
        os.remove(os.path.join(out, 'save'))
        xcf = os.path.join(out, 'gui-%s.xcf' % name)
        Gimp.file_save(Gimp.RunMode.NONINTERACTIVE, image, Gio.File.new_for_path(xcf), None)
        with open(os.path.join(out, 'saved'), 'w') as f:
            f.write(xcf + '\n')
    time.sleep(0.3)

pdb = Gimp.get_pdb()
quit_proc = pdb.lookup_procedure('gimp-quit')
quit_config = quit_proc.create_config()
quit_config.set_property('force', True)
quit_proc.run(quit_config)
