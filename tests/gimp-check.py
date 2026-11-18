# Runs inside GIMP (tests/gimp-check.sh): applies each operation to a
# synthetic star field as a non-destructive filter (Gimp.DrawableFilter),
# checks that it is one, saves the image as XCF with the filter and loads
# it again (the settings must come back), merges the filters and compares
# the pixels with the same operation run in plain GEGL in this process,
# on RGB images of every precision, on a grayscale image and on a layer
# smaller than the image with offsets. Writes the scene and GIMP's float
# results as TIFF for the comparison with the gegl command line in
# gimp-check.sh. Prints PASS or FAIL per check. (After gegl-presence's
# tests/gimp-check.py.)
#
# Copyright 2026 David
# SPDX-License-Identifier: GPL-3.0-or-later
import array
import json
import math
import os
import random

import gi
gi.require_version('Gimp', '3.0')
gi.require_version('Gegl', '0.4')
from gi.repository import Gimp, Gegl, Gio, GLib  # noqa: E402

OUT = os.environ['ASTRO_CHECK_OUT']
W, H = 256, 160
FMT = "R'G'B'A float"
failed = []


def result(name, ok, msg=''):
    print('%s  gimp_%s%s' % ('PASS' if ok else 'FAIL', name, (': ' + msg) if msg else ''))
    if not ok:
        failed.append(name)


# the enums of the operations, in the order of their values; GIMP's
# config takes them by nick, GEGL by number
ENUMS = {
    'astro:background': {'correction': ['subtract', 'divide'],
                         'show': ['corrected', 'model']},
    'astro:stretch': {'curve': ['ghs', 'arcsinh'],
                      'colour': ['preserve', 'per-channel']},
    'astro:star-mask': {'show': ['mask', 'overlay']},
    'astro:star-reduce': {'method': ['erosion', 'selection']},
    'astro:star-color': {},
}

TITLE = 'Remove Sky Background'
CASES = [
    ('background', 'astro:background', TITLE, {'degree': 2}),
    ('background_divide', 'astro:background', TITLE,
     {'correction': 'divide', 'tolerance': 1.5}),
    ('background_model', 'astro:background', TITLE,
     {'show': 'model', 'degree': 4, 'samples': 10}),
    ('background_no_level', 'astro:background', TITLE,
     {'keep-level': False, 'degree': 1}),
    ('background_exclude_stars', 'astro:background', TITLE,
     {'exclude-stars': True, 'samples': 8}),
    ('stretch_ghs', 'astro:stretch', 'Astro Stretch',
     {'stretch': 4.0, 'local-intensity': 2.0, 'symmetry-point': 0.1,
      'protect-shadows': 0.02, 'protect-highlights': 0.9}),
    ('stretch_arcsinh', 'astro:stretch', 'Astro Stretch',
     {'curve': 'arcsinh', 'beta': 50.0, 'colour': 'per-channel', 'black-point': 0.01}),
    ('stretch_auto', 'astro:stretch', 'Astro Stretch',
     {'auto-points': True, 'stretch': 3.0, 'local-intensity': 1.0}),
    ('stretch_inverse', 'astro:stretch', 'Astro Stretch',
     {'stretch': 2.0, 'local-intensity': -1.0, 'symmetry-point': 0.2, 'inverse': True}),
    ('star_mask', 'astro:star-mask', 'Star Mask', {'threshold': 5.0}),
    ('star_mask_overlay', 'astro:star-mask', 'Star Mask',
     {'show': 'overlay', 'grow': 2.0, 'max-size': 10.0}),
    ('star_reduce', 'astro:star-reduce', 'Reduce Stars', {'amount': 80.0}),
    ('star_reduce_erosion', 'astro:star-reduce', 'Reduce Stars',
     {'method': 'erosion', 'iterations': 2, 'luminance-only': False, 'feather': 2.5}),
    ('star_color', 'astro:star-color', 'Star Colour', {'saturation': 120.0}),
    ('star_color_no_cores', 'astro:star-color', 'Star Colour',
     {'fill-cores': False, 'background': 0.0, 'threshold': 6.0}),
]


def scene():
    """a star field: a sky gradient, Gaussian stars, a little noise; the
    top rows half transparent"""
    rnd = random.Random(1)
    sky = []
    for y in range(H):
        for x in range(W):
            u, v = 2 * (x + 0.5) / W - 1, 2 * (y + 0.5) / H - 1
            base = 0.08 + 0.04 * u + 0.03 * v + 0.02 * u * u
            sky.append([base * 1.2, base, base * 0.85])
    for _ in range(150):
        sx, sy = rnd.uniform(0, W), rnd.uniform(0, H)
        sig = rnd.uniform(0.8, 2.0)
        peak = min(0.9, 0.03 * rnd.uniform(1e-3, 1) ** -0.6)
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
    for i, p in enumerate(sky):
        alpha = 0.5 if i < 8 * W else 1.0
        a.extend([min(max(q + (rnd.random() - 0.5) * 0.01, 0), 1) for q in p] + [alpha])
    return a


def gegl_props(op, props):
    out = {}
    for k, v in props.items():
        if k in ENUMS[op]:
            out[k] = ENUMS[op][k].index(v)
        else:
            out[k] = v
    return out


def gegl_result(source, op, props, width, height):
    """the operation on the buffer source in plain GEGL, as R'G'B'A float"""
    rect = Gegl.Rectangle.new(0, 0, width, height)
    g = Gegl.Node()
    src = g.create_child('gegl:buffer-source')
    src.set_property('buffer', source)
    node = g.create_child(op)
    for k, v in gegl_props(op, props).items():
        node.set_property(k, v)
    out = Gegl.Buffer.new(FMT, 0, 0, width, height)
    sink = g.create_child('gegl:write-buffer')
    sink.set_property('buffer', out)
    src.link(node)
    node.link(sink)
    sink.process()
    return out


def quantized(buf, image, like):
    """buf's pixels stored in a layer of the same type as like, and read
    back as R'G'B'A float"""
    rect = Gegl.Rectangle.new(0, 0, W, H)
    q = Gimp.Layer.new(image, 'q', W, H, like.type(), 100, Gimp.LayerMode.NORMAL)
    image.insert_layer(q, None, 1)
    qb = q.get_buffer()
    qb.set(rect, FMT, buf.get(rect, 1.0, FMT, Gegl.AbyssPolicy.NONE))
    qb.flush()
    out = pixels(q)
    image.remove_layer(q)
    return out


def save_tiff(path, data):
    rect = Gegl.Rectangle.new(0, 0, W, H)
    buf = Gegl.Buffer.new(FMT, 0, 0, W, H)
    buf.set(rect, FMT, data.tobytes())
    g = Gegl.Node()
    src = g.create_child('gegl:buffer-source')
    src.set_property('buffer', buf)
    save = g.create_child('gegl:tiff-save')
    save.set_property('path', path)
    src.link(save)
    save.process()


def pixels(drawable, fmt=FMT):
    rect = Gegl.Rectangle.new(0, 0, drawable.get_width(), drawable.get_height())
    return array.array('f', drawable.get_buffer().get(rect, 1.0, fmt, Gegl.AbyssPolicy.NONE))


def max_diff(a, b):
    if len(a) != len(b):
        return float('inf')
    return max(abs(p - q) for p, q in zip(a, b))


def new_image(precision, data, base=Gimp.ImageBaseType.RGB, alpha=True):
    image = Gimp.Image.new_with_precision(W, H, base, precision)
    if base == Gimp.ImageBaseType.RGB:
        kind = Gimp.ImageType.RGBA_IMAGE if alpha else Gimp.ImageType.RGB_IMAGE
    else:
        kind = Gimp.ImageType.GRAYA_IMAGE if alpha else Gimp.ImageType.GRAY_IMAGE
    layer = Gimp.Layer.new(image, 'scene', W, H, kind, 100, Gimp.LayerMode.NORMAL)
    image.insert_layer(layer, None, 0)
    rect = Gegl.Rectangle.new(0, 0, W, H)
    buf = layer.get_buffer()
    buf.set(rect, FMT, data.tobytes())
    buf.flush()
    layer.update(0, 0, W, H)
    return image, layer


def add_filter(drawable, op, title, props):
    f = Gimp.DrawableFilter.new(drawable, op, title)
    cfg = f.get_config()
    for k, v in props.items():
        cfg.set_property(k, v)
    f.update()
    drawable.append_filter(f)
    return f


def same_setting(got, want):
    if isinstance(want, float):
        # XCF keeps the filter's numbers as 32 bit floats
        return abs(got - want) < 1e-5
    if isinstance(want, str) and want.startswith('#'):
        a = got.get_rgba()
        b = Gegl.Color.new(want).get_rgba()
        return all(abs(x - y) < 1e-6 for x, y in zip(a, b))
    return got == want


def tol_for(precision):
    name = precision.value_nick
    if name.startswith('u8'):
        return 1.0 / 255 + 1e-6
    if name.startswith('u16'):
        return 2.0 / 65535
    if name.startswith('half'):
        return 2e-3
    return 2e-5


Gegl.init(None)
data = scene()
save_tiff(os.path.join(OUT, 'scene.tif'), data)
settings_out = {}

PRECISIONS = [Gimp.Precision.U8_NON_LINEAR, Gimp.Precision.U16_NON_LINEAR,
              Gimp.Precision.U32_LINEAR, Gimp.Precision.HALF_LINEAR,
              Gimp.Precision.FLOAT_NON_LINEAR, Gimp.Precision.FLOAT_LINEAR,
              Gimp.Precision.DOUBLE_NON_LINEAR]

for label, op, title, props in CASES:
    settings_out[label] = {'op': op, 'props': props}
    for precision in PRECISIONS + ['gray']:
        if precision == 'gray':
            name = label + '_gray_8_bit'
            image, layer = new_image(Gimp.Precision.U8_NON_LINEAR, data,
                                     Gimp.ImageBaseType.GRAY)
            prec = Gimp.Precision.U8_NON_LINEAR
        else:
            name = label + '_' + precision.value_nick.replace('-', '_')
            image, layer = new_image(precision, data)
            prec = precision
        # plain GEGL on the layer's pixels, stored in the layer's format
        want_buf = gegl_result(layer.get_buffer(), op, props, W, H)
        want = quantized(want_buf, image, layer)
        f = add_filter(layer, op, title, props)
        names = [x.get_operation_name() for x in layer.get_filters()]
        ok_filter = names == [op]

        # the filter survives an XCF save and load, with its settings
        xcf = os.path.join(OUT, name + '.xcf')
        Gimp.file_save(Gimp.RunMode.NONINTERACTIVE, image, Gio.File.new_for_path(xcf), None)
        loaded = Gimp.file_load(Gimp.RunMode.NONINTERACTIVE, Gio.File.new_for_path(xcf))
        ll = loaded.get_layers()[0]
        lf = ll.get_filters()
        same = len(lf) == 1 and lf[0].get_operation_name() == op
        bad = []
        if same:
            cfg = lf[0].get_config()
            for k, v in props.items():
                got = cfg.get_property(k)
                if not same_setting(got, v):
                    same = False
                    bad.append('%s=%r' % (k, got))

        layer.merge_filters()
        got = pixels(layer)
        d = max_diff(got, want)
        ll.merge_filters()
        d2 = max_diff(pixels(ll), want)
        tol = tol_for(prec)
        ok = ok_filter and same and d <= tol and d2 <= tol
        result(name, ok, 'non-destructive %s, settings after XCF %s, against plain GEGL '
               '%.2g, after XCF %.2g' % (ok_filter, 'kept' if same else 'LOST ' + ', '.join(bad),
                                         d, d2))
        if precision == Gimp.Precision.FLOAT_NON_LINEAR:
            save_tiff(os.path.join(OUT, 'gimp-' + label + '.tif'), got)
        image.delete()
        loaded.delete()

# a layer smaller than the image, with offsets: the same as plain GEGL
# on the layer's own pixels
for label, op, title, props in CASES[::2]:
    image = Gimp.Image.new_with_precision(W + 60, H + 40, Gimp.ImageBaseType.RGB,
                                          Gimp.Precision.FLOAT_LINEAR)
    layer = Gimp.Layer.new(image, 'scene', W, H, Gimp.ImageType.RGBA_IMAGE, 100,
                           Gimp.LayerMode.NORMAL)
    image.insert_layer(layer, None, 0)
    layer.set_offsets(37, 21)
    buf = layer.get_buffer()
    buf.set(Gegl.Rectangle.new(0, 0, W, H), FMT, data.tobytes())
    buf.flush()
    want = quantized(gegl_result(layer.get_buffer(), op, props, W, H), image, layer)
    add_filter(layer, op, title, props)
    layer.merge_filters()
    d = max_diff(pixels(layer), want)
    result(label + '_layer_with_offsets', d <= 2e-5,
           'layer at 37, 21 in a larger image, against plain GEGL %.2g' % d)
    image.delete()

with open(os.path.join(OUT, 'cases.json'), 'w') as fh:
    json.dump(settings_out, fh)
print('%d failed' % len(failed))
with open(os.path.join(OUT, 'gimp-check.status'), 'w') as fh:
    fh.write('%d\n' % len(failed))
