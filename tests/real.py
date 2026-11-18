#!/usr/bin/env python3
# The operations on the real images of tests/images (CC0, fetched into
# tests/output/images by tests/images/fetch.py) and on synthetic scenes:
# numbers to decide defaults with, and JPEGs to look at. Not a pass or
# fail test: real images have no known truth, so these are measurements.
# Run through tests/real.sh (GEGL_PATH with the modules of
# tests/output/build-check, a throwaway home). Needs numpy and scipy and
# the GEGL introspection data. Writes tests/output/real/<part>.txt and
# JPEGs there.
#
#   real.py background   sky residual spread on the star-masked sky,
#                        before and after, and before / after JPEGs
#   real.py tolerance    the tolerance study: synthetic scenes and the
#                        real images with a known gradient added
#   real.py stars        star mask counts and overlay JPEGs
#   real.py reduce       Reduce Stars: rings around stars, crops
#   real.py stretch      Astro Stretch, linear and encoded, JPEGs
#
# The JPEGs are 8 bit sRGB; the images are decoded to linear light (the
# sRGB curve) as GIMP does for the filters. These images are finished
# (stretched, gradients mostly removed by their authors), so their
# linear values are not raw data: background removal on them should
# change little, and a gradient added to them has a known truth.
#
# Copyright 2026 David
# SPDX-License-Identifier: GPL-3.0-or-later
import json
import os
import sys

import gi
import numpy as np
from scipy import ndimage

gi.require_version('Gegl', '0.4')
from gi.repository import Gegl  # noqa: E402

Gegl.init(None)
here = os.path.dirname(os.path.abspath(__file__))
images_dir = os.path.join(here, 'output', 'images')
out_dir = os.path.join(here, 'output', 'real')
os.makedirs(out_dir, exist_ok=True)
FMT = 'RGBA float'
LUM = np.array([0.2126, 0.7152, 0.0722], np.float32)   # sRGB primaries
log_file = None


def say(*args):
    line = ' '.join(str(a) for a in args)
    print(line)
    if log_file:
        log_file.write(line + '\n')
        log_file.flush()


def to_buffer(arr, fmt=FMT):
    h, w = arr.shape[:2]
    buf = Gegl.Buffer.new(fmt, 0, 0, w, h)
    buf.set(Gegl.Rectangle.new(0, 0, w, h), fmt, np.ascontiguousarray(arr, np.float32).tobytes())
    return buf


def from_buffer(buf, w, h, fmt=FMT):
    data = buf.get(Gegl.Rectangle.new(0, 0, w, h), 1.0, fmt, Gegl.AbyssPolicy.NONE)
    return np.frombuffer(data, np.float32).reshape(h, w, 4).copy()


def load(name, half=False):
    """the image as linear RGBA float (sRGB primaries), cut to the
    manifest's "crop" (x, y, w, h) when it has one, halved by 2 x 2 means
    with half"""
    g = Gegl.Node()
    src = g.create_child('gegl:load')
    src.set_property('path', os.path.join(images_dir, name))
    r = src.get_bounding_box()
    out = Gegl.Buffer.new(FMT, 0, 0, r.width, r.height)
    sink = g.create_child('gegl:write-buffer')
    sink.set_property('buffer', out)
    src.link(sink)
    sink.process()
    a = from_buffer(out, r.width, r.height)
    crop = next((im.get('crop') for im in manifest() if im['file'] == name), None)
    if crop:
        x, y, w, h = crop
        a = a[y:y + h, x:x + w].copy()
    if half:
        h, w = a.shape[0] // 2 * 2, a.shape[1] // 2 * 2
        a = a[:h, :w].reshape(h // 2, 2, w // 2, 2, 4).mean(axis=(1, 3))
    return a


def run(op, props, arr, fmt=FMT):
    h, w = arr.shape[:2]
    buf = to_buffer(arr, fmt)
    g = Gegl.Node()
    src = g.create_child('gegl:buffer-source')
    src.set_property('buffer', buf)
    node = g.create_child(op)
    for k, v in props.items():
        node.set_property(k, v)
    out = Gegl.Buffer.new(fmt, 0, 0, w, h)
    sink = g.create_child('gegl:write-buffer')
    sink.set_property('buffer', out)
    src.link(node)
    node.link(sink)
    sink.process()
    return from_buffer(out, w, h, fmt)


def save_jpeg(arr, name, width=1600, encoded=False):
    """arr (linear RGBA float, or R'G'B'A with encoded) as an sRGB JPEG,
    scaled down to width by area means"""
    a = np.clip(np.nan_to_num(arr), 0, 1).astype(np.float32)
    if width and a.shape[1] > width:
        f = int(np.ceil(a.shape[1] / width))
        h, w = a.shape[0] // f * f, a.shape[1] // f * f
        a = a[:h, :w].reshape(h // f, f, w // f, f, 4).mean(axis=(1, 3))
    a[..., 3] = 1
    h, w = a.shape[:2]
    g = Gegl.Node()
    src = g.create_child('gegl:buffer-source')
    src.set_property('buffer', to_buffer(a, "R'G'B'A float" if encoded else FMT))
    sink = g.create_child('gegl:jpg-save')
    sink.set_property('path', os.path.join(out_dir, name))
    sink.set_property('quality', 92)
    src.link(sink)
    sink.process()


def side_by_side(a, b, gap=8):
    h = max(a.shape[0], b.shape[0])
    pad = lambda x: np.pad(x, ((0, h - x.shape[0]), (0, 0), (0, 0)))
    sep = np.ones((h, gap, 4), np.float32)
    return np.concatenate([pad(a), sep, pad(b)], axis=1)


def lum(a):
    return a[..., :3] @ LUM


def manifest():
    return json.load(open(os.path.join(here, 'images', 'manifest.json')))


def star_mask(a):
    return run('astro:star-mask', {}, a)[..., 0]


def sky_boxes(a, mask, boxes=24):
    """median luminance of the unmasked sky in boxes (boxes along the
    longer side, square); boxes with under half of their pixels free of
    stars are left out"""
    L = lum(a)
    h, w = L.shape
    size = max(h, w) / boxes
    meds = []
    for by in range(int(np.ceil(h / size))):
        for bx in range(int(np.ceil(w / size))):
            ys = slice(int(by * size), min(h, int((by + 1) * size)))
            xs = slice(int(bx * size), min(w, int((bx + 1) * size)))
            free = mask[ys, xs] <= 0
            if free.mean() >= 0.5:
                meds.append(np.median(L[ys, xs][free]))
    return np.array(meds)


def srgb(x):
    """the sRGB curve: linear to encoded"""
    x = np.clip(x, 0, None)
    return np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(x, 1 / 2.4) - 0.055)


def spread(meds):
    """P5 to P95 and P25 to P75 of the box medians, in 8 bit levels of the
    encoded (sRGB) luminance, which is how evenness looks on screen"""
    p5, p25, p75, p95 = 255 * srgb(np.percentile(meds, [5, 25, 75, 95]))
    return p95 - p5, p75 - p25


# ---- background ----------------------------------------------------------

def part_background():
    say('Remove Sky Background on the real images (half size, linear): spread')
    say('of the medians of the star-free sky in 24 boxes along the longer side,')
    say('P5 to P95 (and P25 to P75), in 8 bit levels of the encoded luminance;')
    say('before, and after the defaults and other settings. The images are')
    say('finished, so much of the spread is nebulae, dust and the Milky Way,')
    say('not gradient: a filter that flattens those has gone too far.')
    say()
    settings = [('tolerance 2', {'tolerance': 2.0}), ('tolerance 5 (default)', {}),
                ('degree 4, tolerance 5', {'degree': 4, 'tolerance': 5.0}),
                ('stars left out', {'exclude-stars': True})]
    say('%-32s %-16s' % ('image', 'before') + ''.join('%-24s' % s[0] for s in settings))
    for im in manifest():
        a = load(im['file'], half=True)
        mask = star_mask(a)
        meds = sky_boxes(a, mask)
        sky = np.median(meds)
        row = '%-32s %5.1f (%4.1f)     ' % (im['file'], *spread(meds))
        outs = {}
        for label, props in settings:
            o = run('astro:background', props, a)
            outs[label] = o
            row += '%5.1f (%4.1f)           ' % spread(sky_boxes(o, mask))
        say(row)
        base = os.path.splitext(im['file'])[0]
        # before | tolerance 2 | tolerance 5, the same display stretch for all
        # (the sky median to 0.25 in the encoded values; finished images
        # are already bright, so a mild boost to see the sky's evenness)
        k = boost_factor(sky)
        save_jpeg(side_by_side(side_by_side(a * k, outs['tolerance 2'] * k),
                               outs['tolerance 5 (default)'] * k), base + '-background.jpg', 3000)
        # the image | the model at tolerance 2 | at tolerance 5 (degree 2)
        m2 = run('astro:background', {'show': 'model', 'tolerance': 2.0}, a)
        m5 = run('astro:background', {'show': 'model', 'tolerance': 5.0}, a)
        save_jpeg(side_by_side(side_by_side(a * k, m2 * k), m5 * k),
                  base + '-background-model.jpg', 3000)
    say()


def boost_factor(sky):
    """a factor that brings the linear sky median to 0.05 linear (0.25
    in sRGB), at most 4 (only to look at the sky)"""
    return float(np.clip(0.05 / max(sky, 1e-6), 1.0, 4.0))


# ---- tolerance study -----------------------------------------------------

def grid(h, w):
    y, x = np.mgrid[0:h, 0:w]
    return 2 * (x + 0.5) / w - 1, 2 * (y + 0.5) / h - 1, x, y


def gradients(h, w):
    """known sky gradients, each with a range of 1 (min 0, max 1)"""
    u, v, x, y = grid(h, w)
    out = {}
    q = 0.5 * u + 0.3 * v + 0.25 * u * u - 0.15 * u * v + 0.2 * v * v
    out['quadratic'] = q
    dx = (u + 1.3) * 0.5
    dy = (v - 1.3) * 0.5 * h / w
    out['dome'] = np.exp(-np.sqrt(dx * dx + dy * dy) / 0.3)
    out['band'] = 0.17 * u + np.exp(-(h - y - 0.5) / (0.15 * h))
    for k in out:
        g = out[k]
        out[k] = ((g - g.min()) / (g.max() - g.min())).astype(np.float32)
    return out


# light pollution's colour (sodium and LED glow is warm), per channel
GLOW = np.array([1.0, 0.8, 0.55], np.float32)


def model_error(a, m0, g, amp, props):
    """the model of a + glow minus m0, the model of a, against the glow:
    max and mean of the absolute error over the channels, % of the glow's
    range per channel"""
    added = a.copy()
    added[..., :3] += amp * g[..., None] * GLOW
    m1 = run('astro:background', dict(props, show='model'), added)
    e = np.abs(m1[..., :3] - m0[..., :3] - amp * g[..., None] * GLOW) / (amp * GLOW)
    return 100 * e.max(), 100 * e.mean()


def synthetic_scene(kind, h=400, w=600, seed=1):
    """the scenes of the RBF probes (PLAN op 1 item 2) with stars, and
    the truth"""
    rng = np.random.default_rng(seed)
    u, v, x, y = grid(h, w)
    extra = 0
    if kind == 'dome':
        dx = (u + 1.3) * 0.5
        dy = (v - 1.3) * 0.5 * h / w
        sky = 0.045 + 0.07 * np.exp(-np.sqrt(dx * dx + dy * dy) / 0.3)
    elif kind == 'band':
        sky = 0.05 + 0.01 * u + 0.06 * np.exp(-(h - y - 0.5) / (0.15 * h))
    elif kind == 'quadratic':
        sky = 0.060 + 0.030 * u + 0.025 * v + 0.012 * u * u - 0.008 * u * v + 0.010 * v * v
    elif kind in ('nebula', 'big-nebula'):
        sky = 0.050 + 0.025 * u + 0.020 * v + 0.010 * u * u - 0.006 * u * v + 0.008 * v * v
        s = 55.0 if kind == 'nebula' else 110.0
        extra = 0.06 * np.exp(-((x - 0.62 * w) ** 2 + (y - 0.45 * h) ** 2) / (2 * s * s))
    elif kind == 'galaxy':
        # a bright elongated galaxy across the middle (like M31 in a frame)
        sky = 0.050 + 0.02 * u - 0.015 * v + 0.01 * u * v
        xr = (x - 0.5 * w) * np.cos(0.6) + (y - 0.5 * h) * np.sin(0.6)
        yr = -(x - 0.5 * w) * np.sin(0.6) + (y - 0.5 * h) * np.cos(0.6)
        extra = 0.25 * np.exp(-np.sqrt((xr / 140.0) ** 2 + (yr / 35.0) ** 2) * 3)
    img = np.zeros((h, w, 4), np.float32)
    for c in range(3):
        img[..., c] = sky + extra
    # stars: Gaussians, few bright, many faint
    n = 600
    sx, sy = rng.uniform(0, w, n), rng.uniform(0, h, n)
    peak = 0.003 * 10 ** rng.uniform(0.5, 2.3, n)
    sig = rng.uniform(0.8, 2.0, n)
    for i in range(n):
        x0, x1 = int(max(0, sx[i] - 8)), int(min(w, sx[i] + 9))
        y0, y1 = int(max(0, sy[i] - 8)), int(min(h, sy[i] + 9))
        yy, xx = np.mgrid[y0:y1, x0:x1]
        img[y0:y1, x0:x1, :3] += (peak[i] * np.exp(-((xx - sx[i]) ** 2 + (yy - sy[i]) ** 2)
                                                 / (2 * sig[i] ** 2)))[..., None]
    img[..., :3] += rng.normal(0, 0.003, (h, w, 3)).astype(np.float32)
    img[..., 3] = 1
    return img, sky.astype(np.float32), np.asarray(extra, np.float32)


TOLERANCES = [1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0]


def part_tolerance():
    say('Tolerance study for Remove Sky Background.')
    say()
    say('A. Synthetic scenes (600 x 400, 600 stars, noise 0.003): model error')
    say('against the true sky, max / mean, % of the sky\'s range; for the')
    say('nebula and galaxy scenes the sky is the gradient without the object')
    say('(error = how far the model bends into the object).')
    say()
    say('%-12s %-3s ' % ('scene', 'deg') + ''.join('%-14s' % ('tol %g' % t) for t in TOLERANCES))
    for kind in ['quadratic', 'dome', 'band', 'nebula', 'big-nebula', 'galaxy']:
        img, sky, extra = synthetic_scene(kind)
        amp = sky.max() - sky.min()
        for deg in (2, 3, 4):
            row = '%-12s %-3d ' % (kind, deg)
            for t in TOLERANCES:
                m = run('astro:background', {'show': 'model', 'degree': deg, 'tolerance': t}, img)
                e = np.abs(m[..., 1] - sky) / amp
                row += '%5.1f / %4.2f  ' % (100 * e.max(), 100 * e.mean())
            say(row)
    say()
    say('B. Real images (half size, linear) with a known warm glow added (range')
    say('= the sky median, at least 0.02 linear; colour 1 : 0.8 : 0.55): the')
    say('model of image + glow')
    say('minus the model of the image, against the glow; max / mean error, %')
    say('of the glow\'s range. And "bend": the range (P1 to P99) of the model')
    say('of the image alone, in 8 bit levels of the encoded luminance (as on')
    say('screen); the images are finished')
    say('with flat skies, so a large bend is the surface following nebulae or')
    say('galaxies (the cost of a high tolerance).')
    say()
    for deg in (2, 4):
        say('degree %d' % deg)
        say('%-30s %-10s ' % ('image', 'case') + ''.join('%-14s' % ('tol %g' % t) for t in TOLERANCES))
        for im in manifest():
            a = load(im['file'], half=True)
            mask = star_mask(a)
            sky = float(np.median(sky_boxes(a, mask)))
            # (M31's sky is clipped to black: at least 0.02 linear)
            amp = max(sky, 0.02)
            h, w = a.shape[:2]
            grads = gradients(h, w)
            rows = {k: '%-30s %-10s ' % (im['file'][:30], k) for k in ['quadratic', 'dome', 'band', 'bend']}
            for t in TOLERANCES:
                props = {'degree': deg, 'tolerance': t}
                m0 = run('astro:background', dict(props, show='model'), a)
                lm = lum(m0)
                bend = 255 * (srgb(np.percentile(lm, 99)) - srgb(np.percentile(lm, 1)))
                for k, g in grads.items():
                    mx, mn = model_error(a, m0, g, amp, props)
                    rows[k] += '%5.1f / %4.2f  ' % (mx, mn)
                rows['bend'] += '%6.1f        ' % bend
            for k in rows:
                say(rows[k])
        say()


# ---- star mask -----------------------------------------------------------

def part_stars():
    say('Star Mask (defaults) on the real images (full size): regions found,')
    say('the share of the frame under the mask, and per megapixel. Overlay')
    say('crops (1:1, 800 x 600 from the centre and a corner) to look for')
    say('false detections by eye.')
    say()
    say('%-32s %-10s %-9s %-10s %-8s' % ('image', 'size', 'stars', 'per MP', 'masked'))
    for im in manifest():
        a = load(im['file'])
        h, w = a.shape[:2]
        m = star_mask(a)
        lab, n = ndimage.label(m > 0.5)
        say('%-32s %-10s %-9d %-10.0f %5.1f %%' % (im['file'], '%dx%d' % (w, h), n,
                                                  n / (w * h / 1e6), 100 * (m > 0.5).mean()))
        ov = run('astro:star-mask', {'show': 'overlay'}, a)
        base = os.path.splitext(im['file'])[0]
        cy, cx = h // 2 - 300, w // 2 - 400
        crops = [a[cy:cy + 600, cx:cx + 800], ov[cy:cy + 600, cx:cx + 800],
                 a[:600, :800], ov[:600, :800]]
        top = side_by_side(crops[0], crops[1])
        bottom = side_by_side(crops[2], crops[3])
        save_jpeg(np.concatenate([top, np.ones((8, top.shape[1], 4), np.float32), bottom]),
                  base + '-stars.jpg', 0)
    say()


# ---- reduce stars --------------------------------------------------------

def ring_dips(L, Lo, lab, order, objs):
    """per star: the radial profile (medians in 1 px rings about the
    light's centroid) of the input L and the result Lo, against the sky
    (the input's median 4 to 8 px beyond the region's farthest pixel), in
    noise sigmas (1.4826 MAD there); the lowest ring from the centre to 3
    px beyond the region, for the input and for the result"""
    h, w = L.shape
    dips_in, dips_out = [], []
    for k in order:
        sl = objs[k - 1]
        y0, y1 = max(0, sl[0].start - 12), min(h, sl[0].stop + 12)
        x0, x1 = max(0, sl[1].start - 12), min(w, sl[1].stop + 12)
        region = lab[y0:y1, x0:x1] == k
        if region.sum() > 2000:
            continue
        li, lo = L[y0:y1, x0:x1], Lo[y0:y1, x0:x1]
        yy, xx = np.mgrid[y0:y1, x0:x1]
        wgt = np.where(region, li, 0)
        cy, cx = (wgt * yy).sum() / wgt.sum(), (wgt * xx).sum() / wgt.sum()
        r = np.hypot(yy - cy, xx - cx)
        rmax = r[region].max()
        # the sky ring: no other star's pixels
        sky_px = (r >= rmax + 4) & (r <= rmax + 8) & (lab[y0:y1, x0:x1] == 0)
        if sky_px.sum() < 16:
            continue
        sky = np.median(li[sky_px])
        noise = 1.4826 * np.median(np.abs(li[sky_px] - sky)) + 1e-6
        others = (lab[y0:y1, x0:x1] != 0) & ~region
        lows = []
        for img in (li, lo):
            low = np.inf
            for rr in range(0, int(rmax + 3) + 1):
                ring = (r >= rr - 0.5) & (r < rr + 0.5) & ~others
                if ring.sum() >= (6 if rr > 1 else 1):
                    low = min(low, (np.median(img[ring]) - sky) / noise)
            lows.append(low)
        dips_in.append(lows[0])
        dips_out.append(lows[1])
    return np.array(dips_in), np.array(dips_out)


def part_reduce():
    say('Reduce Stars on the real images (full size): dark rings. For each of')
    say('the 300 star regions of the mask with the most light: the radial')
    say('profile (medians in 1 px rings about the star, the centre pixel too) against the sky 4 to')
    say('8 px beyond the region, in noise sigmas; "lowest" is the lowest ring')
    say('from the centre to 3 px beyond the region. A dark ring shows as a')
    say('clearly negative lowest ring. The input\'s own value is given too: these')
    say('are finished JPEGs, some with dark rings from their own processing.')
    say()
    say('%-32s %-22s %-6s %-22s %-22s %-8s %s' % ('image', 'settings', 'stars',
                                              'lowest: median input', 'median / worst result',
                                              '< -1', 'deeper'))
    settings = [('defaults', {}), ('erosion 100 %', {'method': 'erosion', 'amount': 100.0}),
                ('erosion 100 %, 2 it.', {'method': 'erosion', 'amount': 100.0, 'iterations': 2})]
    for im in manifest():
        a = load(im['file'])
        h, w = a.shape[:2]
        m = star_mask(a)
        lab, n = ndimage.label(m > 0)
        L = lum(a)
        flux = ndimage.sum(L, lab, np.arange(1, n + 1))
        order = np.argsort(-flux)[:300] + 1
        objs = ndimage.find_objects(lab)
        base = os.path.splitext(im['file'])[0]
        outs = {}
        for label, props in settings:
            o = run('astro:star-reduce', props, a)
            outs[label] = o
            di, do = ring_dips(L, lum(o), lab, order, objs)
            say('%-32s %-22s %-6d %-22.2f %6.2f / %-13.2f %-8s %.2f' % (
                im['file'], label, len(do), np.median(di), np.median(do), do.min(),
                '%d (%d)' % ((do < -1).sum(), (di < -1).sum()), (di - do).max()))
        # 4x crops about the brightest star whose region is under 40 px:
        # input | defaults | erosion 100 %; and 1:1 centre crops
        for k in order:
            sl = objs[k - 1]
            if (sl[0].stop - sl[0].start) < 40:
                break
        cy = (sl[0].start + sl[0].stop) // 2
        cx = (sl[1].start + sl[1].stop) // 2
        y0, x0 = min(max(0, cy - 60), h - 120), min(max(0, cx - 80), w - 160)
        zoom = lambda x: np.repeat(np.repeat(x[y0:y0 + 120, x0:x0 + 160], 4, 0), 4, 1)
        save_jpeg(side_by_side(side_by_side(zoom(a), zoom(outs['defaults'])),
                               zoom(outs['erosion 100 %'])), base + '-reduce-crop.jpg', 0)
        cut = lambda x: x[h // 2 - 300:h // 2 + 300, w // 2 - 400:w // 2 + 400]
        save_jpeg(side_by_side(cut(a), cut(outs['defaults'])), base + '-reduce.jpg', 0)
    say('("< -1": stars whose lowest ring in the result is below -1 sigma, and')
    say(' in brackets the same count for the input; "deeper": the most any')
    say(' star\'s lowest ring went below the input\'s own, in sigmas)')
    say()


# ---- stretch -------------------------------------------------------------

def encoded(a):
    """linear RGBA float (sRGB) to its encoded R'G'B' values"""
    h, w = a.shape[:2]
    return from_buffer(to_buffer(a), w, h, "R'G'B'A float")


def part_stretch():
    say('Astro Stretch on the real images (half size). Each finished image is')
    say('made dark and "linear" first by the inverse GHS on its encoded values')
    say('(standing in for a linear stack, whose numbers a file holds as they')
    say('are), then stretched forward with the same settings: with "encoded')
    say('values", as Siril would treat those numbers, and in linear light (the')
    say('default). Encoded must give the original back; linear differs (on')
    say('these numbers, held with the sRGB curve, it comes out darker: it')
    say('sees them as much smaller linear values, far below the symmetry')
    say('point). Sky: the median of the encoded luminance, 8 bit levels.')
    say()
    ghs = {'stretch': 4.0, 'local-intensity': 2.0, 'symmetry-point': 0.08,
           'protect-highlights': 0.95}
    say('GHS settings: %s' % ghs)
    say('%-32s %-16s %-10s %-10s %-12s %-12s' % ('image', 'max |back-orig|', 'sky orig',
                                                 'sky dark', 'sky encoded', 'sky linear'))
    sky = lambda x: 255 * float(np.median(x[..., :3] @ LUM))
    for im in manifest():
        a = load(im['file'], half=True)
        dark = run('astro:stretch', dict(ghs, inverse=True, encoded=True), a)
        back_enc = run('astro:stretch', dict(ghs, encoded=True), dark)
        back_lin = run('astro:stretch', ghs, dark)
        e, d, be, bl = encoded(a), encoded(dark), encoded(back_enc), encoded(back_lin)
        say('%-32s %-16.2e %-10.1f %-10.1f %-12.1f %-12.1f' % (
            im['file'], np.abs(be[..., :3] - e[..., :3]).max(), sky(e), sky(d), sky(be), sky(bl)))
        base = os.path.splitext(im['file'])[0]
        # the dark input | encoded | linear
        save_jpeg(side_by_side(side_by_side(dark, back_enc), back_lin), base + '-stretch.jpg', 2400)
    say()


def main():
    global log_file
    parts = {'background': part_background, 'tolerance': part_tolerance, 'stars': part_stars,
             'reduce': part_reduce, 'stretch': part_stretch}
    which = sys.argv[1:] or list(parts)
    missing = [im['file'] for im in manifest()
               if not os.path.exists(os.path.join(images_dir, im['file']))]
    if missing:
        sys.exit('missing images (run tests/images/fetch.py): ' + ', '.join(missing))
    for p in which:
        log_file = open(os.path.join(out_dir, p + '.txt'), 'w')
        parts[p]()
        log_file.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
