#!/usr/bin/env python3
# Downloads the real astronomical images listed in manifest.json from
# Wikimedia Commons into tests/output/images, as the renditions Commons
# serves for a 2560 pixel request (3840 px wide on 2026-11-18; those
# marked "full_size": the original, which is smaller),
# and checks each against the SHA-256 in the manifest. All are CC0: the
# licence was read on each file's Commons page (manifest.json names the
# page, author, licence and the licence templates found there, with the
# Flickr licence review). They stay out of git; this script makes the
# same set anywhere. tests/real.sh runs the filters on them.
#
#   tests/images/fetch.py           fetch what is missing, check all
#   tests/images/fetch.py --update  write the SHA-256 of files that have
#                                   none in the manifest yet
#
# (After gegl-photo-restore/tests/images/fetch.py: Commons' media servers
# answer Python's own HTTP client with 429 and serve curl, and the 2560
# pixel renditions are more reliable than large originals.)
#
# Copyright 2026 David
# SPDX-License-Identifier: GPL-3.0-or-later
import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.parse

here = os.path.dirname(os.path.abspath(__file__))
manifest_path = os.path.join(here, 'manifest.json')
dest = os.path.join(os.path.dirname(here), 'output', 'images')
# a generic client name, no contact details
UA = 'astro-test-fetch/0.1 (test image fetcher)'


def get(url, timeout=300):
    """the body of url, through curl; waits and tries again when asked to
    slow down"""
    for attempt in range(10):
        r = subprocess.run(['curl', '-sS', '-L', '--max-time', str(timeout), '-A', UA,
                            '-w', '\n%{http_code}', url], capture_output=True, check=True)
        body, _, code = r.stdout.rpartition(b'\n')
        if code == b'200':
            return body
        if code not in (b'429', b'503') or attempt == 9:
            sys.exit('%s: HTTP %s' % (url, code.decode()))
        print('  Commons asks to wait, %d s' % (60 * (attempt + 1)))
        time.sleep(60 * (attempt + 1))


def file_url(image):
    """the original (with "full_size") or Commons' 2560 pixel wide rendition"""
    params = dict(action='query', format='json', titles=image['title'], prop='imageinfo', iiprop='url')
    if not image.get('full_size'):
        params['iiurlwidth'] = 2560
    data = json.loads(get('https://commons.wikimedia.org/w/api.php?' + urllib.parse.urlencode(params), 60))
    info = next(iter(data['query']['pages'].values()))['imageinfo'][0]
    return info.get('thumburl', info['url'])


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    update = '--update' in sys.argv
    images = json.load(open(manifest_path))
    os.makedirs(dest, exist_ok=True)
    bad = 0
    for image in images:
        path = os.path.join(dest, image['file'])
        if not os.path.exists(path):
            data = get(file_url(image))
            with open(path + '.part', 'wb') as f:
                f.write(data)
            os.rename(path + '.part', path)
            print('%-32s %6d kB  %s (%s)' % (image['file'], len(data) // 1024,
                                             image['author'][:40], image['license']))
            time.sleep(3)
        digest = sha256(path)
        if not image.get('sha256'):
            if update:
                image['sha256'] = digest
            else:
                print('%-32s no SHA-256 in the manifest (run with --update)' % image['file'])
                bad += 1
        elif digest != image['sha256']:
            print('%-32s SHA-256 differs from the manifest: %s' % (image['file'], digest))
            bad += 1
    if update:
        with open(manifest_path, 'w') as f:
            json.dump(images, f, indent=1, ensure_ascii=False)
            f.write('\n')
    print('%d images in %s, %s' % (len(images), os.path.relpath(dest, os.getcwd()),
                                   'all checked' if not bad else '%d not right' % bad))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
