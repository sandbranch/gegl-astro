#!/bin/sh
# The operations on the real CC0 images of tests/images (fetched into
# tests/output/images by tests/images/fetch.py) and on synthetic scenes
# (tests/real.py): measurements and JPEGs in tests/output/real, not a
# pass or fail test. Uses the modules of tests/output/build-check (or
# $BUILD; tests/check.sh builds it), natively with a throwaway home, and
# lists the user's GIMP folders before and after (they must not change).
# Needs numpy and scipy (python3-numpy, python3-scipy).
#
#   tests/real.sh                       every part
#   tests/real.sh tolerance stretch     some parts (see tests/real.py)
#
# Copyright 2026 David
# SPDX-License-Identifier: GPL-3.0-or-later
here=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=SCRIPTDIR/common.sh
. "$here/common.sh"
build=${BUILD:-tests/output/build-check}
case $build in /*) ;; *) build=$top/$build ;; esac
mod="$out/real-modules"
modules="astro-background astro-stretch astro-star-mask astro-star-reduce astro-star-color"
for m in $modules; do
    [ -f "$build/$m.so" ] || { echo "no $build/$m.so: run tests/check.sh first" >&2; exit 2; }
done
python3 -I "$here/images/fetch.py" || exit 2
rm -rf "$mod"
mkdir -p "$mod"
# GEGL loads every file in a module folder: only the modules
for m in $modules; do cp "$build/$m.so" "$mod/"; done

snapshot_begin real
GEGL_PATH="$mod:$(pkg-config --variable=pluginsdir gegl-0.4)"
export GEGL_PATH
pass_env="GEGL_PATH"
run_prefix="--timeout=3600" in_gimp python3 -I "$here/real.py" "$@"
status=$?
snapshot_end real || status=1
exit $status
