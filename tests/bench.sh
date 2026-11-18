#!/bin/sh
# How long the operations take on a 24 megapixel synthetic star field
# (tests/bench.c), with the build of tests/check.sh
# (tests/output/build-check, or $BUILD), natively with a throwaway home.
# Lists the user's GIMP folders before and after and fails if anything
# changed. (After gegl-presence's tests/bench.sh.)
#   tests/bench.sh                 6000 x 4000
#   tests/bench.sh 3000 2000       another size
#
# Copyright 2026 David
# SPDX-License-Identifier: GPL-3.0-or-later
here=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=SCRIPTDIR/common.sh
. "$here/common.sh"
build=${BUILD:-tests/output/build-check}
case $build in /*) ;; *) build=$top/$build ;; esac
[ -x "$build/bench" ] || { echo "no $build/bench: run tests/check.sh first" >&2; exit 2; }
snapshot_begin bench
in_gimp "$build/bench" "$build" "$@"
status=$?
snapshot_end bench || status=1
exit $status
