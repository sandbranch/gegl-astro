#!/bin/sh
# Each filter's dialog in GIMP on a Broadway display, driven with a
# headless Chrome: opens the filter by its title from Search Actions (/),
# checks that a dialog opened, moves it to the top right and leaves a
# screenshot in tests/output/gui/<name>-dialog.png (copy one to docs/ by
# hand when it should go into the README). Run tests/gimp-check.sh first
# (it puts the modules in place). Lists the user's folders before and
# after and fails if anything changed.
#
#   tests/gui/dialogs.sh                        all
#   ASTRO_DIALOGS="stretch" tests/gui/dialogs.sh  some
#     (background stretch star-mask star-reduce)
#
# Prints PASS or FAIL per step and exits non-zero if one failed.
# (After gegl-equalizers' tests/gui/shots.sh.)
#
# Copyright 2026 David
# SPDX-License-Identifier: GPL-3.0-or-later
here=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=SCRIPTDIR/common.sh
. "$here/common.sh"
status=0
pass () { echo "PASS  gui_$1"; }
fail () { echo "FAIL  gui_$1"; status=1; }

# the dialog: the canvas made last that is neither the whole page nor
# small (Broadway draws each window as a canvas with its shadow, the
# title bar's middle 36 pixels below the canvas's top)
dialog_rect () {
    $cdp "$view" "eval:(() => { const c = [...document.querySelectorAll('canvas')]
      .map(e => e.getBoundingClientRect())
      .filter(r => r.width < 1200 && r.height > 150 && r.height < 780).pop();
      return c ? [c.left, c.top, c.width, c.height].map(Math.round).join(',') : '' })()" \
      2>/dev/null | tail -1
}

# moves the dialog to the top right corner of the page
move_dialog () {
    r=$(dialog_rect)
    [ -n "$r" ] || return 1
    l=${r%%,*}; r=${r#*,}; t=${r%%,*}; r=${r#*,}; w=${r%%,*}
    vw=$((w - 52))
    vx=$((1278 - vw))
    vy=120
    gx=$((l + w / 2)); gy=$((t + 36))
    $cdp "$view" down:$gx,$gy move:$(( (gx + vx + vw / 2) / 2 )),$(( (gy + vy + 10) / 2 )) \
      move:$((vx + vw / 2)),$((vy + 10)) up:$((vx + vw / 2 + 1)),$((vy + 10)) wait:1200 >/dev/null
}

# on Broadway the search field needs a click before it takes the keys
SEARCH_FIELD=360,176

# search_filter "<title>" <name>: opens it from Search Actions
search_filter () {
    $cdp "$view" key:/ wait:1500 click:$SEARCH_FIELD wait:500 text:"$1" wait:1500 \
      shot:"$gui_out/$2-search.png" key:Enter wait:6000 >/dev/null
    move_dialog || { fail "$2: no dialog for \"$1\""; return 1; }
    pass "$2: the dialog of \"$1\" opened"
}

title_of () {
    case $1 in
        background) echo "Remove Sky Background" ;;
        stretch) echo "Astro Stretch" ;;
        star-mask) echo "Star Mask" ;;
        star-reduce) echo "Reduce Stars" ;;
    esac
}

session () {
    name=$1
    title=$(title_of "$name")
    [ -n "$title" ] || { fail "$name: no such filter"; return; }
    ASTRO_GUI_NAME=$name
    export ASTRO_GUI_NAME
    tries=0
    until start_gimp "$name"; do
        stop_all
        tries=$((tries + 1))
        [ $tries -lt 3 ] || { fail "$name: GIMP did not start (log: $gui_out/gimp-$name.log)"; return; }
    done
    if search_filter "$title" "$name"; then
        $cdp "$view" wait:6000 shot:"$gui_out/$name-dialog.png" >/dev/null &&
          pass "$name: screenshot $gui_out/$name-dialog.png"
    fi
    stop_all
}

snapshot_begin gui
for s in ${ASTRO_DIALOGS:-background stretch star-mask star-reduce}; do
    session "$s"
done
snapshot_end gui || status=1
[ $status = 0 ] && echo "gui: all passed" || echo "gui: FAILED (screenshots and logs in $gui_out)"
exit $status
