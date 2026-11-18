#!/usr/bin/env python3
# Checks astro:stretch against the published equations, written anew here
# without looking at the operation's code: the Generalised Hyperbolic
# Stretch as given in the GHS documentation (ghsastro.co.uk, Dave Payne and
# Mike Cranfield; docs/research.md 2.6) and the arcsinh stretch of Lupton
# et al. (2004). tests/check.c writes the operation's results on ramps
# (stretch-op.txt in the build folder); this script computes the same
# values in double precision and prints PASS or FAIL per case and
# "N passed, M failed". The limit is 1e-6 (the operation's output is a
# 32 bit float).
#
#   tests/stretch-reference.py <build folder>/stretch-op.txt
#
# The inverse transform is checked against a bisection of the forward
# equations, not against inverse formulas.
#
# Copyright 2026 David
# SPDX-License-Identifier: GPL-3.0-or-later
import math
import sys

LIMIT = 1e-6


def base(b, D, x):
    """the base transform T(x) for local intensity b (x >= 0)"""
    if b == -1:
        return math.log(1 + D * x)
    if b < 0:
        return (1 - (1 - b * D * x) ** ((b + 1) / b)) / (D * (b + 1))
    if b == 0:
        return 1 - math.exp(-D * x)
    return 1 - (1 + b * D * x) ** (-1 / b)


def slope(b, D, x):
    """T'(x), by hand from the formulas above"""
    if b == -1:
        return D / (1 + D * x)
    if b < 0:
        return (1 - b * D * x) ** (1 / b)
    if b == 0:
        return D * math.exp(-D * x)
    return D * (1 + b * D * x) ** (-(1 + b) / b)


def ghs(x, b, stretch, lp, sp, hp):
    """the whole transform: mirrored about SP, straight below LP and above
    HP, normalised so 0 and 1 stay"""
    D = math.expm1(stretch)
    if D == 0:
        return x
    # the operation keeps LP <= SP <= HP
    lp = min(lp, sp)
    hp = max(hp, sp)

    def t1(v):
        if v < lp:
            return -base(b, D, sp - lp) + slope(b, D, sp - lp) * (v - lp)
        if v < sp:
            return -base(b, D, sp - v)
        if v < hp:
            return base(b, D, v - sp)
        return base(b, D, hp - sp) + slope(b, D, hp - sp) * (v - hp)

    lo, hi = t1(0.0), t1(1.0)
    return (t1(x) - lo) / (hi - lo)


def ghs_inverse(y, *args):
    """the x in 0 to 1 with ghs(x) = y, by bisection"""
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if ghs(mid, *args) < y:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def black(x, bp):
    return min(max((x - bp) / (1 - bp), 0.0), 1.0)


def main():
    cases = []
    with open(sys.argv[1]) as fh:
        for line in fh:
            w = line.split()
            if w[0] == 'case':
                cases.append((w[1:], []))
            else:
                cases[-1][1].append((float(w[0]), float(w[1])))
    passed = failed = 0
    for head, pairs in cases:
        label, kind = head[0], head[1]
        if kind == 'arcsinh':
            beta, bp = float(head[2]), float(head[3])

            def want(x):
                v = black(x, bp)
                return v if beta == 0 else math.asinh(beta * v) / math.asinh(beta)
            desc = 'arcsinh beta %g, black point %g' % (beta, bp)
        else:
            b, stretch, lp, sp, hp, bp = map(float, head[2:8])
            inverse = head[8] == '1'

            def want(x):
                v = black(x, bp)
                if inverse:
                    return ghs_inverse(v, b, stretch, lp, sp, hp)
                return ghs(v, b, stretch, lp, sp, hp)
            desc = 'GHS%s b %g, ln(D+1) %g, LP %g, SP %g, HP %g, black point %g' % (
                ' inverse' if inverse else '', b, stretch, lp, sp, hp, bp)
        worst, at = 0.0, 0.0
        for x, y in pairs:
            d = abs(y - want(x))
            if d > worst:
                worst, at = d, x
        ok = worst <= LIMIT
        passed += ok
        failed += not ok
        print('%s  stretch_reference_case_%s: %s: %d values, largest difference %.2g '
              '(at %.6g)' % ('PASS' if ok else 'FAIL', label, desc, len(pairs), worst, at))
    print('reference: %d passed, %d failed' % (passed, failed))
    return 1 if failed or not cases else 0


if __name__ == '__main__':
    sys.exit(main())
