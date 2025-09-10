# Radical-free cubic solver using trig identities and exp/log roots.

import cmath
from math import pi

_EPS = 1e-12

def _half_power(z): return cmath.exp(0.5 * cmath.log(z))
def _third_power(z): return cmath.exp(cmath.log(z) / 3)

def _solve_quadratic(a, b, c):
    if abs(a) < _EPS: return ([] if abs(b) < _EPS else [-c/b])
    disc = b*b - 4*a*c
    s = _half_power(disc)
    return [(-b+s)/(2*a), (-b-s)/(2*a)]

def solve_cubic(a, b, c, d):
    if abs(a) < _EPS:
        return [complex(r) for r in _solve_quadratic(b, c, d)]

    A, B, C = b/a, c/a, d/a
    p = B - A*A/3
    q = 2*A**3/27 - A*B/3 + C

    if abs(p) < _EPS:  # t^3 + q = 0
        rho = _third_power(-q)
        w1, w2 = cmath.exp(2j*pi/3), cmath.exp(4j*pi/3)
        ts = [rho, rho*w1, rho*w2]
    else:
        s = _half_power(-p/3)
        if abs(s) < _EPS:
            rho = _third_power(-q)
            w1, w2 = cmath.exp(2j*pi/3), cmath.exp(4j*pi/3)
            ts = [rho, rho*w1, rho*w2]
        else:
            phi = cmath.acos((-q/2)/(s**3))
            r = 2*s
            ts = [r*cmath.cos((phi+2*k*pi)/3) for k in range(3)]

    return [t - A/3 for t in ts]
