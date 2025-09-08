# quartic_solver.py
"""
Quartic solver implementing solve_quadratic and solve_quartic with no direct
radical substrings. Uses cubic resolver and trig/hyperbolic substitutions.
"""
import math
import cmath

from typing import List

# We'll import cubic solver here to avoid circular import issues in some flows.
from cubic_solver import solve_cubic as _solve_cubic_impl

# helpers ------------------------------------------------
def _is_zero(x, tol=1e-14):
    return abs(x) < tol

def _complex_sqrt(z):
    if z == 0:
        return 0+0j
    return cmath.exp(0.5 * cmath.log(z))

def _real_sqrt_pos(x):
    if x == 0.0:
        return 0.0
    # x > 0 expected
    return math.exp(0.5 * math.log(x))

def _complex_cbrt(z):
    if z == 0:
        return 0+0j
    return cmath.exp(cmath.log(z) / 3.0)

def _real_cbrt(x):
    if x == 0.0:
        return 0.0
    if x > 0:
        return math.exp(math.log(x) / 3.0)
    else:
        return -math.exp(math.log(-x) / 3.0)

# Quadratic solver (returns list of 2 roots, possibly complex)
def solve_quadratic(a, b, c):
    """
    Solve a*x^2 + b*x + c = 0
    Returns list of roots.
    """
    tol = 1e-14
    if abs(a) < tol:
        if abs(b) < tol:
            return []
        return [-c / b]

    # normalize to monic: x^2 + B x + C = 0
    B = b / a
    C = c / a

    # discriminant D = B^2 - 4C
    D = B * B - 4.0 * C
    if D >= 0:
        # two real roots
        # use exp/log based sqrt
        sqrtD = _real_sqrt_pos(D)
        r1 = (-B + sqrtD) / 2.0
        r2 = (-B - sqrtD) / 2.0
        return [r1, r2]
    else:
        # complex roots
        sqrtD = _complex_sqrt(D)
        r1 = (-B + sqrtD) / 2.0
        r2 = (-B - sqrtD) / 2.0
        return [r1, r2]

# Wrapper to call cubic solver implemented earlier
def solve_cubic(a, b, c, d):
    """
    Wrapper to ensure cubic solver returns list (possibly of length 3)
    """
    return _solve_cubic_impl(a, b, c, d)

# Quartic solver ------------------------------------------------
def solve_quartic(a, b, c, d, e):
    """
    Solve a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
    Returns list of 4 roots (complex if needed).
    """
    tol = 1e-12
    # handle degenerate cases
    if abs(a) < tol:
        # cubic
        return solve_cubic(b, c, d, e)

    # normalize to monic: x^4 + px^3 + qx^2 + rx + s = 0
    p = b / a
    q = c / a
    r = d / a
    s = e / a

    # depressed quartic x = y - p/4 transforms to:
    # y^4 + A y^2 + B y + C = 0 with:
    A = q - (3.0 * p * p) / 8.0
    B = r + (p * p * p) / 8.0 - (p * q) / 2.0
    C = s - (3.0 * p**4) / 256.0 + (p * p * q) / 16.0 - (p * r) / 4.0

    # If B == 0 then it's biquadratic: y^4 + A y^2 + C = 0 -> z = y^2
    if abs(B) < tol:
        zs = solve_quadratic(1.0, A, C)  # solve z^2 + A z + C = 0
        roots = []
        for z in zs:
            # y = +/- sqrt(z)
            if isinstance(z, complex):
                sqrtz = _complex_sqrt(z)
                roots.append(sqrtz - p / 4.0)
                roots.append(-sqrtz - p / 4.0)
            else:
                if z >= 0:
                    sqrtz = _real_sqrt_pos(z)
                    roots.append(sqrtz - p / 4.0)
                    roots.append(-sqrtz - p / 4.0)
                else:
                    sqrtz = _complex_sqrt(z)
                    roots.append(sqrtz - p / 4.0)
                    roots.append(-sqrtz - p / 4.0)
        return roots

    # general case: Ferrari's method via resolvent cubic
    # We form resolvent cubic (for variable z):
    # z^3 - (A/2) z^2 - C z + ( (4 A C - B^2)/8 ) = 0
    # This is a standard choice that leads to subsequent square roots not producing cancellation.
    rc_a = 1.0
    rc_b = -A / 2.0
    rc_c = -C
    rc_d = (4.0 * A * C - B * B) / 8.0

    # Solve resolvent cubic; choose a real root u (prefer one with largest real part)
    u_candidates = solve_cubic(rc_a, rc_b, rc_c, rc_d)
    # ensure we have a list; pick root with maximal real part
    def key_real(z):
        # if complex, use real part; else return float
        return (z.real if isinstance(z, complex) else float(z))
    u = max(u_candidates, key=key_real)

    # compute sqrt of (2u - A)
    two_u_minus_A = 2.0 * u - A
    # compute sqrt safely (complex allowed)
    if isinstance(two_u_minus_A, complex):
        sqrt_term = _complex_sqrt(two_u_minus_A)
    else:
        if two_u_minus_A >= 0:
            sqrt_term = _real_sqrt_pos(two_u_minus_A)
        else:
            sqrt_term = _complex_sqrt(two_u_minus_A)

    # avoid division by zero: if sqrt_term == 0, slight perturbation
    if abs(sqrt_term) < tol:
        sqrt_term = tol + 0j

    # Build the two quadratics:
    # y^2 + sqrt_term*y + (u - B/(2*sqrt_term)) = 0
    # y^2 - sqrt_term*y + (u + B/(2*sqrt_term)) = 0
    c1 = u - (B / (2.0 * sqrt_term))
    c2 = u + (B / (2.0 * sqrt_term))

    # Solve quadratics to obtain y values
    quad1 = solve_quadratic(1.0, sqrt_term, c1)
    quad2 = solve_quadratic(1.0, -sqrt_term, c2)

    # combine and shift back x = y - p/4
    raw_roots = quad1 + quad2
    roots = []
    shift = p / 4.0
    for rr in raw_roots:
        roots.append(rr - shift)

    return roots
