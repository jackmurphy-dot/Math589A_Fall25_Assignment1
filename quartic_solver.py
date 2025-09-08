# quartic_solver.py
import math
import cmath
from cubic_solver import solve_cubic as _solve_cubic_impl

def _complex_sqrt(z):
    if z == 0:
        return 0+0j
    return cmath.exp(0.5*cmath.log(z))

def solve_quadratic(a, b, c):
    tol = 1e-14
    if abs(a) < tol:
        if abs(b) < tol:
            return []
        return [-c/b]
    B = b / a
    C = c / a
    D = B*B - 4*C
    if isinstance(D, complex) or D < 0:
        sqrtD = _complex_sqrt(D)
        r1 = (-B + sqrtD)/2
        r2 = (-B - sqrtD)/2
    else:
        sqrtD = math.exp(0.5*math.log(D))
        r1 = (-B + sqrtD)/2
        r2 = (-B - sqrtD)/2
    return [r1, r2]

def solve_cubic(a, b, c, d):
    return _solve_cubic_impl(a, b, c, d)

def solve_quartic(a, b, c, d, e):
    tol = 1e-12
    if abs(a) < tol:
        return solve_cubic(b, c, d, e)

    # normalize
    p = b/a
    q = c/a
    r = d/a
    s = e/a

    # depressed quartic
    A = q - 3*p*p/8
    B = r + p*p*p/8 - p*q/2
    C = s - 3*p**4/256 + p*p*q/16 - p*r/4

    if abs(B) < tol:
        # biquadratic
        zs = solve_quadratic(1, A, C)
        roots = []
        for z in zs:
            if isinstance(z, complex) or z < 0:
                sqrtz = _complex_sqrt(z)
            else:
                sqrtz = math.exp(0.5*math.log(z))
            roots.append(sqrtz - p/4)
            roots.append(-sqrtz - p/4)
        return roots

    # general case: Ferrari's method
    rc_a = 1
    rc_b = -A/2
    rc_c = -C
    rc_d = (4*A*C - B*B)/8

    u_candidates = solve_cubic(rc_a, rc_b, rc_c, rc_d)
    # pick largest real part
    def key_real(z):
        return z.real if isinstance(z, complex) else float(z)
    u = max(u_candidates, key=key_real)

    two_u_minus_A = 2*u - A
    if isinstance(two_u_minus_A, complex) or two_u_minus_A < 0:
        sqrt_term = _complex_sqrt(two_u_minus_A)
    else:
        sqrt_term = math.exp(0.5*math.log(two_u_minus_A))

    if abs(sqrt_term) < tol:
        sqrt_term += tol + 0j

    c1 = u - (B/(2*sqrt_term))
    c2 = u + (B/(2*sqrt_term))

    quad1 = solve_quadratic(1, sqrt_term, c1)
    quad2 = solve_quadratic(1, -sqrt_term, c2)

    roots = []
    for rr in quad1 + quad2:
        roots.append(rr - p/4)

    return roots
