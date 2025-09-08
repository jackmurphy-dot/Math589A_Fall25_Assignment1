# quartic_solver.py
import math
import cmath
from cubic_solver import solve_cubic

def solve_quadratic(a, b, c):
    tol = 1e-14
    if abs(a) < tol:
        if abs(b) < tol:
            return []
        return [-c/b]
    # Depressed quadratic: x = y - b/(2a)
    p = -b/(2*a)
    q = c/a - (b*b)/(4*a*a)
    discriminant = -q
    roots = []
    if discriminant >= 0:
        sqrt_disc = math.sqrt(discriminant)
        roots.append(p + sqrt_disc)
        roots.append(p - sqrt_disc)
    else:
        sqrt_disc = cmath.exp(0.5*cmath.log(-discriminant))*1j
        roots.append(p + sqrt_disc)
        roots.append(p - sqrt_disc)
    return roots

def solve_quartic(a, b, c, d, e):
    tol = 1e-12
    if abs(a) < tol:
        return solve_cubic(b, c, d, e)

    # Normalize
    p = b/a
    q = c/a
    r = d/a
    s = e/a

    # Depressed quartic x = y - p/4
    A = -3*p*p/8 + q
    B = p*p*p/8 - p*q/2 + r
    C = -3*p**4/256 + p*p*q/16 - p*r/4 + s

    roots = []

    if abs(B) < tol:
        # Biquadratic
        zs = solve_quadratic(1, A, C)
        for z in zs:
            if isinstance(z, complex) or z < 0:
                sqrt_z = cmath.exp(0.5*cmath.log(z))
            else:
                sqrt_z = math.sqrt(z)
            roots.append(sqrt_z - p/4)
            roots.append(-sqrt_z - p/4)
        return roots

    # General Ferrari method
    # Resolvent cubic
    rc_a = 1
    rc_b = -A/2
    rc_c = -C
    rc_d = (4*A*C - B*B)/8

    u_candidates = solve_cubic(rc_a, rc_b, rc_c, rc_d)
    # Pick real solution if possible
    u = max(u_candidates, key=lambda x: x.real if isinstance(x, complex) else x)

    two_u_minus_A = 2*u - A
    if isinstance(two_u_minus_A, complex) or two_u_minus_A < 0:
        sqrt_term = cmath.exp(0.5*cmath.log(two_u_minus_A))
    else:
        sqrt_term = math.sqrt(two_u_minus_A)

    c1 = u - B/(2*sqrt_term)
    c2 = u + B/(2*sqrt_term)

    quad1 = solve_quadratic(1, sqrt_term, c1)
    quad2 = solve_quadratic(1, -sqrt_term, c2)

    for rr in quad1 + quad2:
        roots.append(rr - p/4)

    return roots
