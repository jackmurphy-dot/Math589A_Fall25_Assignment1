# cubic_solver.py
import cmath
import math

def _complex_cbrt(z):
    """Compute complex cube root using exp/log"""
    if z == 0:
        return 0+0j
    return cmath.exp(cmath.log(z)/3)

def solve_cubic(a, b, c, d):
    tol = 1e-14
    if abs(a) < tol:
        # Degenerate cubic => quadratic
        from quartic_solver import solve_quadratic
        return solve_quadratic(b, c, d)

    # Normalize to monic cubic: x^3 + A x^2 + B x + C = 0
    A = b / a
    B = c / a
    C = d / a

    # Depressed cubic: x = y - A/3
    A3 = A / 3
    p = B - A*A3
    q = 2*A3*A*A3 - A3*B + C

    roots = []

    discriminant = (q/2)**2 + (p/3)**3

    if abs(discriminant) < tol:
        discriminant = 0

    if discriminant > 0:
        # One real, two complex
        sqrt_disc = cmath.exp(0.5*cmath.log(discriminant))
        u = _complex_cbrt(-q/2 + sqrt_disc)
        v = _complex_cbrt(-q/2 - sqrt_disc)
        y1 = u + v
        roots.append(y1 - A3)
        omega = complex(-0.5, math.sqrt(3)/2)
        roots.append(omega*u + omega.conjugate()*v - A3)
        roots.append(omega.conjugate()*u + omega*v - A3)
    elif discriminant == 0:
        u = _complex_cbrt(-q/2)
        roots.append(2*u - A3)
        roots.append(-u - A3)
        roots.append(-u - A3)
    else:
        # Three real roots
        r = 2*math.sqrt(-p/3)
        phi = math.acos(-q/(2*((-p/3)**1.5)))
        roots.append(r*math.cos(phi/3) - A3)
        roots.append(r*math.cos((phi+2*math.pi)/3) - A3)
        roots.append(r*math.cos((phi+4*math.pi)/3) - A3)

    return roots
