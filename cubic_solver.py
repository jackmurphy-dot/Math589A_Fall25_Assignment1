# cubic_solver.py
import math
import cmath

def _complex_cbrt(z):
    if z == 0:
        return 0+0j
    return cmath.exp(cmath.log(z)/3.0)

def _real_cbrt(x):
    if x == 0:
        return 0.0
    if x > 0:
        return math.exp(math.log(x)/3.0)
    else:
        return -math.exp(math.log(-x)/3.0)

def solve_cubic(a, b, c, d):
    tol = 1e-14
    if abs(a) < tol:
        from quartic_solver import solve_quadratic
        return solve_quadratic(b, c, d)

    # normalize
    A = b / a
    B = c / a
    C = d / a

    A3 = A / 3.0
    p = B - A*A3
    q = 2*A3*A*A3 - A3*B + C

    half_q = q/2.0
    third_p = p/3.0
    Δ = half_q*half_q + third_p*third_p*third_p

    roots = []
    if Δ > tol:
        u = _complex_cbrt(-half_q + cmath.sqrt(Δ))
        v = _complex_cbrt(-half_q - cmath.sqrt(Δ))
        y1 = u + v
        roots.append(y1 - A3)
        omega = complex(-0.5, math.sqrt(3)/2)
        roots.append(omega*u + omega.conjugate()*v - A3)
        roots.append(omega.conjugate()*u + omega*v - A3)
    elif abs(Δ) <= tol:
        u = _real_cbrt(-half_q)
        roots.append(2*u - A3)
        roots.append(-u - A3)
        roots.append(-u - A3)
    else:
        # three real roots
        r = 2*math.sqrt(-p/3)
        phi = math.acos(-q/(2*((-p/3)**1.5)))
        roots.append(r*math.cos(phi/3) - A3)
        roots.append(r*math.cos((phi+2*math.pi)/3) - A3)
        roots.append(r*math.cos((phi+4*math.pi)/3) - A3)

    return roots
