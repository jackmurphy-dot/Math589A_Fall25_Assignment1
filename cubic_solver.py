# cubic_solver.py
import math
import cmath

def solve_cubic(a, b, c, d):
    tol = 1e-14
    if abs(a) < tol:
        from quartic_solver import solve_quadratic
        return solve_quadratic(b, c, d)

    # Normalize
    A = b/a
    B = c/a
    C = d/a

    # Depressed cubic: x = y - A/3
    p = B - A*A/3
    q = 2*A**3/27 - A*B/3 + C

    roots = []

    # Casus irreducibilis (3 real roots)
    discriminant = (q/2)**2 + (p/3)**3
    if abs(discriminant) < tol:
        discriminant = 0

    if discriminant > 0:
        # One real, two complex
        u = cmath.exp(cmath.log(-q/2 + cmath.sqrt(discriminant))/3)
        v = cmath.exp(cmath.log(-q/2 - cmath.sqrt(discriminant))/3)
        y1 = u + v
        roots.append(y1 - A/3)
        omega = complex(-0.5, math.sqrt(3)/2)
        roots.append(omega*u + omega.conjugate()*v - A/3)
        roots.append(omega.conjugate()*u + omega*v - A/3)
    elif discriminant == 0:
        u = cmath.exp(cmath.log(-q/2)/3)
        roots.append(2*u - A/3)
        roots.append(-u - A/3)
        roots.append(-u - A/3)
    else:
        r = 2*math.sqrt(-p/3)
        phi = math.acos(-q/(2*((-p/3)**1.5)))
        roots.append(r*math.cos(phi/3) - A/3)
        roots.append(r*math.cos((phi + 2*math.pi)/3) - A/3)
        roots.append(r*math.cos((phi + 4*math.pi)/3) - A/3)

    return roots
