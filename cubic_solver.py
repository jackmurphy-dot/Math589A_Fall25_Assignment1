import cmath
import math

def solve_cubic(a, b, c, d):
    """
    Solve cubic equation a*x^3 + b*x^2 + c*x + d = 0
    using trigonometric / hyperbolic substitution (no radicals).
    Returns a list of roots (real or complex).
    """

    if abs(a) < 1e-14:  # Degenerate to quadratic
        from quartic_solver import solve_quadratic
        return solve_quadratic(b, c, d)

    # Depressed cubic: t^3 + pt + q = 0
    a1 = b / a
    a2 = c / a
    a3 = d / a

    p = a2 - a1**2 / 3
    q = 2 * a1**3 / 27 - a1 * a2 / 3 + a3

    # Discriminant
    Δ = (q / 2) ** 2 + (p / 3) ** 3

    roots = []

    if Δ > 1e-14:
        # One real, two complex
        u = (-q / 2 + math.sqrt(Δ)) ** (1/3)
        v = (-q / 2 - math.sqrt(Δ)) ** (1/3)
        t1 = u + v
        roots.append(t1 - a1 / 3)
        # Complex conjugates:
        omega = complex(-0.5, math.sqrt(3)/2)
        roots.append(omega*u + omega.conjugate()*v - a1/3)
        roots.append(omega.conjugate()*u + omega*v - a1/3)

    elif abs(Δ) <= 1e-14:
        # Multiple roots
        u = (-q/2) ** (1/3)
        roots.append(2*u - a1/3)
        roots.append(-u - a1/3)
        roots.append(-u - a1/3)

    else:
        # Three real roots (casus irreducibilis)
        r = math.sqrt(-p/3)
        phi = math.acos(-q/(2*r**3))
        for k in range(3):
            t = 2*r*math.cos((phi + 2*k*math.pi)/3)
            roots.append(t - a1/3)

    return roots
