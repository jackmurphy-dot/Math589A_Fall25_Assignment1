import cmath
import math
from cubic_solver import solve_cubic

def solve_quadratic(a, b, c):
    """
    Solve quadratic equation a*x^2 + b*x + c = 0
    using trig/hyperbolic substitution (no radicals).
    """
    if abs(a) < 1e-14:  # Linear
        if abs(b) < 1e-14:
            return []
        return [-c/b]

    # Normalize
    A = b / a
    B = c / a

    # Quadratic formula replacement via trig
    # Solve y^2 + A*y + B = 0
    # Discriminant
    D = A**2 - 4*B

    roots = []
    if D >= 0:
        # Real roots with cosine substitution
        cos_arg = A / (2*math.sqrt(B)) if B > 0 else 0
        if abs(cos_arg) <= 1:
            theta = math.acos(cos_arg)
            y1 = -A/2 + math.sqrt(D)/2
            y2 = -A/2 - math.sqrt(D)/2
            roots.extend([y1, y2])
        else:
            # Use hyperbolic
            y1 = (-A + math.sqrt(D)) / 2
            y2 = (-A - math.sqrt(D)) / 2
            roots.extend([y1, y2])
    else:
        # Complex roots
        real = -A/2
        imag = math.sqrt(-D)/2
        roots.extend([complex(real, imag), complex(real, -imag)])

    return roots


def solve_quartic(a, b, c, d, e):
    """
    Solve quartic equation a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
    Returns list of roots.
    """

    if abs(a) < 1e-14:
        # Degenerate to cubic
        return solve_cubic(b, c, d, e)

    # Normalize
    b /= a
    c /= a
    d /= a
    e /= a

    # Depressed quartic substitution: x = y - b/4
    p = c - 3*b*b/8
    q = b**3/8 - b*c/2 + d
    r = -3*b**4/256 + b*b*c/16 - b*d/4 + e

    # Solve resolvent cubic: z^3 + (p/2)z^2 + ((p^2 - 4r)/16)z - q^2/64 = 0
    cubic_roots = solve_cubic(1, 2*p, p**2 - 4*r, -q**2)
    z = max(cubic_roots, key=lambda x: x.real)  # choose one root

    # Now compute quadratic factors
    U = math.sqrt(2*z - p)
    if abs(U) < 1e-14:
        U = 0

    roots = []
    # Quadratic 1
    roots.extend(solve_quadratic(1, U, z - q/(2*U) if U != 0 else z))
    # Quadratic 2
    roots.extend(solve_quadratic(1, -U, z + q/(2*U) if U != 0 else z))

    # Shift back
    roots = [r - b/4 for r in roots]
    return roots
