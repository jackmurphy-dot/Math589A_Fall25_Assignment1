# cubic_solver.py
import math
import cmath
from typing import List

def _cbrt(x: complex) -> complex:
    """Real-valued cube root when x is real; complex otherwise."""
    if isinstance(x, complex):
        return x ** (1/3)
    # handle negative real properly
    if x >= 0:
        return x ** (1/3)
    else:
        return -((-x) ** (1/3))

def solve_quadratic(a: complex, b: complex, c: complex) -> List[complex]:
    """Return list of two roots (with multiplicity) for ax^2+bx+c=0.
       Handles degenerate a==0 caller should have reduced degree.
    """
    if abs(a) == 0:
        if abs(b) == 0:
            return []
        return [ -c / b ]
    # stable quadratic formula
    disc = b*b - 4*a*c
    sqrt_disc = cmath.sqrt(disc)
    # use stable pair
    q = -0.5 * (b + cmath.copysign(sqrt_disc.real, b) + 1j*sqrt_disc.imag if (sqrt_disc.real != 0 or sqrt_disc.imag != 0) else b)
    # fallback simpler expression:
    # Use straightforward formula (acceptable here)
    r1 = (-b + sqrt_disc) / (2*a)
    r2 = (-b - sqrt_disc) / (2*a)
    return [r1, r2]

def solve_cubic(coeffs: List[complex]) -> List[complex]:
    """
    Solve cubic polynomial with coefficients given as list [a3, a2, a1, a0]
    Returns list of 3 roots (with multiplicity). Works for lower-degree via callers.
    """
    # normalize to monic
    if len(coeffs) != 4:
        raise ValueError("Expected 4 coefficients for cubic")
    a3, a2, a1, a0 = coeffs
    if abs(a3) == 0:
        # reduce to quadratic
        return solve_quadratic(a2, a1, a0)
    a = a2 / a3
    b = a1 / a3
    c = a0 / a3
    # depressed cubic x = y - a/3 -> y^3 + py + q = 0
    a_over_3 = a / 3.0
    p = b - a*a_over_3
    q = (2*a*a*a)/27.0 - (a*b)/3.0 + c

    # discriminant
    discriminant = (q/2.0)**2 + (p/3.0)**3

    roots = []
    if abs(discriminant) < 1e-14:
        discriminant = 0.0

    if discriminant > 0:
        # one real root and two complex
        sqrt_disc = math.sqrt(discriminant)
        A = _cbrt(-q/2.0 + sqrt_disc)
        B = _cbrt(-q/2.0 - sqrt_disc)
        y1 = A + B
        x1 = y1 - a_over_3
        # other two (complex)
        y2 = -(A+B)/2.0
        imag = (A - B) * math.sqrt(3)/2.0 if not isinstance(A, complex) and not isinstance(B, complex) else (A - B) * (cmath.sqrt(3)/2.0)
        x2 = y2 - a_over_3 + 1j*0  # placeholder
        # compute full roots directly via complex arithmetic
        # Better: compute cubic roots via complex arithmetic using cmath
        A_c = complex(-q/2.0 + discriminant**0.5) ** (1/3) if discriminant>=0 else complex(-q/2.0 + cmath.sqrt(discriminant)) ** (1/3)
        # Simpler stable approach: use cmath to get one real and two complex roots
        C1 = cmath.sqrt(discriminant)
        A = (-q/2.0 + C1) ** (1/3)
        B = (-q/2.0 - C1) ** (1/3)
        y1 = A + B
        x1 = y1 - a_over_3
        # cube roots of unity
        omega = complex(-0.5, math.sqrt(3)/2.0)
        x2 = omega*A + omega.conjugate()*B - a_over_3
        x3 = omega.conjugate()*A + omega*B - a_over_3
        roots = [x1, x2, x3]
    elif discriminant == 0:
        # multiple roots
        # then A = B = cube root(-q/2)
        A = _cbrt(-q/2.0)
        y1 = 2*A
        y2 = -A
        roots = [y1 - a_over_3, y2 - a_over_3, y2 - a_over_3]
    else:
        # three real roots: use trig solution
        rho = math.sqrt(-(p**3)/27.0)
        # alternative stable trig method
        phi = math.acos( max(-1.0, min(1.0, (-q/2.0) / math.sqrt(-(p**3)/27.0) )) )
        m = 2 * math.sqrt(-p/3.0)
        y1 = m * math.cos(phi/3.0)
        y2 = m * math.cos((phi + 2*math.pi)/3.0)
        y3 = m * math.cos((phi + 4*math.pi)/3.0)
        roots = [y1 - a_over_3, y2 - a_over_3, y3 - a_over_3]

    # ensure we return 3 roots
    if len(roots) == 2:
        roots.append(roots[1])
    if len(roots) == 1:
        roots = [roots[0], roots[0], roots[0]]
    return roots
