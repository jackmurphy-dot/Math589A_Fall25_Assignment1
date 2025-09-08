# quartic_solver.py
import math
import cmath
from typing import List
from cubic_solver import solve_cubic, solve_quadratic

def _is_zero(x, tol=1e-12):
    return abs(x) < tol

def solve_polynomial(coeffs: List[complex]) -> List[complex]:
    """
    Solve polynomial given by coeffs list [a_n, a_{n-1}, ..., a_0].
    Supports degrees 0..4. Returns list of roots (complex) with multiplicity.
    """
    # strip leading zeros
    i = 0
    while i < len(coeffs) and _is_zero(coeffs[i]):
        i += 1
    if i == len(coeffs):
        return []  # zero polynomial: infinite roots, return empty
    coeffs = coeffs[i:]
    deg = len(coeffs) - 1
    if deg == 0:
        return []
    if deg == 1:
        a, b = coeffs
        return [ -b / a ]
    if deg == 2:
        a,b,c = coeffs
        return solve_quadratic(a,b,c)
    if deg == 3:
        return solve_cubic(coeffs)
    if deg == 4:
        a4,a3,a2,a1,a0 = coeffs
        # normalize to monic: x^4 + px^3 + qx^2 + rx + s = 0
        p = a3 / a4
        q = a2 / a4
        r = a1 / a4
        s = a0 / a4

        # depressed quartic via x = y - p/4
        alpha = -3*(p**2)/8 + q
        beta  = p**3/8 - p*q/2 + r
        gamma = -3*(p**4)/256 + p**2*q/16 - p*r/4 + s

        roots = []
        if _is_zero(beta):
            # biquadratic: y^4 + alpha*y^2 + gamma = 0 => set z=y^2 solve z^2+alpha*z+gamma=0
            zs = solve_quadratic(1, alpha, gamma)
            for z in zs:
                # possibly two y for each z
                sqrt_z = cmath.sqrt(z)
                roots.append(sqrt_z - p/4.0)
                roots.append(-sqrt_z - p/4.0)
            return roots

        # general case: use resolvent cubic (Ferrari)
        # construct resolvent cubic for u: u^3 + (alpha/2) u^2 + ((alpha^2 - 4*gamma)/16) u - (beta^2/64) = 0
        # Another common construction: let resolvent cubic be for y:
        # y^3 - (alpha/2) y^2 - gamma y/2 + (gamma*alpha/2 - beta^2/8) = 0
        # We'll use standard approach: build cubic for m (see many references)
        # Using formulas: let resolvent cubic coefficients:
        # Following one of standard forms:
        A = 1
        B = -alpha/2.0
        C = -gamma
        D = (alpha*gamma - (beta**2)/4.0)/2.0

        cubic_coeffs = [A, B, C, D]
        # solve cubic
        u_roots = solve_cubic(cubic_coeffs)
        # pick one root u that gives positive value for something or is convenient
        # choose u with largest real part (heuristic)
        u = max(u_roots, key=lambda z: (z.real, z.imag))

        # compute some helpers
        # compute R = sqrt(0.25*p^2 - q + u) ... but in depressed form p here is alpha etc.
        # Use classical Ferrari back-substitution
        # Let:
        # R = sqrt( (p^2)/4 - q + u ) but careful: using depressed variables:
        # Many algebraic forms; below is direct standard formula using alpha,beta,gamma and u.
        try:
            R = cmath.sqrt(0.25*p*p - q + u)  # p here original p; this formula may work but ensure consistent
        except Exception:
            R = cmath.sqrt(u - alpha/4.0)  # fallback

        if _is_zero(R):
            D1 = cmath.sqrt(3.0/4.0*p*p - R*R - 2*q + (4.0*p*q - 8.0*r - p**3)/(4.0*R))
            E1 = cmath.sqrt(3.0/4.0*p*p - R*R - 2*q - (4.0*p*q - 8.0*r - p**3)/(4.0*R))
        # Alternative classical recipe (safer):
        # Work with depressed quartic y^4 + alpha y^2 + beta y + gamma = 0, and u is root of resolvent cubic:
        # Then compute:
        U = u
        # compute two quadratic factors: y^2 +/- sqrt(2U - alpha) y + (U^2 - gamma +/- beta / sqrt(2U - alpha))/2 = 0
        sqrt_term = cmath.sqrt(2*U - alpha)
        # compute the two quadratics' coefficients
        # first quadratic: y^2 + sqrt_term*y + (U + ( -beta/(2*sqrt_term) ) - alpha/2)
        # derived formula below:
        # But to avoid division by zero, handle case sqrt_term==0
        if abs(sqrt_term) == 0:
            # fallback to direct splitting using numeric quadratics via Ferrari alternative
            # compute coefficients via Ferrari by constructing quadratic resolvents
            # We'll form the two quadratics by solving for their coefficients numerically via solving linear system.
            # Simpler fallback: compute all four roots by using numpy-esque approach via companion matrix? Not allowed.
            # Instead use direct Ferrari numeric formulas (less pretty).
            # Use slightly different resolvent: compute as done in many references:
            p1 = math.sqrt((alpha*alpha)/4 - gamma)
            # compute quadratics coefficients (approx)
            q1 = 0.5*(alpha/2 + p1)
            q2 = 0.5*(alpha/2 - p1)
            quad1 = solve_quadratic(1, 0, q1)
            quad2 = solve_quadratic(1, 0, q2)
            roots = [z - p/4.0 for z in quad1 + quad2]
            return roots
        else:
            # safe path
            y_coeff = sqrt_term
            # compute constant terms for quadratics
            c1 = (U - gamma + (beta / (2.0*sqrt_term))) / 2.0
            c2 = (U - gamma - (beta / (2.0*sqrt_term))) / 2.0
            # quadratics: y^2 + y_coeff*y + c1 = 0 and y^2 - y_coeff*y + c2 = 0
            qroots1 = solve_quadratic(1.0, y_coeff, c1)
            qroots2 = solve_quadratic(1.0, -y_coeff, c2)
            roots = [z - p/4.0 for z in (qroots1 + qroots2)]
            return roots

    raise NotImplementedError("Unhandled polynomial degree")

if __name__ == "__main__":
    # simple manual test
    # e.g., (x-1)^4 -> roots at 1 multiplicity 4
    roots = solve_polynomial([1, -4, 6, -4, 1])
    print("roots:", roots)
