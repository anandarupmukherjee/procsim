"""Small numerical helpers (numpy only, so the core also runs under Pyodide)."""
import math


class ConvergenceError(RuntimeError):
    pass


def brent(f, a, b, tol=1e-9, maxiter=200):
    """Brent's method for a root of f in [a, b]; f(a) and f(b) must differ in sign."""
    fa, fb = f(a), f(b)
    if fa == 0:
        return a
    if fb == 0:
        return b
    if fa * fb > 0:
        raise ConvergenceError(f"Root not bracketed: f({a})={fa:.4g}, f({b})={fb:.4g}")
    c, fc = a, fa
    d = e = b - a
    for _ in range(maxiter):
        if fb * fc > 0:
            c, fc = a, fa
            d = e = b - a
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        tol1 = 2 * 1e-15 * abs(b) + 0.5 * tol
        m = 0.5 * (c - b)
        if abs(m) <= tol1 or fb == 0:
            return b
        if abs(e) >= tol1 and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                p, q = 2 * m * s, 1 - s
            else:
                q, r = fa / fc, fb / fc
                p = s * (2 * m * q * (q - r) - (b - a) * (r - 1))
                q = (q - 1) * (r - 1) * (s - 1)
            if p > 0:
                q = -q
            p = abs(p)
            if 2 * p < min(3 * m * q - abs(tol1 * q), abs(e * q)):
                e, d = d, p / q
            else:
                d = e = m
        else:
            d = e = m
        a, fa = b, fb
        b = b + d if abs(d) > tol1 else b + math.copysign(tol1, m)
        fb = f(b)
    raise ConvergenceError("Brent did not converge")


def bracket_root(f, x0, step, lo, hi, maxiter=60):
    """Expand outward from x0 until f changes sign; return (a, b)."""
    f0 = f(x0)
    a, b = x0, x0
    fa = fb = f0
    for _ in range(maxiter):
        a2, b2 = max(lo, a - step), min(hi, b + step)
        fa2, fb2 = f(a2), f(b2)
        if f0 * fa2 <= 0:
            return a2, a
        if f0 * fb2 <= 0:
            return b, b2
        a, b, fa, fb = a2, b2, fa2, fb2
        step *= 1.6
        if a <= lo and b >= hi:
            break
    raise ConvergenceError(f"Could not bracket root in [{lo}, {hi}]")
