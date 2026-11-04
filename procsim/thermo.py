"""Thermodynamic property packages for the volatile (vapour-liquid) part of a stream.

Two packages share one interface:
  * PengRobinson: PR cubic EOS, van der Waals mixing rules with binary kij,
    fugacity-based K-values and enthalpy departure functions.
  * Ideal: Raoult's law with Lee-Kesler vapour pressure; liquid enthalpy from
    the ideal-gas enthalpy minus a Pitzer-type latent heat.

Enthalpies are on an ideal-gas formation basis (elements at 298.15 K), so heats of
reaction appear automatically in reactor energy balances.
"""
import math
import numpy as np

from .components import R, component, kij
from .numerics import brent, bracket_root, ConvergenceError

SQ2 = math.sqrt(2.0)


def rachford_rice(z, K):
    """Solve sum z(K-1)/(1+b(K-1)) = 0 for b, allowing negative-flash values.

    Returns None when no root exists (all K>1 or all K<1)."""
    m = z > 1e-300
    z, K = z[m], K[m]
    Km, Kx = K.min(), K.max()
    if Kx <= 1.0 + 1e-12 or Km >= 1.0 - 1e-12:
        return None
    lo = 1.0 / (1.0 - Kx)
    hi = 1.0 / (1.0 - Km)

    def g(b):
        return float(np.sum(z * (K - 1.0) / (1.0 + b * (K - 1.0))))
    eps = 1e-10 * (hi - lo)
    a, c = lo + eps, hi - eps
    while g(a) < 0:
        eps *= 1e-3
        a2 = lo + eps
        if a2 == a or a2 <= lo:
            return a  # root sits on the lower asymptote (trace component)
        a = a2
    eps = 1e-10 * (hi - lo)
    while g(c) > 0:
        eps *= 1e-3
        c2 = hi - eps
        if c2 == c or c2 >= hi:
            return c
        c = c2
    return brent(g, a, c, tol=1e-14)


class FlashResult:
    __slots__ = ("T", "P", "beta", "x", "y", "phase")

    def __init__(self, T, P, beta, x, y, phase):
        self.T, self.P, self.beta, self.x, self.y, self.phase = T, P, beta, x, y, phase


class ThermoPackage:
    name = "base"

    def __init__(self, names):
        self.names = list(names)
        self.comps = [component(n) for n in self.names]
        self.n = len(self.comps)
        self.tc = np.array([c.tc for c in self.comps])
        self.pc = np.array([c.pc for c in self.comps])
        self.w = np.array([c.omega for c in self.comps])

    # ---- to be provided by subclasses -------------------------------------
    def K(self, T, P, x, y):
        raise NotImplementedError

    def h(self, T, P, x, phase):
        """Molar enthalpy [J/mol] of a phase of composition x ('L', 'V' or 'auto')."""
        raise NotImplementedError

    # ---- shared -----------------------------------------------------------
    def h_ig(self, T, x):
        return float(sum(xi * c.h_ig(T) for xi, c in zip(x, self.comps)))

    def wilson_K(self, T, P):
        return self.pc / P * np.exp(5.373 * (1.0 + self.w) * (1.0 - self.tc / T))

    def flash_TP(self, T, P, z):
        z = np.asarray(z, float)
        z = z / z.sum()
        if self.n == 1:
            return self._flash_pure(T, P, z)
        K = self.wilson_K(T, P)
        beta = 0.5
        x, y = z.copy(), z.copy()
        for it in range(300):
            b = rachford_rice(z, K)
            if b is None:
                # all K on one side: single phase (pick by sign)
                beta = 1.0 if K.min() >= 1.0 else 0.0
                x = y = z.copy()
                # refine K once more against the single phase
                Knew = self.K(T, P, z, z)
                if np.allclose(Knew, K, rtol=1e-8):
                    break
                K = Knew
                if rachford_rice(z, K) is None:
                    break
                continue
            beta = b
            d = 1.0 + beta * (K - 1.0)
            x = z / d
            y = K * x
            x, y = x / x.sum(), y / y.sum()
            Knew = self.K(T, P, x, y)
            err = np.max(np.abs(np.log(Knew / K)))
            K = Knew
            if err < 1e-10:
                break
            if np.max(np.abs(np.log(K))) < 1e-4:  # trivial solution
                beta = 1.0 if self._prefers_vapour(T, P, z) else 0.0
                x = y = z.copy()
                break
        if beta <= 0.0:
            return FlashResult(T, P, 0.0, z.copy(), z.copy(), "L")
        if beta >= 1.0:
            return FlashResult(T, P, 1.0, z.copy(), z.copy(), "V")
        return FlashResult(T, P, float(beta), x, y, "VL")

    def _prefers_vapour(self, T, P, z):
        return float(np.sum(z * self.tc)) < T

    def _flash_pure(self, T, P, z):
        try:
            Tb = self.bubble_T(z, P)
        except ConvergenceError:
            Tb = None
        if Tb is None or T > Tb:
            return FlashResult(T, P, 1.0, z, z, "V")
        return FlashResult(T, P, 0.0, z, z, "L")

    def enthalpy(self, fr):
        """Enthalpy [J/mol of volatile] for a FlashResult."""
        if fr.phase == "V":
            return self.h(fr.T, fr.P, fr.y, "V")
        if fr.phase == "L":
            return self.h(fr.T, fr.P, fr.x, "L")
        return fr.beta * self.h(fr.T, fr.P, fr.y, "V") + (1 - fr.beta) * self.h(fr.T, fr.P, fr.x, "L")

    def _boundary_T(self, z, P, target, T0):
        """Temperature where the flash vapour fraction crosses `target` (robust for mixtures)."""
        def f(T):
            return self.flash_TP(T, P, z).beta - target
        a, b = bracket_root(f, T0, max(2.0, 0.02 * T0), 20.0, 3000.0)
        return brent(f, a, b, tol=1e-5)

    def _wilson_T(self, z, P, bubble=True):
        def gw(T):
            K = self.wilson_K(T, P)
            return math.log(np.sum(z * K)) if bubble else -math.log(np.sum(z / K))
        try:
            a, b = bracket_root(gw, 300.0, 20.0, 20.0, 2500.0)
            return brent(gw, a, b, tol=1e-6)
        except ConvergenceError:
            return 300.0

    def bubble_T(self, x, P, T0=None):
        x = np.asarray(x, float) / np.sum(x)
        if np.count_nonzero(x > 1e-12) > 1:
            return self._boundary_T(x, P, 1e-7, T0 or self._wilson_T(x, P, True))

        def g(T):
            y = x.copy()
            K = self.wilson_K(T, P)
            for _ in range(50):
                y = K * x / np.sum(K * x)
                Kn = self.K(T, P, x, y)
                if np.max(np.abs(Kn - K)) < 1e-10:
                    K = Kn
                    break
                K = Kn
            return math.log(np.sum(K * x))
        return self._solve_T(g, x, P, T0)

    def dew_T(self, y, P, T0=None):
        y = np.asarray(y, float) / np.sum(y)
        if np.count_nonzero(y > 1e-12) > 1:
            return self._boundary_T(y, P, 1 - 1e-7, T0 or self._wilson_T(y, P, False))

        def g(T):
            x = y.copy()
            K = self.wilson_K(T, P)
            for _ in range(50):
                x = y / K / np.sum(y / K)
                Kn = self.K(T, P, x, y)
                if np.max(np.abs(Kn - K)) < 1e-10:
                    K = Kn
                    break
                K = Kn
            return -math.log(np.sum(y / K))
        return self._solve_T(g, y, P, T0)

    def _solve_T(self, g, z, P, T0):
        if T0 is None:
            # Wilson estimate
            def gw(T):
                return math.log(np.sum(z * self.wilson_K(T, P)))
            try:
                a, b = bracket_root(gw, 300.0, 20.0, 20.0, 2500.0)
                T0 = brent(gw, a, b, tol=1e-6)
            except ConvergenceError:
                T0 = 300.0
        a, b = bracket_root(g, T0, max(2.0, 0.02 * T0), 20.0, 3000.0)
        return brent(g, a, b, tol=1e-7)


class Ideal(ThermoPackage):
    name = "Ideal"

    def psat(self, T):
        Tr = T / self.tc
        f0 = 5.92714 - 6.09648 / Tr - 1.28862 * np.log(Tr) + 0.169347 * Tr ** 6
        f1 = 15.2518 - 15.6875 / Tr - 13.4721 * np.log(Tr) + 0.43577 * Tr ** 6
        return self.pc * np.exp(np.minimum(f0 + self.w * f1, 50.0))

    def K(self, T, P, x, y):
        return self.psat(T) / P

    def wilson_K(self, T, P):
        return self.psat(T) / P

    def dhvap(self, T):
        Tr = np.minimum(T / self.tc, 1.0)
        return R * self.tc * (7.08 * (1 - Tr) ** 0.354 + 10.95 * self.w * (1 - Tr) ** 0.456)

    def h(self, T, P, x, phase):
        x = np.asarray(x, float)
        hig = self.h_ig(T, x)
        if phase == "V":
            return hig
        return hig - float(np.sum(x * self.dhvap(T)))


class PengRobinson(ThermoPackage):
    name = "Peng-Robinson"

    def __init__(self, names):
        super().__init__(names)
        self.kappa = 0.37464 + 1.54226 * self.w - 0.26992 * self.w ** 2
        self.ac = 0.45724 * R ** 2 * self.tc ** 2 / self.pc
        self.b = 0.07780 * R * self.tc / self.pc
        self.k = np.array([[kij(a, b) for b in self.names] for a in self.names])

    def _a(self, T):
        sTr = np.sqrt(T / self.tc)
        alpha = (1 + self.kappa * (1 - sTr)) ** 2
        ai = self.ac * alpha
        dai = -self.ac * self.kappa * (1 + self.kappa * (1 - sTr)) / np.sqrt(T * self.tc)
        sa = np.sqrt(ai)
        aij = (1 - self.k) * np.outer(sa, sa)
        daij = (1 - self.k) * 0.5 * (np.outer(dai / sa, sa) + np.outer(sa, dai / sa))
        return aij, daij

    @staticmethod
    def _roots(A, B):
        coeffs = [1.0, -(1 - B), A - 3 * B ** 2 - 2 * B, -(A * B - B ** 2 - B ** 3)]
        r = np.roots(coeffs)
        r = np.real(r[np.abs(np.imag(r)) < 1e-9])
        r = r[r > B]
        if r.size == 0:
            raise ConvergenceError("No physical PR root")
        return r

    def _state(self, T, P, x, phase):
        x = np.asarray(x, float)
        aij, daij = self._a(T)
        am = float(x @ aij @ x)
        bm = float(x @ self.b)
        A = am * P / (R * T) ** 2
        B = bm * P / (R * T)
        roots = self._roots(A, B)
        if phase == "L":
            Z = roots.min()
        elif phase == "V":
            Z = roots.max()
        else:
            Z = self._min_g_root(roots, A, B)
        return x, aij, daij, am, bm, A, B, Z

    @staticmethod
    def _min_g_root(roots, A, B):
        if roots.size == 1:
            return roots[0]
        def g(Z):  # residual Gibbs (dimensionless, up to constants)
            return Z - 1 - math.log(Z - B) - A / (2 * SQ2 * B) * math.log((Z + (1 + SQ2) * B) / (Z + (1 - SQ2) * B))
        return min(roots, key=g)

    def lnphi(self, T, P, x, phase):
        x, aij, _, am, bm, A, B, Z = self._state(T, P, x, phase)
        bi = self.b / bm
        sa = (aij @ x) / am
        L = math.log((Z + (1 + SQ2) * B) / (Z + (1 - SQ2) * B))
        return bi * (Z - 1) - math.log(Z - B) - A / (2 * SQ2 * B) * (2 * sa - bi) * L

    def K(self, T, P, x, y):
        return np.exp(self.lnphi(T, P, x, "L") - self.lnphi(T, P, y, "V"))

    def h(self, T, P, x, phase):
        x, aij, daij, am, bm, A, B, Z = self._state(T, P, x, phase)
        dam = float(x @ daij @ x)
        L = math.log((Z + (1 + SQ2) * B) / (Z + (1 - SQ2) * B))
        hdep = R * T * (Z - 1) + (T * dam - am) / (2 * SQ2 * bm) * L
        return self.h_ig(T, x) + hdep

    def flash_TP(self, T, P, z):
        fr = super().flash_TP(T, P, z)
        if fr.phase in ("L", "V"):
            fr.phase = fr.phase  # single phase: use min-Gibbs root for properties
        return fr

    def enthalpy(self, fr):
        if fr.phase in ("L", "V"):
            return self.h(fr.T, fr.P, fr.x, "auto")
        return super().enthalpy(fr)

    def _prefers_vapour(self, T, P, z):
        x, aij, _, am, bm, A, B, Z = self._state(T, P, z, "auto")
        return Z > 0.3


PACKAGES = {"Peng-Robinson": PengRobinson, "PR": PengRobinson, "Ideal": Ideal}


def make_package(name, components):
    try:
        return PACKAGES[name](components)
    except KeyError:
        raise KeyError(f"Unknown property package '{name}'. Use one of {sorted(set(PACKAGES))}")
