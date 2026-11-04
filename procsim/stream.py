"""Material streams: fluid components (volatile + dissolved solutes) and particulate minerals."""
import math
import numpy as np

from .components import component, mineral, T_REF
from .thermo import make_package
from .numerics import brent, bracket_root, ConvergenceError

# Sieve series: 50 mm downwards in sqrt(2) steps, 24 sieves + pan.
# Class i (0..23) is material retained on sieve i: size in (sieve_i, sieve_{i-1}];
# class 0 is oversize (> 50 mm); class 24 is the pan (< smallest sieve).
SIEVES_UM = np.array([50000.0 * 2 ** (-i / 2) for i in range(24)])
NBINS = len(SIEVES_UM) + 1
# Representative (geometric-mean) size of each class, um
MID_UM = np.append(SIEVES_UM * 2 ** 0.25, SIEVES_UM[-1] / 2)


def passing_at_sieves(mass_by_bin):
    """Cumulative fraction passing each sieve."""
    m = np.asarray(mass_by_bin, float)
    tot = m.sum()
    if tot <= 0:
        return np.zeros(len(SIEVES_UM))
    cum = np.cumsum(m[::-1])[::-1]  # cum[i] = sum of classes i..end
    return cum[1:] / tot            # passing sieve i = classes i+1..end


def p_size(mass_by_bin, frac=0.8):
    """Size [um] at which `frac` of the mass passes (log-linear interpolation)."""
    pas = passing_at_sieves(mass_by_bin)
    if pas.size == 0 or np.sum(mass_by_bin) <= 0:
        return float("nan")
    if pas[0] < frac:
        return float(SIEVES_UM[0])
    for i in range(len(pas) - 1):
        hi, lo = pas[i], pas[i + 1]
        if hi >= frac >= lo:
            if hi == lo:
                return float(SIEVES_UM[i + 1])
            t = (frac - lo) / (hi - lo)
            return float(math.exp(math.log(SIEVES_UM[i + 1]) + t * (math.log(SIEVES_UM[i]) - math.log(SIEVES_UM[i + 1]))))
    return float(SIEVES_UM[-1])


def rrsb_fractions(p80_um, spread=1.0):
    """Mass fractions per class for a Rosin-Rammler distribution with given P80."""
    d63 = p80_um / (math.log(5.0)) ** (1.0 / spread)
    passing = 1.0 - np.exp(-(SIEVES_UM / d63) ** spread)
    frac = np.empty(NBINS)
    frac[0] = 1.0 - passing[0]                 # retained on the top sieve
    frac[1:-1] = passing[:-1] - passing[1:]    # retained on sieve i
    frac[-1] = passing[-1]                     # pan
    frac = np.clip(frac, 0, None)
    return frac / frac.sum()


class Context:
    """Shared definition of what a flowsheet's streams carry."""

    def __init__(self, components, minerals=(), package="Peng-Robinson"):
        self.names = list(components)
        self.comps = [component(n) for n in self.names]
        self.vol = np.array([i for i, c in enumerate(self.comps) if not c.solute], int)
        self.sol = np.array([i for i, c in enumerate(self.comps) if c.solute], int)
        self.mw = np.array([c.mw for c in self.comps]) / 1000.0  # kg/mol
        self.minerals = list(minerals)
        self.mins = [mineral(m) for m in self.minerals]
        self.min_cp = np.array([m.cp for m in self.mins]) if self.mins else np.zeros(0)
        self.package_name = package
        self.pkg = make_package(package, [self.names[i] for i in self.vol]) if len(self.vol) else None

    def index(self, name):
        return self.names.index(name)

    def mindex(self, name):
        return self.minerals.index(name)


class Stream:
    def __init__(self, ctx, T=298.15, P=101325.0, F=None, S=None, name=""):
        self.ctx = ctx
        self.name = name
        self.T, self.P = float(T), float(P)
        self.F = np.zeros(len(ctx.names)) if F is None else np.asarray(F, float).copy()
        self.S = np.zeros((len(ctx.minerals), NBINS)) if S is None else np.asarray(S, float).copy()
        self.beta = 0.0
        self.x = self.y = None
        self.phase = "L"
        self.H = 0.0
        self._solved = False

    # ---- construction --------------------------------------------------
    @classmethod
    def from_spec(cls, ctx, spec, name=""):
        """spec: {T [K], P [Pa], flows {comp: mol/s}, solids {mineral: {mass: kg/s, p80: um, spread}}}"""
        s = cls(ctx, spec.get("T", 298.15), spec.get("P", 101325.0), name=name)
        for k, v in (spec.get("flows") or {}).items():
            s.F[ctx.index(k)] = float(v)
        for k, v in (spec.get("solids") or {}).items():
            if isinstance(v, (int, float)):
                v = {"mass": v}
            frac = rrsb_fractions(float(v.get("p80", 10000.0)), float(v.get("spread", 1.0)))
            s.S[ctx.mindex(k)] = float(v.get("mass", 0.0)) * frac
        s.flash_TP()
        return s

    def copy(self, name=None):
        s = Stream(self.ctx, self.T, self.P, self.F, self.S, self.name if name is None else name)
        s.beta, s.x, s.y, s.phase, s.H, s._solved = self.beta, self.x, self.y, self.phase, self.H, self._solved
        return s

    # ---- totals --------------------------------------------------------
    @property
    def n_vol(self):
        return float(self.F[self.ctx.vol].sum()) if len(self.ctx.vol) else 0.0

    @property
    def molar_flow(self):
        return float(self.F.sum())

    @property
    def fluid_mass(self):
        return float(self.F @ self.ctx.mw)

    @property
    def solids_mass(self):
        return float(self.S.sum())

    @property
    def mass(self):
        return self.fluid_mass + self.solids_mass

    # ---- thermodynamics ------------------------------------------------
    def _h_nonvol(self, T):
        h = 0.0
        ctx = self.ctx
        for i in ctx.sol:
            c = ctx.comps[i]
            h += self.F[i] * (c.hf + c.cp_liq * (T - T_REF))
        if self.S.size:
            h += float(self.S.sum(axis=1) @ ctx.min_cp) * (T - T_REF)
        return h

    def _flash_vol(self, T):
        ctx = self.ctx
        n = self.n_vol
        if n <= 1e-15 or ctx.pkg is None:
            return None, 0.0
        fr = ctx.pkg.flash_TP(T, self.P, self.F[ctx.vol] / n)
        return fr, n * ctx.pkg.enthalpy(fr)

    def flash_TP(self, T=None, P=None):
        if T is not None:
            self.T = float(T)
        if P is not None:
            self.P = float(P)
        fr, hv = self._flash_vol(self.T)
        if fr is None:
            self.beta, self.x, self.y, self.phase = 0.0, None, None, "L"
        else:
            self.beta, self.x, self.y, self.phase = fr.beta, fr.x, fr.y, fr.phase
        self.H = hv + self._h_nonvol(self.T)
        self._solved = True
        return self

    def enthalpy_at(self, T):
        fr, hv = self._flash_vol(T)
        return hv + self._h_nonvol(T)

    def flash_PH(self, H, P=None, T_guess=None):
        if P is not None:
            self.P = float(P)
        T0 = T_guess or self.T or 300.0
        scale = max(abs(H), 1.0)

        def f(T):
            return (self.enthalpy_at(T) - H) / scale
        try:
            a, b = bracket_root(f, T0, 5.0, 30.0, 3000.0)
            T = brent(f, a, b, tol=1e-6)
        except ConvergenceError as e:
            raise ConvergenceError(f"PH flash failed for stream {self.name}: {e}")
        return self.flash_TP(T)

    def vapour_fraction(self):
        return self.beta

    # ---- phase split ---------------------------------------------------
    def split_phases(self):
        """Return (vapour, liquid) streams at the stream's T, P. Solids and solutes go to the liquid."""
        ctx = self.ctx
        V = Stream(ctx, self.T, self.P)
        L = Stream(ctx, self.T, self.P, self.F, self.S)
        if self.beta > 0 and self.x is not None:
            n = self.n_vol
            if self.beta >= 1.0:
                V.F[ctx.vol] = self.F[ctx.vol]
                L.F[ctx.vol] = 0.0
            else:
                V.F[ctx.vol] = self.beta * n * self.y
                L.F[ctx.vol] = (1 - self.beta) * n * self.x
        V.flash_TP()
        L.flash_TP()
        return V, L

    # ---- reporting -----------------------------------------------------
    def psd_p80(self):
        tot = self.S.sum(axis=0)
        return p_size(tot, 0.8) if tot.sum() > 0 else None

    def summary(self):
        ctx = self.ctx
        d = {
            "name": self.name,
            "T_K": self.T, "T_C": self.T - 273.15, "P_bar": self.P / 1e5,
            "vapour_fraction": self.beta, "phase": self.phase,
            "molar_flow_mol_s": self.molar_flow,
            "mass_flow_kg_s": self.mass,
            "enthalpy_kW": self.H / 1e3,
            "flows_mol_s": {n: float(v) for n, v in zip(ctx.names, self.F)},
            "mole_fractions": {n: float(v / self.molar_flow) if self.molar_flow > 0 else 0.0
                               for n, v in zip(ctx.names, self.F)},
        }
        if ctx.minerals:
            d["solids_kg_s"] = {m: float(self.S[i].sum()) for i, m in enumerate(ctx.minerals)}
            d["solids_total_kg_s"] = self.solids_mass
            d["solids_pct_w"] = 100 * self.solids_mass / self.mass if self.mass > 0 else 0.0
            d["P80_um"] = self.psd_p80()
            d["psd_passing"] = passing_at_sieves(self.S.sum(axis=0)).tolist() if self.solids_mass > 0 else None
        return d
