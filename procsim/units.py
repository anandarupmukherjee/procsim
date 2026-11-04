"""Unit operation models.

Every unit declares its ports and parameter schema, so the flowsheet engine and the
web UI are driven from the same definitions. Several streams may be connected to
one inlet port; they are mixed adiabatically before the unit is calculated.

User-facing units: temperature degC, pressure bar, duty/power kW, size um, time min.
"""
import math
import re
import numpy as np

from .stream import Stream, MID_UM, NBINS, SIEVES_UM, passing_at_sieves, rrsb_fractions, p_size
from .numerics import brent, ConvergenceError

K0 = 273.15
UNIT_TYPES = {}


def register(cls):
    UNIT_TYPES[cls.TYPE] = cls
    return cls


def P_(name, kind, default, units="", help="", options=None, when=None):
    d = {"name": name, "kind": kind, "default": default, "units": units, "help": help}
    if options is not None:
        d["options"] = options
    if when is not None:
        d["when"] = when
    return d


def mix(streams, ctx, P=None, name="mix"):
    streams = [s for s in streams if s is not None]
    if not streams:
        raise ValueError("No inlet streams")
    if len(streams) == 1 and P is None:
        return streams[0].copy(name)
    out = Stream(ctx, name=name)
    out.F = sum(s.F for s in streams)
    out.S = sum(s.S for s in streams)
    H = sum(s.H for s in streams)
    Pm = min(s.P for s in streams) if P is None else P
    Tg = sum(s.T * s.mass for s in streams) / max(sum(s.mass for s in streams), 1e-30)
    out.P = Pm
    out.T = Tg
    if out.molar_flow <= 0 and out.solids_mass <= 0:
        return out.flash_TP()
    return out.flash_PH(H, Pm, T_guess=Tg)


def parse_table(text, value_keys):
    """Parse lines like 'chalcopyrite: k=1.2, Rmax=0.95' or 'chalcopyrite 1.2 0.95'."""
    out = {}
    for line in (text or "").replace(";", "\n").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if ":" in line:
            key, rest = line.split(":", 1)
        else:
            parts = line.split()
            key, rest = parts[0], " ".join(parts[1:])
        key = key.strip()
        vals = {}
        named = re.findall(r"([A-Za-z_]+)\s*=\s*([-+0-9.eE]+)", rest)
        if named:
            for k, v in named:
                vals[k] = float(v)
        else:
            nums = [float(v) for v in re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?", rest)]
            for k, v in zip(value_keys, nums):
                vals[k] = v
        out[key] = vals
    return out


class Unit:
    TYPE = "Unit"
    CATEGORY = "General"
    LABEL = ""
    DESCRIPTION = ""
    INPORTS = ["in"]
    OUTPORTS = ["out"]
    PARAMS = []

    def __init__(self, uid, params=None):
        self.id = uid
        p = {d["name"]: d["default"] for d in self.PARAMS}
        p.update(params or {})
        self.p = p
        self.results = {}

    @classmethod
    def schema(cls):
        return {"type": cls.TYPE, "category": cls.CATEGORY, "label": cls.LABEL or cls.TYPE,
                "description": cls.DESCRIPTION, "inports": cls.INPORTS, "outports": cls.OUTPORTS,
                "params": cls.PARAMS}

    def num(self, k):
        v = self.p.get(k)
        return None if v in (None, "") else float(v)

    def feed(self, inlets, ctx, port="in"):
        return mix(inlets.get(port, []), ctx, name=f"{self.id}.{port}")

    def run(self, inlets, ctx):
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Fluid units
# ---------------------------------------------------------------------------
@register
class Mixer(Unit):
    TYPE, CATEGORY, LABEL = "Mixer", "Fluid", "Mixer"
    DESCRIPTION = "Adiabatic mixing of any number of streams."
    PARAMS = [P_("P_out", "number", "", "bar", "Outlet pressure (blank = lowest inlet)")]

    def run(self, inlets, ctx):
        out = mix(inlets.get("in", []), ctx, P=None if self.num("P_out") is None else self.num("P_out") * 1e5)
        self.results = {"Outlet T [C]": out.T - K0}
        return {"out": out}


@register
class Splitter(Unit):
    TYPE, CATEGORY, LABEL = "Splitter", "Fluid", "Splitter"
    DESCRIPTION = "Splits a stream into two of identical composition and state."
    OUTPORTS = ["out1", "out2"]
    PARAMS = [P_("frac1", "number", 0.5, "-", "Fraction of the inlet sent to out1")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        a = min(max(self.num("frac1"), 0.0), 1.0)
        o1, o2 = f.copy(), f.copy()
        o1.F, o1.S, o1.H = f.F * a, f.S * a, f.H * a
        o2.F, o2.S, o2.H = f.F * (1 - a), f.S * (1 - a), f.H * (1 - a)
        self.results = {"Split to out1": a}
        return {"out1": o1, "out2": o2}


@register
class Heater(Unit):
    TYPE, CATEGORY, LABEL = "Heater", "Fluid", "Heater / cooler"
    DESCRIPTION = "Sets outlet temperature (or applies a duty) with an optional pressure drop."
    PARAMS = [P_("mode", "select", "temperature", "", "Specification", ["temperature", "duty", "vapour fraction"]),
              P_("T_out", "number", 100.0, "C", "Outlet temperature", when={"mode": "temperature"}),
              P_("duty", "number", 0.0, "kW", "Heat added (negative = cooling)", when={"mode": "duty"}),
              P_("vf", "number", 1.0, "-", "Outlet vapour fraction (0 = bubble, 1 = dew)", when={"mode": "vapour fraction"}),
              P_("dP", "number", 0.0, "bar", "Pressure drop")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        out = f.copy(f"{self.id}.out")
        P = f.P - (self.num("dP") or 0.0) * 1e5
        mode = self.p["mode"]
        if mode == "temperature":
            out.flash_TP(self.num("T_out") + K0, P)
        elif mode == "duty":
            out.flash_PH(f.H + self.num("duty") * 1e3, P)
        else:
            vf = min(max(self.num("vf"), 0.0), 1.0)
            out.P = P
            n = out.n_vol
            if n <= 0:
                raise ValueError("Vapour fraction spec needs volatile components")
            z = out.F[ctx.vol] / n
            if vf <= 0.0:
                T = ctx.pkg.bubble_T(z, P)
            elif vf >= 1.0:
                T = ctx.pkg.dew_T(z, P)
            else:
                Tb, Td = ctx.pkg.bubble_T(z, P), ctx.pkg.dew_T(z, P)
                T = brent(lambda t: ctx.pkg.flash_TP(t, P, z).beta - vf, Tb + 1e-6, Td - 1e-6, tol=1e-6)
            out.flash_TP(T, P)
        Q = out.H - f.H
        self.results = {"Duty [kW]": Q / 1e3, "Outlet T [C]": out.T - K0, "Outlet vapour fraction": out.beta}
        self.duty = Q
        return {"out": out}


@register
class Valve(Unit):
    TYPE, CATEGORY, LABEL = "Valve", "Fluid", "Valve"
    DESCRIPTION = "Isenthalpic pressure let-down (Joule-Thomson)."
    PARAMS = [P_("P_out", "number", 1.01325, "bar", "Outlet pressure")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        out = f.copy(f"{self.id}.out").flash_PH(f.H, self.num("P_out") * 1e5)
        self.results = {"Outlet T [C]": out.T - K0, "Outlet vapour fraction": out.beta}
        return {"out": out}


@register
class Flash(Unit):
    TYPE, CATEGORY, LABEL = "Flash", "Fluid", "Flash drum"
    DESCRIPTION = "Two-phase vapour-liquid separator. Solids and dissolved species leave with the liquid."
    OUTPORTS = ["vapour", "liquid"]
    PARAMS = [P_("mode", "select", "T and P", "", "Specification", ["T and P", "adiabatic", "duty"]),
              P_("T", "number", 25.0, "C", "Drum temperature", when={"mode": "T and P"}),
              P_("P", "number", 1.01325, "bar", "Drum pressure"),
              P_("duty", "number", 0.0, "kW", "Heat added", when={"mode": "duty"})]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        s = f.copy(f"{self.id}.drum")
        P = self.num("P") * 1e5
        if self.p["mode"] == "T and P":
            s.flash_TP(self.num("T") + K0, P)
        elif self.p["mode"] == "adiabatic":
            s.flash_PH(f.H, P)
        else:
            s.flash_PH(f.H + self.num("duty") * 1e3, P)
        V, L = s.split_phases()
        Q = s.H - f.H
        self.results = {"Duty [kW]": Q / 1e3, "T [C]": s.T - K0, "Vapour fraction": s.beta}
        return {"vapour": V, "liquid": L}


@register
class ComponentSplitter(Unit):
    TYPE, CATEGORY, LABEL = "ComponentSplitter", "Fluid", "Component separator"
    DESCRIPTION = "Black-box separator: specify the fraction of each component (or mineral) sent to out1."
    OUTPORTS = ["out1", "out2"]
    PARAMS = [P_("splits", "text", "water: 0.99", "", "Lines 'name: fraction to out1'; unlisted species go to out2"),
              P_("default", "number", 0.0, "-", "Fraction to out1 for unlisted species")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        tab = {k: v.get("frac", v.get("f", next(iter(v.values()), 0.0)) if v else 0.0)
               for k, v in parse_table(self.p["splits"], ["frac"]).items()}
        d = self.num("default") or 0.0
        o1, o2 = Stream(ctx, f.T, f.P), Stream(ctx, f.T, f.P)
        for i, n in enumerate(ctx.names):
            a = tab.get(n, d)
            o1.F[i], o2.F[i] = f.F[i] * a, f.F[i] * (1 - a)
        for i, m in enumerate(ctx.minerals):
            a = tab.get(m, d)
            o1.S[i], o2.S[i] = f.S[i] * a, f.S[i] * (1 - a)
        o1.flash_TP()
        o2.flash_TP()
        Q = o1.H + o2.H - f.H
        self.results = {"Duty to keep T [kW]": Q / 1e3}
        return {"out1": o1, "out2": o2}


def _species_index(ctx, token):
    token = token.strip()
    for i, c in enumerate(ctx.comps):
        if token.lower() in (c.name.lower(), c.formula.lower()):
            return i
    raise KeyError(f"Species '{token}' is not in the component list")


def parse_reactions(text, ctx):
    """Lines: 'CO + 2 H2 -> CH3OH : CO : 0.3'  (stoichiometry : key component : conversion)."""
    rxns = []
    for line in (text or "").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.split(":")]
        if len(parts) != 3:
            raise ValueError(f"Reaction line needs 'equation : key : conversion' -> {line}")
        eq, key, conv = parts
        if "->" not in eq:
            raise ValueError(f"Reaction needs '->': {eq}")
        lhs, rhs = eq.split("->")
        nu = np.zeros(len(ctx.names))
        for side, sign in ((lhs, -1), (rhs, 1)):
            for term in side.split("+"):
                term = term.strip()
                m = re.match(r"^([0-9]*\.?[0-9]+)?\s*(.+)$", term)
                coef = float(m.group(1)) if m.group(1) else 1.0
                nu[_species_index(ctx, m.group(2))] += sign * coef
        k = _species_index(ctx, key)
        if nu[k] >= 0:
            raise ValueError(f"Key component {key} must be a reactant")
        rxns.append({"eq": eq.strip(), "nu": nu, "key": k, "conv": float(conv)})
    return rxns


@register
class ConversionReactor(Unit):
    TYPE, CATEGORY, LABEL = "ConversionReactor", "Fluid", "Conversion reactor"
    DESCRIPTION = ("Stoichiometric reactor. Reactions are applied in sequence, each to a fractional "
                   "conversion of its key reactant. Heat of reaction comes from formation enthalpies.")
    PARAMS = [P_("reactions", "text", "CO + 2 H2 -> CH3OH : CO : 0.3", "",
                 "One per line: 'equation : key reactant : conversion'"),
              P_("mode", "select", "isothermal", "", "Energy specification", ["isothermal", "outlet temperature", "adiabatic"]),
              P_("T_out", "number", 250.0, "C", "Outlet temperature", when={"mode": "outlet temperature"}),
              P_("dP", "number", 0.0, "bar", "Pressure drop")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        F = f.F.copy()
        extents = {}
        for r in parse_reactions(self.p["reactions"], ctx):
            xi = r["conv"] * F[r["key"]] / -r["nu"][r["key"]]
            F = F + xi * r["nu"]
            if np.any(F < -1e-9):
                raise ValueError(f"Reaction '{r['eq']}' consumes more of a reactant than is available")
            F = np.maximum(F, 0.0)
            extents[r["eq"]] = xi
        out = Stream(ctx, f.T, f.P - (self.num("dP") or 0) * 1e5, F, f.S, f"{self.id}.out")
        mode = self.p["mode"]
        if mode == "isothermal":
            out.flash_TP(f.T)
        elif mode == "outlet temperature":
            out.flash_TP(self.num("T_out") + K0)
        else:
            out.flash_PH(f.H)
        Q = out.H - f.H
        self.results = {"Duty [kW]": Q / 1e3, "Outlet T [C]": out.T - K0,
                        **{f"Extent: {k} [mol/s]": v for k, v in extents.items()}}
        return {"out": out}


@register
class ShortcutColumn(Unit):
    TYPE, CATEGORY, LABEL = "ShortcutColumn", "Fluid", "Distillation (shortcut)"
    DESCRIPTION = ("Fenske-Underwood-Gilliland column with Kirkbride feed location and a total condenser. "
                   "Non-key components are distributed by Fenske; volatilities are re-evaluated at "
                   "top and bottom conditions.")
    OUTPORTS = ["distillate", "bottoms"]
    PARAMS = [P_("LK", "component", "", "", "Light key"),
              P_("HK", "component", "", "", "Heavy key"),
              P_("LK_rec", "number", 0.99, "-", "Light-key recovery in distillate"),
              P_("HK_rec", "number", 0.99, "-", "Heavy-key recovery in bottoms"),
              P_("R_factor", "number", 1.3, "-", "Reflux ratio / minimum reflux"),
              P_("P", "number", 1.01325, "bar", "Column pressure")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        pkg = ctx.pkg
        if pkg is None:
            raise ValueError("Column needs volatile components")
        P = self.num("P") * 1e5
        vol = ctx.vol
        names = [ctx.names[i] for i in vol]
        lk, hk = names.index(self.p["LK"]), names.index(self.p["HK"])
        Fv = f.F[vol].copy()
        Ft = Fv.sum()
        z = Fv / Ft
        rl, rh = self.num("LK_rec"), self.num("HK_rec")

        # feed thermal condition q at column pressure
        Tb_f, Td_f = pkg.bubble_T(z, P), pkg.dew_T(z, P)
        hL = pkg.h(Tb_f, P, z, "L")
        hV = pkg.h(Td_f, P, z, "V")
        hF = (f.H - f._h_nonvol(f.T)) / Ft
        q = (hV - hF) / (hV - hL)

        def Keq(T, comp):
            fr = pkg.flash_TP(T, P, comp)
            if fr.phase == "VL":
                return fr.y / np.maximum(fr.x, 1e-300)
            return pkg.wilson_K(T, P)

        # volatilities, iterated between top and bottom conditions
        Tmid = 0.5 * (Tb_f + Td_f)
        alpha = Keq(Tmid, z)
        alpha = alpha / alpha[hk]
        for _ in range(6):
            Nmin = math.log((rl / (1 - rl)) * (rh / (1 - rh))) / math.log(alpha[lk])
            dhk_bhk = (1 - rh) / rh
            ratio = dhk_bhk * alpha ** Nmin  # d_i / b_i
            ratio[lk], ratio[hk] = rl / (1 - rl), dhk_bhk
            d = Fv * ratio / (1 + ratio)
            b = Fv - d
            xD, xB = d / d.sum(), b / b.sum()
            Ttop, Tbot = pkg.dew_T(xD, P), pkg.bubble_T(xB, P)
            Kt = Keq(Ttop - 0.01, xD)
            Kb = Keq(Tbot + 0.01, xB)
            a_new = np.sqrt((Kt / Kt[hk]) * (Kb / Kb[hk]))
            if np.max(np.abs(a_new / alpha - 1)) < 1e-6:
                alpha = a_new
                break
            alpha = a_new
        if alpha[lk] <= 1.0:
            raise ValueError("Light key is not more volatile than heavy key at column conditions")
        Nmin = math.log((rl / (1 - rl)) * (rh / (1 - rh))) / math.log(alpha[lk])

        # Underwood
        def und(th):
            return float(np.sum(alpha * z / (alpha - th))) - (1 - q)
        theta = brent(und, 1.0 + 1e-9, alpha[lk] - 1e-9)
        Rmin = float(np.sum(alpha * xD / (alpha - theta))) - 1
        Rmin = max(Rmin, 1e-6)
        R = self.num("R_factor") * Rmin
        X = (R - Rmin) / (R + 1)
        Y = 1 - math.exp((1 + 54.4 * X) / (11 + 117.2 * X) * (X - 1) / math.sqrt(X)) if X > 0 else 0.0
        N = (Y + Nmin) / (1 - Y)
        D, B = d.sum(), b.sum()
        kirk = ((z[hk] / z[lk]) * (xB[lk] / xD[hk]) ** 2 * (B / D)) ** 0.206
        NR = N / (1 + kirk)

        # products
        Tdb = pkg.bubble_T(xD, P)
        dist = Stream(ctx, Tdb, P, name=f"{self.id}.distillate")
        dist.F[vol] = d
        dist.flash_TP(Tdb, P)
        bott = Stream(ctx, Tbot, P, name=f"{self.id}.bottoms")
        bott.F = f.F.copy()
        bott.F[vol] = b
        bott.S = f.S.copy()
        bott.flash_TP(Tbot, P)
        Vtop = (R + 1) * D
        Qc = Vtop * (pkg.h(Ttop, P, xD, "V") - pkg.h(Tdb, P, xD, "L"))
        Qr = dist.H + bott.H + Qc - f.H
        self.results = {"Minimum stages (Fenske)": Nmin, "Minimum reflux (Underwood)": Rmin,
                        "Reflux ratio": R, "Theoretical stages": N,
                        "Feed stage from top (Kirkbride)": NR, "Feed q": q,
                        "Condenser duty [kW]": -Qc / 1e3, "Reboiler duty [kW]": Qr / 1e3,
                        "Top T [C]": Ttop - K0, "Bottom T [C]": Tbot - K0,
                        "LK/HK relative volatility": float(alpha[lk])}
        return {"distillate": dist, "bottoms": bott}


# ---------------------------------------------------------------------------
# Mineral processing units
# ---------------------------------------------------------------------------
def _recombine(S, frac_to_a):
    """Split solids matrix by per-(mineral, class) fraction."""
    a = S * frac_to_a
    return a, S - a


@register
class Comminution(Unit):
    TYPE, CATEGORY, LABEL = "Comminution", "Minerals", "Crusher / mill"
    DESCRIPTION = ("Size reduction to a target P80. Product size follows a Rosin-Rammler curve "
                   "(never coarser than the feed); specific energy from Bond's law.")
    PARAMS = [P_("P80", "number", 150.0, "um", "Target product P80"),
              P_("spread", "number", 1.0, "-", "Rosin-Rammler spread n of the product"),
              P_("work_index", "number", "", "kWh/t", "Bond work index (blank = mass-weighted mineral values)")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        out = f.copy(f"{self.id}.out")
        tot = f.S.sum(axis=0)
        if tot.sum() <= 0:
            raise ValueError("Comminution feed has no solids")
        F80 = p_size(tot)
        target = rrsb_fractions(self.num("P80"), self.num("spread"))
        for i in range(len(ctx.minerals)):
            m = f.S[i]
            if m.sum() <= 0:
                continue
            pas_f = passing_at_sieves(m)
            pas_p = passing_at_sieves(target)
            pas = np.maximum(pas_f, pas_p)
            frac = np.empty(NBINS)
            frac[0] = 1 - pas[0]
            frac[1:-1] = pas[:-1] - pas[1:]
            frac[-1] = pas[-1]
            out.S[i] = m.sum() * np.clip(frac, 0, None)
        P80 = p_size(out.S.sum(axis=0))
        wi = self.num("work_index")
        if wi is None:
            ms = f.S.sum(axis=1)
            wi = float(sum(ms[i] * ctx.mins[i].work_index for i in range(len(ms))) / ms.sum())
        E = 10 * wi * (1 / math.sqrt(P80) - 1 / math.sqrt(F80)) if F80 > P80 else 0.0  # kWh/t
        tph = tot.sum() * 3.6
        out.flash_TP()
        self.results = {"Feed F80 [um]": F80, "Product P80 [um]": P80, "Work index [kWh/t]": wi,
                        "Specific energy [kWh/t]": E, "Solids [t/h]": tph, "Power [kW]": E * tph}
        return {"out": out}


@register
class Classifier(Unit):
    TYPE, CATEGORY, LABEL = "Classifier", "Minerals", "Hydrocyclone / screen"
    DESCRIPTION = ("Size classification with a Whiten partition curve: coarse recovery = Rf + (1-Rf)E(d/d50c). "
                   "Use a sharpness of 2-4 for cyclones, 8-15 for screens.")
    OUTPORTS = ["coarse", "fine"]
    PARAMS = [P_("d50c", "number", 120.0, "um", "Corrected cut size"),
              P_("sharpness", "number", 3.0, "-", "Whiten sharpness alpha"),
              P_("Rf", "number", 0.3, "-", "Fraction of fluid (and fines bypass) to coarse product")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        a, d50, Rf = self.num("sharpness"), self.num("d50c"), min(max(self.num("Rf"), 0), 1)
        x = MID_UM / d50
        ax = np.minimum(a * x, 600.0)
        E = (np.exp(ax) - 1) / (np.exp(ax) + math.exp(min(a, 600.0)) - 2)
        E = np.clip(np.nan_to_num(E, nan=1.0), 0, 1)
        Ec = Rf + (1 - Rf) * E
        sc, sf = _recombine(f.S, Ec[None, :])
        c = Stream(ctx, f.T, f.P, f.F * Rf, sc, f"{self.id}.coarse").flash_TP()
        fi = Stream(ctx, f.T, f.P, f.F * (1 - Rf), sf, f"{self.id}.fine").flash_TP()
        self.results = {"Coarse solids [t/h]": c.solids_mass * 3.6, "Fine solids [t/h]": fi.solids_mass * 3.6,
                        "Fine P80 [um]": fi.psd_p80(), "Coarse P80 [um]": c.psd_p80()}
        return {"coarse": c, "fine": fi}


DEFAULT_KINETICS = """chalcopyrite: k=1.2, Rmax=0.95
sphalerite: k=0.9, Rmax=0.92
galena: k=1.4, Rmax=0.95
pentlandite: k=0.8, Rmax=0.90
pyrite: k=0.25, Rmax=0.60
quartz: k=0.02, Rmax=0.10"""


@register
class FlotationBank(Unit):
    TYPE, CATEGORY, LABEL = "FlotationBank", "Minerals", "Flotation bank"
    DESCRIPTION = ("Bank of N perfectly mixed cells with first-order kinetics per mineral: "
                   "R = Rmax[1-(1+k tau/N)^-N] times a size factor, plus entrainment of fines with water.")
    OUTPORTS = ["concentrate", "tailings"]
    PARAMS = [P_("kinetics", "text", DEFAULT_KINETICS, "", "Lines 'mineral: k=rate [1/min], Rmax=fraction'"),
              P_("tau", "number", 12.0, "min", "Total residence time"),
              P_("cells", "number", 6, "-", "Number of cells in series"),
              P_("water_rec", "number", 0.15, "-", "Fraction of water (fluid) reporting to concentrate"),
              P_("entrainment", "number", 0.6, "-", "Entrainment factor for particles below 40 um"),
              P_("size_effect", "select", "on", "", "Reduce floatability of very fine and coarse particles", ["on", "off"])]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        tab = parse_table(self.p["kinetics"], ["k", "Rmax"])
        tau, N = self.num("tau"), max(int(self.num("cells")), 1)
        Rw = min(max(self.num("water_rec"), 0), 1)
        ent = self.num("entrainment")
        if self.p["size_effect"] == "on":
            sf = np.exp(-0.5 * (np.log(MID_UM / 50.0) / 1.7) ** 2)
        else:
            sf = np.ones(NBINS)
        entf = np.where(MID_UM < 40.0, ent * Rw * (1 - MID_UM / 40.0) + 0.0, 0.0)
        frac = np.zeros_like(f.S)
        recs = {}
        for i, m in enumerate(ctx.minerals):
            kv = tab.get(m, {"k": 0.0, "Rmax": 0.0})
            k, Rmax = kv.get("k", 0.0), kv.get("Rmax", 1.0)
            Rk = Rmax * (1 - (1 + k * tau / N) ** (-N)) * sf
            frac[i] = np.clip(Rk + (1 - Rk) * entf, 0, 1)
        sc, st = _recombine(f.S, frac)
        conc = Stream(ctx, f.T, f.P, f.F * Rw, sc, f"{self.id}.concentrate").flash_TP()
        tail = Stream(ctx, f.T, f.P, f.F * (1 - Rw), st, f"{self.id}.tailings").flash_TP()
        cm = conc.solids_mass
        for i, m in enumerate(ctx.minerals):
            fin = f.S[i].sum()
            if fin > 0:
                recs[f"Recovery {m} [%]"] = 100 * sc[i].sum() / fin
                recs[f"Grade {m} in conc [% w]"] = 100 * sc[i].sum() / cm if cm > 0 else 0.0
        self.results = {"Concentrate solids [t/h]": cm * 3.6, "Mass pull [%]": 100 * cm / max(f.solids_mass, 1e-30),
                        **recs}
        return {"concentrate": conc, "tailings": tail}


@register
class LeachTank(Unit):
    TYPE, CATEGORY, LABEL = "LeachTank", "Minerals", "Leach reactor"
    DESCRIPTION = ("Dissolves a specified fraction of each mineral into the liquor as metal ions. "
                   "Non-metal mass is carried as 'dissolved other' so the mass balance closes. "
                   "Operates at a set temperature; the duty is reported.")
    PARAMS = [P_("extraction", "text", "chalcopyrite: 0.90\npyrite: 0.10", "", "Lines 'mineral: fraction dissolved'"),
              P_("T", "number", 80.0, "C", "Operating temperature"),
              P_("tau", "number", 240.0, "min", "Residence time (reported only)")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        tab = parse_table(self.p["extraction"], ["x"])
        out = f.copy(f"{self.id}.out")
        dissolved = {}
        for i, m in enumerate(ctx.minerals):
            v = tab.get(m)
            if not v:
                continue
            x = min(max(next(iter(v.values())), 0), 1)
            mass = f.S[i].sum() * x
            if mass <= 0:
                continue
            out.S[i] = f.S[i] * (1 - x)
            mol = mass / (ctx.mins[i].mw / 1000)
            metal_mass = 0.0
            for sp, nu in ctx.mins[i].leach_products.items():
                j = ctx.index(sp)
                out.F[j] += mol * nu
                metal_mass += mol * nu * ctx.mw[j]
            j = ctx.index("dissolved other")
            out.F[j] += (mass - metal_mass) / ctx.mw[j]
            dissolved[m] = 100 * x
        out.flash_TP(self.num("T") + K0)
        Q = out.H - f.H
        res = {"Duty [kW]": Q / 1e3, "Residence time [min]": self.num("tau")}
        for m, x in dissolved.items():
            res[f"Extraction {m} [%]"] = x
        for sp in ("Cu2+", "Zn2+", "Ni2+", "Li+", "Pb2+"):
            if sp in ctx.names:
                j = ctx.index(sp)
                dm = (out.F[j] - f.F[j]) * ctx.mw[j]
                if dm > 0:
                    res[f"{sp} dissolved [t/h]"] = dm * 3.6
        self.results = res
        return {"out": out}


@register
class Thickener(Unit):
    TYPE, CATEGORY, LABEL = "Thickener", "Minerals", "Thickener / filter"
    DESCRIPTION = "Solid-liquid separation to a target underflow solids content."
    OUTPORTS = ["underflow", "overflow"]
    PARAMS = [P_("solids_rec", "number", 0.995, "-", "Fraction of solids to underflow"),
              P_("uf_solids", "number", 60.0, "% w", "Underflow solids content")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        r = min(max(self.num("solids_rec"), 0), 1)
        c = min(max(self.num("uf_solids") / 100, 1e-3), 0.999)
        Su = f.S * r
        need = Su.sum() * (1 - c) / c
        fl = f.fluid_mass
        a = min(need / fl, 1.0) if fl > 0 else 0.0
        u = Stream(ctx, f.T, f.P, f.F * a, Su, f"{self.id}.underflow").flash_TP()
        o = Stream(ctx, f.T, f.P, f.F * (1 - a), f.S - Su, f"{self.id}.overflow").flash_TP()
        self.results = {"Underflow solids [% w]": u.summary()["solids_pct_w"] if ctx.minerals else 0,
                        "Liquor to overflow [%]": 100 * (1 - a)}
        return {"underflow": u, "overflow": o}


# ---------------------------------------------------------------------------
# Hydrometallurgy and refining
# ---------------------------------------------------------------------------
def tridiag_solve(lo, di, up, rhs):
    """Thomas algorithm: lo[i] multiplies x[i-1], up[i] multiplies x[i+1]."""
    n = len(di)
    c = [0.0] * n
    d = [0.0] * n
    c[0] = up[0] / di[0]
    d[0] = rhs[0] / di[0]
    for i in range(1, n):
        m = di[i] - lo[i] * c[i - 1]
        c[i] = up[i] / m if i < n - 1 else 0.0
        d[i] = (rhs[i] - lo[i] * d[i - 1]) / m
    x = [0.0] * n
    x[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x


DEFAULT_BETA = """La3+: 1
Ce3+: 6
Pr3+: 12
Nd3+: 17.4
Sm+HREE3+: 150
Th4+: 2000
Ca2+: 0.002"""


@register
class SXCascade(Unit):
    TYPE, CATEGORY, LABEL = "SXCascade", "Hydrometallurgy", "Solvent extraction cascade"
    DESCRIPTION = ("Counter-current mixer-settler bank with an extraction section and a scrub section, solved "
                   "stage by stage with constant distribution ratios D = D(key) x beta / beta(key). Typical "
                   "P507 (PC88A) chloride-system separation factors are the default. Loaded organic is fully "
                   "stripped; organic, scrub and strip acid circulate inside the unit and are not tracked.")
    OUTPORTS = ["extract", "raffinate"]
    PARAMS = [P_("betas", "text", DEFAULT_BETA, "", "Separation factors 'species: beta' (unlisted species stay in the aqueous)"),
              P_("key", "solute", "Nd3+", "", "Key species for the distribution ratio"),
              P_("D_key", "number", 1.0, "-", "Distribution ratio of the key species (organic/aqueous)"),
              P_("N_ext", "number", 30, "-", "Extraction stages"),
              P_("N_scrub", "number", 30, "-", "Scrub stages"),
              P_("OA", "number", 1.0, "-", "Organic flow / aqueous feed flow (by volume)"),
              P_("SA", "number", 0.5, "-", "Scrub acid flow / aqueous feed flow (by volume)"),
              P_("strip_water", "number", 0.3, "-", "Strip liquor volume / aqueous feed volume")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        beta = {k: next(iter(v.values()), 0.0) if v else 0.0 for k, v in parse_table(self.p["betas"], ["beta"]).items()}
        key = self.p["key"]
        if key not in beta or beta[key] <= 0:
            raise ValueError(f"Key species '{key}' needs a separation factor in the list")
        Dk, Ne, Ns = self.num("D_key"), max(int(self.num("N_ext")), 1), max(int(self.num("N_scrub")), 0)
        OA, SA, sw = self.num("OA"), self.num("SA"), min(max(self.num("strip_water"), 0.0), 1.0)
        N = Ne + Ns
        Lf = 1.0 + SA          # aqueous flow in the extraction section (feed basis)
        Ls = SA if SA > 0 else 1e-12  # aqueous flow in the scrub section
        ext = Stream(ctx, f.T, f.P, name=f"{self.id}.extract")
        raf = Stream(ctx, f.T, f.P, name=f"{self.id}.raffinate")
        raf.S = f.S.copy()
        rec = {}
        for i, n in enumerate(ctx.names):
            Fi = f.F[i]
            if ctx.comps[i].solute and n in beta and Fi > 0:
                D = Dk * beta[n] / beta[key]
                e = [D * OA / (Lf if k < Ne else Ls) for k in range(N)]
                lo = [0.0] + [e[k - 1] for k in range(1, N)]
                di = [-(1.0 + e[k]) for k in range(N)]
                up = [1.0] * (N - 1) + [0.0]
                rhs = [0.0] * N
                rhs[Ne - 1] = -Fi
                a = tridiag_solve(lo, di, up, rhs)
                to_ext = min(max(e[N - 1] * a[N - 1], 0.0), Fi)
                ext.F[i] = to_ext
                raf.F[i] = Fi - to_ext
                rec[n] = to_ext / Fi
            elif ctx.comps[i].solute or i not in ctx.vol:
                raf.F[i] = Fi
            else:  # volatile (water): strip liquor takes a share of the feed water
                ext.F[i] = Fi * sw
                raf.F[i] = Fi * (1 - sw)
        ext.flash_TP()
        raf.flash_TP()
        res = {"Stages (extraction + scrub)": N}
        re_ions = [n for n in rec if n not in ("Ca2+", "Th4+")]
        tot = sum(ext.F[ctx.index(n)] for n in re_ions)
        for n, r in rec.items():
            res[f"{n} to extract [%]"] = 100 * r
        if tot > 0:
            res[f"{key} share of RE in extract [%]"] = 100 * ext.F[ctx.index(key)] / tot
        self.results = res
        return {"extract": ext, "raffinate": raf}


def parse_precip(text, ctx):
    rx = []
    for line in (text or "").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.split(":")]
        if len(parts) != 2 or "->" not in parts[0]:
            raise ValueError(f"Precipitation line needs 'ion -> solid : fraction' -> {line}")
        ion, solid = [x.strip() for x in parts[0].split("->")]
        rx.append((ctx.index(ion), ctx.mindex(solid), float(parts[1])))
    return rx


@register
class Precipitator(Unit):
    TYPE, CATEGORY, LABEL = "Precipitator", "Hydrometallurgy", "Precipitation and filtration"
    DESCRIPTION = ("Precipitates a fraction of each dissolved ion as a named solid, then filters. The anion or "
                   "oxygen mass of the solid is drawn from the 'reagent' species, so feed a reagent stream to the "
                   "inlet. Oxalate precipitation followed by calcination can be represented directly as the oxide.")
    OUTPORTS = ["solids", "liquor"]
    PARAMS = [P_("reactions", "text", "Nd3+ -> Nd2O3 : 0.995", "", "Lines 'ion -> solid : fraction precipitated'"),
              P_("T", "number", 60.0, "C", "Operating temperature"),
              P_("moisture", "number", 15.0, "% w", "Liquor content of the filter cake"),
              P_("p80", "number", 20.0, "um", "Precipitate particle size P80")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        out = f.copy()
        jr = ctx.index("reagent") if "reagent" in ctx.names else None
        used = 0.0
        frac = rrsb_fractions(self.num("p80"), 1.2)
        made = {}
        for ion, sol, x in parse_precip(self.p["reactions"], ctx):
            x = min(max(x, 0.0), 1.0)
            m = ctx.mins[sol]
            nu = m.leach_products.get(ctx.names[ion])
            if not nu:
                raise ValueError(f"{m.name} does not contain {ctx.names[ion]}")
            mol_ion = out.F[ion] * x
            mol_solid = mol_ion / nu
            need = mol_solid * m.mw / 1000 - mol_ion * ctx.mw[ion]  # kg/s of anion/oxygen
            out.F[ion] -= mol_ion
            out.S[sol] += mol_solid * m.mw / 1000 * frac
            used += need
            made[m.name] = made.get(m.name, 0.0) + mol_solid * m.mw / 1000
        if used > 1e-12:
            if jr is None:
                raise ValueError("Add the 'reagent' species and a reagent feed stream to supply the precipitant")
            have = out.F[jr] * ctx.mw[jr]
            if used > have * (1 + 1e-9):
                raise ValueError(f"Not enough reagent: needs {used * 3.6:.4g} t/h, inlet carries {have * 3.6:.4g} t/h")
            out.F[jr] = max(out.F[jr] - used / ctx.mw[jr], 0.0)
        mc = min(max(self.num("moisture") / 100, 0.0), 0.99)
        ms = out.solids_mass
        fl = out.fluid_mass
        a = min(ms * mc / (1 - mc) / fl, 1.0) if fl > 0 else 0.0
        T = self.num("T") + K0
        cake = Stream(ctx, T, f.P, out.F * a, out.S, f"{self.id}.solids").flash_TP()
        liq = Stream(ctx, T, f.P, out.F * (1 - a), np.zeros_like(out.S), f"{self.id}.liquor").flash_TP()
        res = {"Duty [kW]": (cake.H + liq.H - f.H) / 1e3, "Reagent consumed [t/h]": used * 3.6}
        for k, v in made.items():
            res[f"{k} precipitated [t/h]"] = v * 3.6
        self.results = res
        return {"solids": cake, "liquor": liq}


@register
class Electrolysis(Unit):
    TYPE, CATEGORY, LABEL = "Electrolysis", "Hydrometallurgy", "Molten-salt electrolysis"
    DESCRIPTION = ("Reduces oxide feed to metal in a fluoride melt. Oxygen leaves as anode gas (reported as O2; "
                   "with carbon anodes it is really CO and CO2). Electrical energy from a specific consumption.")
    OUTPORTS = ["metal", "offgas"]
    PARAMS = [P_("reactions", "text", "Nd2O3 -> Nd metal : 0.98", "", "Lines 'oxide -> metal : fraction reduced'"),
              P_("T", "number", 1050.0, "C", "Cell temperature"),
              P_("energy", "number", 10.0, "kWh/kg", "Electrical energy per kg metal"),
              P_("current_eff", "number", 0.75, "-", "Current efficiency (reported)")]

    def run(self, inlets, ctx):
        f = self.feed(inlets, ctx)
        out = f.copy()
        if "oxygen" not in ctx.names:
            raise ValueError("Add 'oxygen' to the components to carry the anode gas")
        jo = ctx.index("oxygen")
        metal_kg = 0.0
        for line in (self.p["reactions"] or "").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            eq, x = [p.strip() for p in line.split(":")]
            ox, me = [ctx.mindex(t.strip()) for t in eq.split("->")]
            x = min(max(float(x), 0.0), 1.0)
            mo, mm = ctx.mins[ox], ctx.mins[me]
            ion = next(iter(mm.leach_products))
            n_ox, n_me = mo.leach_products.get(ion), mm.leach_products.get(ion)
            if not n_ox or not n_me:
                raise ValueError(f"{mo.name} and {mm.name} must share a metal")
            mass_ox = out.S[ox] * x
            mol_metal = mass_ox.sum() / (mo.mw / 1000) * n_ox
            m_metal = mol_metal / n_me * mm.mw / 1000
            out.S[ox] -= mass_ox
            tot = mass_ox.sum()
            out.S[me] += (mass_ox / tot * m_metal) if tot > 0 else 0.0
            out.F[jo] += (tot - m_metal) / ctx.mw[jo]
            metal_kg += m_metal
        out.flash_TP(self.num("T") + K0)
        V, L = out.split_phases()
        L.name, V.name = f"{self.id}.metal", f"{self.id}.offgas"
        kw = self.num("energy") * metal_kg * 3600
        self.results = {"Metal produced [t/h]": metal_kg * 3.6, "Electrical power [kW]": kw,
                        "Current efficiency": self.num("current_eff"), "Thermal duty [kW]": (out.H - f.H) / 1e3}
        return {"metal": L, "offgas": V}


def catalogue():
    return [cls.schema() for cls in UNIT_TYPES.values()]
