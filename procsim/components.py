"""Pure-component and mineral property data.

Fluid data: critical constants and acentric factors (Poling, Prausnitz & O'Connell,
"The Properties of Gases and Liquids"), ideal-gas heat capacity polynomials
Cp = A + B T + C T^2 + D T^3 [J/mol/K, T in K] (Reid, Prausnitz & Poling, 4th ed.)
and ideal-gas enthalpies of formation at 298.15 K [kJ/mol].

Solutes are non-volatile dissolved species produced by leaching. They never enter
the vapour phase and are carried in the liquid with a constant heat capacity.

Minerals are solids tracked by mass and particle-size class.
"""
from dataclasses import dataclass, field

R = 8.314462618  # J/mol/K
T_REF = 298.15


@dataclass(frozen=True)
class Component:
    name: str
    formula: str
    mw: float          # g/mol
    tc: float = 0.0    # K
    pc: float = 0.0    # Pa
    omega: float = 0.0
    cp: tuple = (0.0, 0.0, 0.0, 0.0)  # ideal gas Cp coefficients
    hf: float = 0.0    # J/mol, ideal gas at 298.15 K
    solute: bool = False
    cp_liq: float = 75.0  # J/mol/K, solutes only

    def cp_ig(self, T):
        A, B, C, D = self.cp
        return A + B * T + C * T ** 2 + D * T ** 3

    def h_ig(self, T):
        """Ideal-gas enthalpy relative to elements at 298.15 K [J/mol]."""
        A, B, C, D = self.cp
        def F(t):
            return A * t + B * t ** 2 / 2 + C * t ** 3 / 3 + D * t ** 4 / 4
        return self.hf + F(T) - F(T_REF)


def _c(name, formula, mw, tc, pc_bar, omega, cp, hf_kj):
    return Component(name, formula, mw, tc, pc_bar * 1e5, omega, cp, hf_kj * 1e3)


COMPONENTS = {c.name: c for c in [
    _c("methane", "CH4", 16.043, 190.56, 45.99, 0.011, (19.25, 5.213e-2, 1.197e-5, -1.132e-8), -74.52),
    _c("ethane", "C2H6", 30.070, 305.32, 48.72, 0.099, (5.409, 1.781e-1, -6.938e-5, 8.713e-9), -83.82),
    _c("propane", "C3H8", 44.097, 369.83, 42.48, 0.152, (-4.224, 3.063e-1, -1.586e-4, 3.215e-8), -104.68),
    _c("n-butane", "C4H10", 58.123, 425.12, 37.96, 0.200, (9.487, 3.313e-1, -1.108e-4, -2.822e-9), -125.79),
    _c("n-pentane", "C5H12", 72.150, 469.70, 33.70, 0.252, (-3.626, 4.873e-1, -2.580e-4, 5.305e-8), -146.76),
    _c("n-hexane", "C6H14", 86.177, 507.60, 30.25, 0.301, (-4.413, 5.820e-1, -3.119e-4, 6.494e-8), -166.92),
    _c("n-heptane", "C7H16", 100.204, 540.20, 27.40, 0.350, (-5.146, 6.762e-1, -3.651e-4, 7.658e-8), -187.78),
    _c("benzene", "C6H6", 78.114, 562.05, 48.95, 0.210, (-33.92, 4.739e-1, -3.017e-4, 7.130e-8), 82.88),
    _c("toluene", "C7H8", 92.141, 591.75, 41.08, 0.264, (-24.35, 5.125e-1, -2.765e-4, 4.911e-8), 50.17),
    _c("water", "H2O", 18.015, 647.14, 220.64, 0.344, (32.24, 1.924e-3, 1.055e-5, -3.596e-9), -241.81),
    _c("nitrogen", "N2", 28.014, 126.20, 33.98, 0.037, (31.15, -1.357e-2, 2.680e-5, -1.168e-8), 0.0),
    _c("oxygen", "O2", 31.999, 154.58, 50.43, 0.022, (28.11, -3.680e-6, 1.746e-5, -1.065e-8), 0.0),
    _c("hydrogen", "H2", 2.016, 33.19, 13.13, -0.216, (27.14, 9.274e-3, -1.381e-5, 7.645e-9), 0.0),
    _c("carbon monoxide", "CO", 28.010, 132.85, 34.94, 0.045, (30.87, -1.285e-2, 2.789e-5, -1.272e-8), -110.53),
    _c("carbon dioxide", "CO2", 44.010, 304.21, 73.83, 0.224, (19.80, 7.344e-2, -5.602e-5, 1.715e-8), -393.51),
    _c("hydrogen sulfide", "H2S", 34.080, 373.40, 89.63, 0.090, (31.94, 1.436e-3, 2.432e-5, -1.176e-8), -20.63),
    _c("sulfur dioxide", "SO2", 64.064, 430.80, 78.84, 0.245, (23.85, 6.699e-2, -4.961e-5, 1.328e-8), -296.81),
    _c("ammonia", "NH3", 17.031, 405.40, 113.53, 0.257, (27.31, 2.383e-2, 1.707e-5, -1.185e-8), -45.90),
    _c("methanol", "CH3OH", 32.042, 512.64, 80.97, 0.565, (21.15, 7.092e-2, 2.587e-5, -2.852e-8), -200.94),
    _c("ethanol", "C2H5OH", 46.069, 513.92, 61.48, 0.649, (9.014, 2.141e-1, -8.390e-5, 1.373e-9), -234.95),
]}

# Non-volatile dissolved species (leach liquor). Reference state: dissolved at 298.15 K, Hf = 0.
for _n, _f, _mw in [("Cu2+", "Cu2+", 63.546), ("Zn2+", "Zn2+", 65.38), ("Pb2+", "Pb2+", 207.2),
                    ("Fe2+", "Fe2+", 55.845), ("Ni2+", "Ni2+", 58.693), ("Li+", "Li+", 6.94),
                    ("dissolved other", "X", 100.0),
                    ("La3+", "La3+", 138.905), ("Ce3+", "Ce3+", 140.116), ("Pr3+", "Pr3+", 140.908),
                    ("Nd3+", "Nd3+", 144.242), ("Sm+HREE3+", "SmX3+", 150.36), ("Th4+", "Th4+", 232.038),
                    ("Ca2+", "Ca2+", 40.078),
                    # Generic precipitant / neutraliser anion mass carrier (carbonate-equivalent basis)
                    ("reagent", "R", 60.0)]:
    COMPONENTS[_n] = Component(_n, _f, _mw, solute=True)


@dataclass(frozen=True)
class Mineral:
    name: str
    formula: str
    mw: float          # g/mol
    density: float     # kg/m3
    cp: float          # J/kg/K
    leach_products: dict = field(default_factory=dict)  # solute -> mol per mol mineral
    work_index: float = 14.0  # Bond ball-mill work index, kWh/t (typical)


MINERALS = {m.name: m for m in [
    Mineral("chalcopyrite", "CuFeS2", 183.53, 4200, 540, {"Cu2+": 1}, 10.0),
    Mineral("pyrite", "FeS2", 119.98, 5010, 520, {"Fe2+": 1}, 9.0),
    Mineral("sphalerite", "ZnS", 97.44, 4050, 470, {"Zn2+": 1}, 11.0),
    Mineral("galena", "PbS", 239.27, 7500, 210, {"Pb2+": 1}, 10.0),
    Mineral("pentlandite", "Ni4.5Fe4.5S8", 771.9, 4800, 600, {"Ni2+": 4.5}, 12.0),
    Mineral("spodumene", "LiAlSi2O6", 186.09, 3150, 850, {"Li+": 1}, 13.0),
    Mineral("magnetite", "Fe3O4", 231.53, 5150, 650, {"Fe2+": 3}, 10.0),
    Mineral("hematite", "Fe2O3", 159.69, 5260, 650, {"Fe2+": 2}, 12.0),
    Mineral("quartz", "SiO2", 60.08, 2650, 740, {}, 16.0),
    # ---- rare earths ------------------------------------------------------------------
    # Bastnaesite (RE)CO3F. RE split from USGS Minerals Yearbook 2020, Table 2 (Mountain Pass, % of TREO),
    # converted to mole fractions; Sm to Lu and Y lumped as Sm+HREE3+. Trace Th (assumed).
    Mineral("bastnaesite", "(Ce,La,Nd,Pr)CO3F", 219.81, 4950, 600,
            {"La3+": 0.3514, "Ce3+": 0.4774, "Pr3+": 0.0415, "Nd3+": 0.1171, "Sm+HREE3+": 0.0126, "Th4+": 0.002}, 12.0),
    # Monazite (RE,Th)PO4. RE split from USGS Minerals Yearbook 2020, Table 2 (Mt Weld Central, % of TREO);
    # Th at 6 mol% is an assumption (no Th figure in the GSCO Data Library).
    Mineral("monazite", "(Ce,La,Nd,Th)PO4", 241.59, 5150, 560,
            {"La3+": 0.2323, "Ce3+": 0.4378, "Pr3+": 0.048, "Nd3+": 0.1708, "Sm+HREE3+": 0.0511, "Th4+": 0.06}, 13.0),
    # Monazite-xenotime concentrate. RE split from Iluka Eneabba refinery presentation 2024, p25 and p27
    # (La 21.8, Ce 45.0, Pr 4.6, Nd 16.6, Sm 3, Gd 1, Tb 0.2, Dy 0.9, Y 6 % of REO); Th assumed as above.
    Mineral("monazite-xenotime", "(Ce,La,Nd,Y,Th)PO4", 242.34, 5000, 560,
            {"La3+": 0.2087, "Ce3+": 0.4078, "Pr3+": 0.0421, "Nd3+": 0.1539, "Sm+HREE3+": 0.1276, "Th4+": 0.06}, 13.0),
    Mineral("apatite", "Ca5(PO4)3F", 504.3, 3190, 700, {"Ca2+": 5}, 9.0),
    Mineral("calcite", "CaCO3", 100.09, 2710, 820, {"Ca2+": 1}, 10.0),
    Mineral("barite", "BaSO4", 233.39, 4480, 460, {}, 6.0),
    # Precipitated intermediates and products (leachable back to ions in acid)
    Mineral("La2(CO3)3", "La2(CO3)3", 457.84, 2600, 600, {"La3+": 2}, 5.0),
    Mineral("Ce2(CO3)3", "Ce2(CO3)3", 460.26, 2600, 600, {"Ce3+": 2}, 5.0),
    Mineral("Pr2(CO3)3", "Pr2(CO3)3", 461.84, 2600, 600, {"Pr3+": 2}, 5.0),
    Mineral("Nd2(CO3)3", "Nd2(CO3)3", 468.51, 2600, 600, {"Nd3+": 2}, 5.0),
    Mineral("Sm2(CO3)3", "Sm2(CO3)3", 480.75, 2600, 600, {"Sm+HREE3+": 2}, 5.0),
    Mineral("La2O3", "La2O3", 325.81, 6510, 330, {"La3+": 2}, 8.0),
    Mineral("CeO2", "CeO2", 172.11, 7220, 360, {"Ce3+": 1}, 8.0),
    Mineral("Pr6O11", "Pr6O11", 1021.44, 6900, 330, {"Pr3+": 6}, 8.0),
    Mineral("Nd2O3", "Nd2O3", 336.48, 7240, 330, {"Nd3+": 2}, 8.0),
    Mineral("Sm2O3", "Sm2O3", 348.72, 8350, 330, {"Sm+HREE3+": 2}, 8.0),
    Mineral("ThO2", "ThO2", 264.04, 10000, 240, {"Th4+": 1}, 8.0),
    Mineral("Nd metal", "Nd", 144.242, 7010, 190, {"Nd3+": 1}, 8.0),
]}

# Default kij for Peng-Robinson (typical literature values; 0 otherwise).
_KIJ = {}
_HC = ["methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "n-heptane", "benzene", "toluene"]
for _h in _HC:
    _KIJ[frozenset(("water", _h))] = 0.48
    _KIJ[frozenset(("carbon dioxide", _h))] = 0.12
    _KIJ[frozenset(("nitrogen", _h))] = 0.08
    _KIJ[frozenset(("hydrogen sulfide", _h))] = 0.07
_KIJ[frozenset(("water", "carbon dioxide"))] = 0.12
_KIJ[frozenset(("water", "nitrogen"))] = 0.48


def kij(a, b):
    return _KIJ.get(frozenset((a, b)), 0.0)


def component(name):
    try:
        return COMPONENTS[name]
    except KeyError:
        raise KeyError(f"Unknown component '{name}'. Available: {sorted(COMPONENTS)}")


def mineral(name):
    try:
        return MINERALS[name]
    except KeyError:
        raise KeyError(f"Unknown mineral '{name}'. Available: {sorted(MINERALS)}")
