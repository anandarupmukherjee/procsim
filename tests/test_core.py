import json
import os
import math

import numpy as np
import pytest

from procsim import Flowsheet
from procsim.stream import Context, Stream, p_size, rrsb_fractions
from procsim.thermo import PengRobinson, Ideal

EX = os.path.join(os.path.dirname(__file__), "..", "procsim", "data", "examples")

# Normal boiling points [K] (NIST WebBook)
NBP = {"propane": 231.04, "n-butane": 272.66, "n-hexane": 341.88, "benzene": 353.24,
       "toluene": 383.75, "water": 373.12, "methanol": 337.69, "ammonia": 239.82}


@pytest.mark.parametrize("pkg", [PengRobinson, Ideal])
@pytest.mark.parametrize("name", list(NBP))
def test_normal_boiling_points(pkg, name):
    T = pkg([name]).bubble_T([1.0], 101325.0)
    assert abs(T - NBP[name]) < 3.0


def test_benzene_toluene_vle():
    t = PengRobinson(["benzene", "toluene"])
    assert abs(t.bubble_T([0.5, 0.5], 101325) - 365.2) < 1.0
    assert abs(t.dew_T([0.5, 0.5], 101325) - 371.9) < 1.0


def test_bubble_below_dew_for_wide_mixture():
    t = PengRobinson(["methane", "ethane", "propane", "n-butane", "n-hexane"])
    z = np.array([0.4, 0.2, 0.15, 0.15, 0.1])
    assert t.bubble_T(z, 20e5) < t.dew_T(z, 20e5)


def test_water_latent_heat():
    t = PengRobinson(["water"])
    dh = t.h(373.15, 101325, [1], "V") - t.h(373.15, 101325, [1], "L")
    assert abs(dh - 40660) / 40660 < 0.06


def test_heat_of_reaction_methanol():
    spec = {"components": ["carbon monoxide", "hydrogen", "methanol"], "package": "Ideal",
            "units": [{"id": "R", "type": "ConversionReactor",
                       "params": {"reactions": "CO + 2 H2 -> CH3OH : CO : 1.0", "mode": "isothermal"}}],
            "streams": [{"id": "f", "to": {"unit": "R", "port": "in"},
                         "feed": {"T": 600.0, "P": 1e4, "flows": {"carbon monoxide": 1, "hydrogen": 2}}},
                        {"id": "p", "from": {"unit": "R", "port": "out"}}]}
    fs = Flowsheet(spec).solve()
    q = fs.units["R"].results["Duty [kW]"]
    # Reaction enthalpy at 600 K is about -100.4 kJ/mol (gas phase)
    assert -103 < q < -97


def test_rrsb_p80():
    for p80 in (75, 150, 1000, 12000):
        assert abs(p_size(rrsb_fractions(p80, 1.0)) - p80) / p80 < 0.03


def test_mixer_conserves_energy():
    # Hydrocarbons: PR with kij = 0 is appropriate. (Polar mixtures such as water-ethanol need an
    # activity-coefficient model; PR with kij = 0 gives a large spurious heat of mixing.)
    ctx = Context(["n-hexane", "n-heptane"])
    a = Stream.from_spec(ctx, {"T": 290, "P": 2e5, "flows": {"n-hexane": 3}})
    b = Stream.from_spec(ctx, {"T": 340, "P": 1.5e5, "flows": {"n-heptane": 1}})
    from procsim.units import mix
    m = mix([a, b], ctx)
    assert abs(m.H - a.H - b.H) < 1e-3 * abs(a.H)
    assert m.P == 1.5e5
    assert 290 < m.T < 340


@pytest.mark.parametrize("name", ["ngl_fractionation", "methanol_loop", "copper_concentrator", "neodymium_bastnaesite",
                                  "neodymium_monazite_mtweld", "neodymium_eneabba_concentrate"])
def test_examples_close_balances(name):
    with open(os.path.join(EX, name + ".json")) as f:
        spec = json.load(f)
    fs = Flowsheet(spec).solve()
    b = fs.balances()
    assert abs(b["mass_closure_pct"]) < 1e-2
    assert abs(b["energy_closure_kW"]) < 1.0
    assert all(lp["converged"] for lp in fs.loops)


def test_column_meets_key_recoveries():
    with open(os.path.join(EX, "ngl_fractionation.json")) as f:
        fs = Flowsheet(json.load(f)).solve()
    ctx = fs.ctx
    i = ctx.index("propane")
    feed, dist = fs.stream("C3+"), fs.stream("LPG")
    assert abs(dist.F[i] / feed.F[i] - 0.98) < 1e-6


def test_comminution_hits_target_and_flotation_conserves():
    with open(os.path.join(EX, "copper_concentrator.json")) as f:
        fs = Flowsheet(json.load(f)).solve()
    assert abs(fs.units["BM"].results["Product P80 [um]"] - 250) / 250 < 0.05
    f_ = fs.stream("Cyclone overflow")
    c, t = fs.stream("Rougher conc"), fs.stream("Rougher tail")
    assert np.allclose(f_.S, c.S + t.S)
    assert np.allclose(f_.F, c.F + t.F)


# ---------------------------------------------------------------------------
# Rare-earth / hydrometallurgy
# ---------------------------------------------------------------------------
def _sx_feed(ctx, flows, water=100.0):
    return Stream.from_spec(ctx, {"T": 313.15, "P": 101325, "flows": dict(flows, water=water)})


def test_single_stage_sx_matches_analytic():
    from procsim.units import SXCascade
    ctx = Context(["water", "Nd3+", "Pr3+"], package="Ideal")
    f = _sx_feed(ctx, {"Nd3+": 1.0, "Pr3+": 1.0})
    u = SXCascade("S", {"betas": "Pr3+: 1\nNd3+: 1.5", "key": "Nd3+", "D_key": 2.0, "N_ext": 1, "N_scrub": 0, "OA": 1.0, "SA": 0.0})
    out = u.run({"in": [f]}, ctx)
    E = 2.0  # D * O/A with no scrub
    assert abs(out["extract"].F[ctx.index("Nd3+")] - E / (1 + E)) < 1e-9
    Ep = 2.0 / 1.5
    assert abs(out["extract"].F[ctx.index("Pr3+")] - Ep / (1 + Ep)) < 1e-9


def test_sx_cascade_matches_dense_solve():
    from procsim.units import tridiag_solve
    rng = np.random.default_rng(1)
    n = 40
    lo, up = rng.uniform(0.1, 2, n), rng.uniform(0.1, 2, n)
    di = -(lo + up + 1)
    lo[0] = up[-1] = 0
    rhs = rng.normal(size=n)
    A = np.diag(di) + np.diag(lo[1:], -1) + np.diag(up[:-1], 1)
    assert np.allclose(tridiag_solve(list(lo), list(di), list(up), list(rhs)), np.linalg.solve(A, rhs))


def test_precipitator_conserves_mass_and_needs_reagent():
    from procsim.units import Precipitator
    from procsim.flowsheet import FlowsheetError
    ctx = Context(["water", "Nd3+", "reagent"], ["Nd2O3"], package="Ideal")
    f = _sx_feed(ctx, {"Nd3+": 1.0, "reagent": 1.0})
    u = Precipitator("P", {"reactions": "Nd3+ -> Nd2O3 : 1.0"})
    out = u.run({"in": [f]}, ctx)
    assert abs(out["solids"].mass + out["liquor"].mass - f.mass) < 1e-12
    assert abs(out["solids"].S.sum() - 0.5 * 0.33648) < 1e-9
    f2 = _sx_feed(ctx, {"Nd3+": 1.0, "reagent": 0.01})
    with pytest.raises(ValueError, match="Not enough reagent"):
        u.run({"in": [f2]}, ctx)


ROUTE_WINDOWS = {
    # name: (feed stream, value mineral, Nd recovery window feed -> oxide)
    "neodymium_bastnaesite": ("Run-of-mine ore", "bastnaesite", (0.5, 0.7)),
    "neodymium_monazite_mtweld": ("Run-of-mine ore", "monazite", (0.45, 0.65)),
    "neodymium_eneabba_concentrate": ("Monazite-xenotime concentrate", "monazite-xenotime", (0.85, 0.95)),
}


@pytest.mark.parametrize("name", list(ROUTE_WINDOWS))
def test_neodymium_route_products(name):
    from procsim.components import MINERALS
    feed_sid, mineral, (lo, hi) = ROUTE_WINDOWS[name]
    with open(os.path.join(EX, name + ".json")) as f:
        fs = Flowsheet(json.load(f)).solve()
    ctx = fs.ctx
    assert abs(fs.balances()["mass_closure_pct"]) < 1e-3
    assert all(lp["converged"] for lp in fs.loops)
    nd = fs.stream("Nd strip liquor")
    re = sum(nd.F[ctx.index(n)] for n in ["La3+", "Ce3+", "Pr3+", "Nd3+", "Sm+HREE3+"])
    assert nd.F[ctx.index("Nd3+")] / re > 0.995
    m = MINERALS[mineral]
    nd_in = fs.stream(feed_sid).S[ctx.mindex(mineral)].sum() / (m.mw / 1000) * m.leach_products["Nd3+"] * 0.144242
    nd_out = fs.stream("Nd oxide").S[ctx.mindex("Nd2O3")].sum() * 2 * 144.242 / 336.48
    assert lo < nd_out / nd_in < hi
    assert fs.stream("Nd metal").S[ctx.mindex("Nd metal")].sum() > 0


def test_sourced_distributions_reproduce_published_oxide_split():
    """Mineral RE mole fractions must reproduce the cited % of TREO for La, Ce, Pr, Nd within 0.1 point."""
    from procsim.components import MINERALS
    ox = {"La3+": 162.905, "Ce3+": 172.115, "Pr3+": 170.24, "Nd3+": 168.24, "Sm+HREE3+": 174.36}
    cited = {"bastnaesite": (34.00, 48.80, 4.20, 11.70),        # USGS MYB 2020 Table 2, Mountain Pass
             "monazite": (23.88, 47.55, 5.16, 18.13),           # USGS MYB 2020 Table 2, Mt Weld Central
             "monazite-xenotime": (21.8, 45.0, 4.6, 16.6)}      # Iluka Eneabba refinery 2024 p25
    for m, ref in cited.items():
        lp = MINERALS[m].leach_products
        mass = {k: v * ox[k] for k, v in lp.items() if k in ox}
        tot = sum(mass.values())
        pct = [100 * mass[k] / tot for k in ("La3+", "Ce3+", "Pr3+", "Nd3+")]
        scale = sum(ref) / sum(pct)  # Sm-to-Lu and Y are lumped at the Sm oxide mass, so compare LREE shares
        assert all(abs(a * scale - b) < 0.1 for a, b in zip(pct, ref)), (m, pct, ref)
