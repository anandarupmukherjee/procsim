"""Generates the bundled example flowsheets (procsim/data/examples/*.json). Positions (x, y) are for the web editor.
Run from the repository root: python examples/build_examples.py"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # run from a checkout without installing
OUT = os.path.join(os.path.dirname(HERE), "procsim", "data", "examples")


def S(id, frm=None, to=None, feed=None, guess=None):
    d = {"id": id, "from": None, "to": None}
    if frm:
        d["from"] = {"unit": frm[0], "port": frm[1]}
    if to:
        d["to"] = {"unit": to[0], "port": to[1]}
    if feed:
        d["feed"] = feed
    if guess:
        d["guess"] = guess
    return d


def U(id, type, x, y, **params):
    return {"id": id, "type": type, "x": x, "y": y, "params": params}


C = 273.15

ngl = {
    "name": "NGL recovery and fractionation",
    "description": "Chill a rich gas, separate the condensate, then de-ethanise and de-propanise the liquids.",
    "components": ["nitrogen", "carbon dioxide", "methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane"],
    "package": "Peng-Robinson",
    "units": [
        U("E-100", "Heater", 200, 160, mode="temperature", T_out=-35.0, dP=0.5),
        U("V-100", "Flash", 380, 160, mode="adiabatic", P=59.5),
        U("V-101", "Valve", 540, 260, P_out=26.0),
        U("T-100", "ShortcutColumn", 700, 260, LK="ethane", HK="propane", LK_rec=0.98, HK_rec=0.98, R_factor=1.3, P=26.0),
        U("T-101", "ShortcutColumn", 880, 360, LK="propane", HK="n-butane", LK_rec=0.98, HK_rec=0.98, R_factor=1.3, P=15.0),
    ],
    "streams": [
        S("Rich gas", to=("E-100", "in"), feed={"T": 30 + C, "P": 60e5, "flows": {
            "nitrogen": 1.0, "carbon dioxide": 1.5, "methane": 70.0, "ethane": 12.0, "propane": 8.0,
            "n-butane": 4.0, "n-pentane": 2.0, "n-hexane": 1.5}}),
        S("Chilled gas", ("E-100", "out"), ("V-100", "in")),
        S("Residue gas", ("V-100", "vapour")),
        S("Condensate", ("V-100", "liquid"), ("V-101", "in")),
        S("Column feed", ("V-101", "out"), ("T-100", "in")),
        S("Ethane overhead", ("T-100", "distillate")),
        S("C3+", ("T-100", "bottoms"), ("T-101", "in")),
        S("LPG", ("T-101", "distillate")),
        S("Natural gasoline", ("T-101", "bottoms")),
    ],
}

methanol = {
    "name": "Methanol synthesis loop",
    "description": "Syngas reacts over a conversion reactor; crude methanol is condensed and unreacted gas recycled with a purge.",
    "components": ["carbon monoxide", "carbon dioxide", "hydrogen", "methane", "nitrogen", "methanol", "water"],
    "package": "Peng-Robinson",
    "units": [
        U("MIX", "Mixer", 160, 200),
        U("PRE", "Heater", 320, 200, mode="temperature", T_out=230.0, dP=0.0),
        U("R-1", "ConversionReactor", 480, 200, reactions="CO + 2 H2 -> CH3OH : CO : 0.35\nCO2 + 3 H2 -> CH3OH + H2O : CO2 : 0.20",
          mode="isothermal", dP=1.0),
        U("COOL", "Heater", 640, 200, mode="temperature", T_out=40.0, dP=0.0),
        U("SEP", "Flash", 800, 200, mode="T and P", T=40.0, P=79.0),
        U("PURGE", "Splitter", 800, 60, frac1=0.04),
        U("RECY", "Valve", 480, 60, P_out=80.0),
    ],
    "streams": [
        S("Syngas", to=("MIX", "in"), feed={"T": 40 + C, "P": 80e5, "flows": {
            "carbon monoxide": 30.0, "carbon dioxide": 5.0, "hydrogen": 74.0, "methane": 1.0, "nitrogen": 0.5}}),
        S("Loop feed", ("MIX", "out"), ("PRE", "in")),
        S("Reactor in", ("PRE", "out"), ("R-1", "in")),
        S("Reactor out", ("R-1", "out"), ("COOL", "in")),
        S("Cooled", ("COOL", "out"), ("SEP", "in")),
        S("Crude methanol", ("SEP", "liquid")),
        S("Loop gas", ("SEP", "vapour"), ("PURGE", "in")),
        S("Purge", ("PURGE", "out1")),
        S("Recycle", ("PURGE", "out2"), ("RECY", "in")),
        S("Recycle 2", ("RECY", "out"), ("MIX", "in")),
    ],
}

ore_tph = 100.0
ore = ore_tph / 3.6
copper = {
    "name": "Copper concentrator with concentrate leach",
    "description": "Closed-circuit ball mill with hydrocyclones, rougher-cleaner flotation, thickening and tank leaching of the concentrate.",
    "components": ["water"],
    "minerals": ["chalcopyrite", "pyrite", "quartz"],
    "package": "Ideal",
    "units": [
        U("SUMP", "Mixer", 160, 220),
        U("BM", "Comminution", 320, 220, P80=250.0, spread=0.9),
        U("CYC", "Classifier", 480, 220, d50c=110.0, sharpness=3.0, Rf=0.30),
        U("ROUGH", "FlotationBank", 640, 320, tau=14.0, cells=7, water_rec=0.12),
        U("CLEAN", "FlotationBank", 800, 220, tau=10.0, cells=5, water_rec=0.35, entrainment=0.3),
        U("TAILMIX", "Mixer", 960, 420),
        U("THK1", "Thickener", 960, 160, solids_rec=0.995, uf_solids=65.0),
        U("LEACH", "LeachTank", 1120, 160, extraction="chalcopyrite: 0.92\npyrite: 0.08", T=85.0, tau=360.0),
        U("CCD", "Thickener", 1280, 160, solids_rec=0.998, uf_solids=55.0),
    ],
    "streams": [
        S("Ore", to=("SUMP", "in"), feed={"T": 20 + C, "P": 101325.0, "flows": {"water": 0.0},
                                          "solids": {"chalcopyrite": {"mass": 0.02 * ore, "p80": 9000},
                                                     "pyrite": {"mass": 0.05 * ore, "p80": 9000},
                                                     "quartz": {"mass": 0.93 * ore, "p80": 12000}}}),
        S("Mill water", to=("SUMP", "in"), feed={"T": 20 + C, "P": 101325.0, "flows": {"water": ore / 0.018015}}),
        S("Mill feed", ("SUMP", "out"), ("BM", "in")),
        S("Mill discharge", ("BM", "out"), ("CYC", "in")),
        S("Cyclone underflow", ("CYC", "coarse"), ("SUMP", "in")),
        S("Cyclone overflow", ("CYC", "fine"), ("ROUGH", "in")),
        S("Rougher conc", ("ROUGH", "concentrate"), ("CLEAN", "in")),
        S("Rougher tail", ("ROUGH", "tailings"), ("TAILMIX", "in")),
        S("Cleaner tail", ("CLEAN", "tailings"), ("TAILMIX", "in")),
        S("Final tailings", ("TAILMIX", "out")),
        S("Cleaner conc", ("CLEAN", "concentrate"), ("THK1", "in")),
        S("Process water", ("THK1", "overflow")),
        S("Thickened conc", ("THK1", "underflow"), ("LEACH", "in")),
        S("Leach discharge", ("LEACH", "out"), ("CCD", "in")),
        S("Pregnant liquor", ("CCD", "overflow")),
        S("Leach residue", ("CCD", "underflow")),
    ],
}

if __name__ == "__main__":
    for name, fs in (("ngl_fractionation", ngl), ("methanol_loop", methanol), ("copper_concentrator", copper)):
        with open(os.path.join(OUT, name + ".json"), "w") as f:
            json.dump(fs, f, indent=1)
        print("wrote", name)


# ---------------------------------------------------------------------------
# Neodymium routes (three feedstocks, one downstream design): see ree_routes.py
# ---------------------------------------------------------------------------
from procsim.ree import ROUTES, build  # noqa: E402

if __name__ == "__main__":
    for fname, cfg in ROUTES.items():
        with open(os.path.join(OUT, fname + ".json"), "w") as f:
            json.dump(build(cfg), f, indent=1)
        print("wrote", fname)
