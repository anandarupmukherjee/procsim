"""Neodymium routes from three feedstocks, built on one downstream design.

  A. Bastnaesite carbonatite ore (Mountain Pass type)
  B. Monazite carbonatite laterite ore (Mt Weld type)
  C. Monazite-xenotime mineral-sands concentrate (Eneabba type), fed straight to the refinery

Every route runs: [comminution and flotation] -> acid digestion and leach -> Th/Fe removal ->
MREC -> HCl redissolution -> SX-1 (Nd/Sm) -> SX-2 (Ce/Pr) -> SX-3 (Pr/Nd) -> Nd oxide -> Nd metal,
with Pr oxide, SEG oxide and LaCe carbonate as co-products.

Sources for the feed figures are given in each config; anything without a source is marked ASSUMED.

Typical use::

    from procsim import Flowsheet
    from procsim.ree import ROUTES, build, report
    cfg = ROUTES["neodymium_monazite_mtweld"]
    fs = Flowsheet(build(cfg)).solve()
    print(report(fs, cfg)["nd_oxide_tph"])
"""
import copy

C = 273.15
W = 18.015e-3  # kg/mol water


def S(id, frm=None, to=None, feed=None):
    d = {"id": id, "from": None, "to": None}
    if frm:
        d["from"] = {"unit": frm[0], "port": frm[1]}
    if to:
        d["to"] = {"unit": to[0], "port": to[1]}
    if feed:
        d["feed"] = feed
    return d


def U(id, type, x, y, **params):
    return {"id": id, "type": type, "x": x, "y": y, "params": params}


def feed(flows, T=20.0):
    return {"T": T + C, "P": 101325.0, "flows": flows}


def kin(table):
    return "\n".join(f"{m}: k={k}, Rmax={r}" for m, (k, r) in table.items())


def lines(table, fmt="{k}: {v}"):
    return "\n".join(fmt.format(k=k, v=v) for k, v in table.items())


RE_IONS = ["La3+", "Ce3+", "Pr3+", "Nd3+", "Sm+HREE3+"]
CARB = {"La3+": "La2(CO3)3", "Ce3+": "Ce2(CO3)3", "Pr3+": "Pr2(CO3)3", "Nd3+": "Nd2(CO3)3", "Sm+HREE3+": "Sm2(CO3)3"}


def build(cfg):
    """Return a flowsheet spec for one route configuration."""
    value = cfg["value_mineral"]
    gangue = cfg["gangue"]
    ore_route = "ore" in cfg
    minerals = [value] + list(gangue) + ["ThO2", "La2(CO3)3", "Ce2(CO3)3", "Pr2(CO3)3", "Nd2(CO3)3",
                                         "Sm2(CO3)3", "CeO2", "Pr6O11", "Nd2O3", "Sm2O3", "Nd metal"]
    for m in cfg.get("extra_minerals", []):
        if m not in minerals:
            minerals.append(m)
    units, streams = [], []

    if ore_route:
        ore = cfg["ore"]
        kg_s = ore["tph"] / 3.6
        fl = cfg["flotation"]
        units += [
            U("CR-100", "Comminution", 150, 140, P80=12000.0, spread=0.9),
            U("SUMP", "Mixer", 330, 160),
            U("BM-100", "Comminution", 450, 140, P80=fl["grind_p80"], spread=0.9),
            U("CY-100", "Classifier", 650, 110, d50c=fl["d50c"], sharpness=3.0, Rf=0.35),
            U("FL-ROUGH", "FlotationBank", 840, 320, tau=fl["tau"][0], cells=8, water_rec=0.12, entrainment=0.5, kinetics=kin(fl["kinetics"][0])),
            U("FL-CL1", "FlotationBank", 1060, 170, tau=fl["tau"][1], cells=6, water_rec=0.30, entrainment=0.3, kinetics=kin(fl["kinetics"][1])),
            U("FL-CL2", "FlotationBank", 1290, 170, tau=fl["tau"][2], cells=5, water_rec=0.35, entrainment=0.2, kinetics=kin(fl["kinetics"][2])),
            U("FL-CL3", "FlotationBank", 1520, 170, tau=fl["tau"][3], cells=5, water_rec=0.40, entrainment=0.1, kinetics=kin(fl["kinetics"][3])),
            U("TAILMIX", "Mixer", 1080, 480),
            U("TK-100", "Thickener", 1760, 170, solids_rec=0.995, uf_solids=65.0),
        ]
        streams += [
            S("Run-of-mine ore", to=("CR-100", "in"), feed={"T": 20 + C, "P": 101325.0, "flows": {"water": 0.0},
              "solids": {m: {"mass": f * kg_s, "p80": 150000, "spread": 0.8} for m, f in ore["composition"].items()}}),
            S("Crushed ore", ("CR-100", "out"), ("SUMP", "in")),
            S("Mill water", to=("SUMP", "in"), feed=feed({"water": 1.6 * kg_s / W})),
            S("Mill feed", ("SUMP", "out"), ("BM-100", "in")),
            S("Mill discharge", ("BM-100", "out"), ("CY-100", "in")),
            S("Cyclone underflow", ("CY-100", "coarse"), ("SUMP", "in")),
            S("Cyclone overflow", ("CY-100", "fine"), ("FL-ROUGH", "in")),
            S("Rougher conc", ("FL-ROUGH", "concentrate"), ("FL-CL1", "in")),
            S("Rougher tail", ("FL-ROUGH", "tailings"), ("TAILMIX", "in")),
            S("Cleaner 1 conc", ("FL-CL1", "concentrate"), ("FL-CL2", "in")),
            S("Cleaner 1 tail", ("FL-CL1", "tailings"), ("FL-ROUGH", "in")),
            S("Cleaner 2 conc", ("FL-CL2", "concentrate"), ("FL-CL3", "in")),
            S("Cleaner 2 tail", ("FL-CL2", "tailings"), ("FL-CL1", "in")),
            S("Cleaner 3 tail", ("FL-CL3", "tailings"), ("FL-CL2", "in")),
            S("Cleaner 3 conc", ("FL-CL3", "concentrate"), ("TK-100", "in")),
            S("Flotation tailings", ("TAILMIX", "out")),
            S("Process water return", ("TK-100", "overflow")),
            S(f"{cfg['concentrate_name']}", ("TK-100", "underflow"), ("LX-100", "in")),
        ]
    else:
        conc = cfg["concentrate"]
        kg_s = conc["tph"] / 3.6
        streams += [
            S(cfg["concentrate_name"], to=("LX-100", "in"), feed={"T": 20 + C, "P": 101325.0, "flows": {"water": 0.35 * kg_s / W},
              "solids": {m: {"mass": f * kg_s, "p80": conc["p80"], "spread": 1.2} for m, f in conc["composition"].items()}}),
        ]

    lx = cfg["leach"]
    re_ext = {m: v for m, v in lx["extraction"].items()}
    rg = cfg["reagents"]
    sx = cfg["sx"]
    prec_mrec = {ion: (CARB[ion], 0.995) for ion in RE_IONS}
    units += [
        U("LX-100", "LeachTank", 200, 680, T=lx["T"], tau=lx["tau"], extraction=lines(re_ext)),
        U("FT-100", "Thickener", 440, 690, solids_rec=0.998, uf_solids=80.0),
        U("PR-100", "Precipitator", 740, 670, T=60.0, moisture=40.0, p80=15.0,
          reactions=lines(cfg["impurity_removal"], "{k} : {v}")),
        U("PR-200", "Precipitator", 1040, 670, T=60.0, moisture=30.0, p80=25.0,
          reactions="\n".join(f"{ion} -> {s} : {x}" for ion, (s, x) in prec_mrec.items()) + "\nCa2+ -> calcite : 0.15"),
        U("LX-200", "LeachTank", 1340, 670, T=60.0, tau=60.0,
          extraction="La2(CO3)3: 1\nCe2(CO3)3: 1\nPr2(CO3)3: 1\nNd2(CO3)3: 1\nSm2(CO3)3: 1\ncalcite: 1"),
        U("SX-1", "SXCascade", 220, 960, key="Nd3+", strip_water=0.1, **sx["SX-1"]),
        U("SX-2", "SXCascade", 590, 960, key="Pr3+", strip_water=0.3, **sx["SX-2"]),
        U("SX-3", "SXCascade", 960, 960, key="Nd3+", strip_water=0.8, **sx["SX-3"]),
        U("PR-300", "Precipitator", 1300, 950, T=60.0, moisture=20.0, p80=10.0,
          reactions="Nd3+ -> Nd2O3 : 0.995\nPr3+ -> Pr6O11 : 0.995\nSm+HREE3+ -> Sm2O3 : 0.995"),
        U("CAL-100", "Flash", 1540, 930, mode="T and P", T=900.0, P=1.01325),
        U("EC-100", "Electrolysis", 1760, 950, T=1050.0, energy=10.0, current_eff=0.75, reactions="Nd2O3 -> Nd metal : 0.98"),
        U("PR-400", "Precipitator", 1180, 1240, T=60.0, moisture=20.0, p80=10.0,
          reactions="Pr3+ -> Pr6O11 : 0.995\nNd3+ -> Nd2O3 : 0.995\nCe3+ -> CeO2 : 0.995"),
        U("PR-500", "Precipitator", 420, 1240, T=60.0, moisture=20.0, p80=10.0,
          reactions="Sm+HREE3+ -> Sm2O3 : 0.995\nNd3+ -> Nd2O3 : 0.995\nPr3+ -> Pr6O11 : 0.995"),
        U("PR-600", "Precipitator", 800, 1240, T=60.0, moisture=30.0, p80=25.0,
          reactions="La3+ -> La2(CO3)3 : 0.995\nCe3+ -> Ce2(CO3)3 : 0.995\nPr3+ -> Pr2(CO3)3 : 0.995"),
    ]
    streams += [
        S("Acid and dilution water", to=("LX-100", "in"), feed=feed({"water": rg["leach_water"]})),
        S("Leach slurry", ("LX-100", "out"), ("FT-100", "in")),
        S("Leach residue", ("FT-100", "underflow")),
        S("Pregnant leach solution", ("FT-100", "overflow"), ("PR-100", "in")),
        S("Neutraliser", to=("PR-100", "in"), feed=feed({"water": 20.0, "reagent": rg["neutraliser"]})),
        S("Th residue", ("PR-100", "solids")),
        S("Purified liquor", ("PR-100", "liquor"), ("PR-200", "in")),
        S("Carbonate precipitant", to=("PR-200", "in"), feed=feed({"water": 600.0, "reagent": rg["mrec"]})),
        S("MREC", ("PR-200", "solids"), ("LX-200", "in")),
        S("Barren liquor", ("PR-200", "liquor")),
        S("HCl redissolution", to=("LX-200", "in"), feed=feed({"water": rg["hcl_water"]})),
        S("RE chloride feed", ("LX-200", "out"), ("SX-1", "in")),
        S("SEG strip liquor", ("SX-1", "extract"), ("PR-500", "in")),
        S("La-Nd raffinate", ("SX-1", "raffinate"), ("SX-2", "in")),
        S("PrNd strip liquor", ("SX-2", "extract"), ("SX-3", "in")),
        S("LaCe raffinate", ("SX-2", "raffinate"), ("PR-600", "in")),
        S("Nd strip liquor", ("SX-3", "extract"), ("PR-300", "in")),
        S("Pr raffinate", ("SX-3", "raffinate"), ("PR-400", "in")),
        S("Oxalic acid (Nd)", to=("PR-300", "in"), feed=feed({"water": 30.0, "reagent": rg["ox_nd"]})),
        S("Nd oxalate cake", ("PR-300", "solids"), ("CAL-100", "in")),
        S("Nd barren liquor", ("PR-300", "liquor")),
        S("Calciner offgas", ("CAL-100", "vapour")),
        S("Nd oxide", ("CAL-100", "liquid"), ("EC-100", "in")),
        S("Nd metal", ("EC-100", "metal")),
        S("Anode gas", ("EC-100", "offgas")),
        S("Oxalic acid (Pr)", to=("PR-400", "in"), feed=feed({"water": 20.0, "reagent": rg["ox_pr"]})),
        S("Pr oxide cake", ("PR-400", "solids")),
        S("Pr barren liquor", ("PR-400", "liquor")),
        S("Oxalic acid (SEG)", to=("PR-500", "in"), feed=feed({"water": 10.0, "reagent": rg["ox_seg"]})),
        S("SEG oxide cake", ("PR-500", "solids")),
        S("SEG barren liquor", ("PR-500", "liquor")),
        S("Carbonate precipitant (LaCe)", to=("PR-600", "in"), feed=feed({"water": 300.0, "reagent": rg["lace"]})),
        S("LaCe carbonate", ("PR-600", "solids")),
        S("LaCe barren liquor", ("PR-600", "liquor")),
    ]
    if not ore_route:  # no comminution or flotation band: lift the sheet
        for u in units:
            u["y"] -= 520
    return {"name": cfg["name"], "description": cfg["description"], "sources": cfg["sources"],
            "components": ["water", "oxygen", "reagent"], "minerals": minerals, "package": "Ideal",
            "units": units, "streams": streams}


# ---------------------------------------------------------------------------
# Route configurations
# ---------------------------------------------------------------------------
BASTNAESITE = {
    "name": "Neodymium: bastnaesite ore to Nd metal",
    "description": ("Mountain Pass-type carbonatite ore (8.0 % REO) through crushing, grinding, rougher and three cleaner "
                    "flotation stages, acid digestion, thorium removal, mixed RE carbonate (MREC), HCl redissolution, three "
                    "solvent-extraction splits (Nd/Sm, Ce/Pr, Pr/Nd), oxalate precipitation with calcination, and "
                    "molten-salt electrolysis to Nd metal."),
    "sources": {
        "RE distribution": "USGS Minerals Yearbook 2020, Rare earths, Table 2, Mountain Pass row (p60.8)",
        "Feed and concentrate grade": "Rozelle, Tarka and Mamula (OSTI 2020), Table 2: feed 7.7-8.4 % REO, concentrate 58.0-63 % REO (p19)",
        "Gangue split, Th, kinetics, SX factors": "ASSUMED",
    },
    "value_mineral": "bastnaesite",
    "gangue": ["calcite", "barite", "quartz"],
    "concentrate_name": "Bastnaesite concentrate",
    "ore": {"tph": 100.0, "composition": {"bastnaesite": 0.105, "calcite": 0.40, "barite": 0.25, "quartz": 0.245}},
    "flotation": {"grind_p80": 110.0, "d50c": 70.0, "tau": [16.0, 12.0, 10.0, 10.0], "kinetics": [
        {"bastnaesite": (0.75, 0.93), "calcite": (0.12, 0.50), "barite": (0.14, 0.50), "quartz": (0.01, 0.10)},
        {"bastnaesite": (0.80, 0.95), "calcite": (0.11, 0.50), "barite": (0.13, 0.50), "quartz": (0.01, 0.10)},
        {"bastnaesite": (0.80, 0.96), "calcite": (0.10, 0.50), "barite": (0.12, 0.50), "quartz": (0.01, 0.10)},
        {"bastnaesite": (0.80, 0.97), "calcite": (0.09, 0.50), "barite": (0.11, 0.50), "quartz": (0.01, 0.10)}]},
    "leach": {"T": 95.0, "tau": 480.0, "extraction": {"bastnaesite": 0.95, "calcite": 1.0}},
    "impurity_removal": {"Th4+ -> ThO2": 0.995},
    "sx": {"SX-1": dict(D_key=0.5, N_ext=12, N_scrub=12, OA=0.5, SA=0.3),
           "SX-2": dict(D_key=2.6, N_ext=30, N_scrub=30, OA=1.3, SA=1.9),
           "SX-3": dict(D_key=2.6, N_ext=70, N_scrub=70, OA=1.6, SA=3.0)},
    "reagents": {"leach_water": 250.0, "neutraliser": 0.2, "mrec": 13.0, "hcl_water": 250.0,
                 "ox_nd": 0.5, "ox_pr": 0.2, "ox_seg": 0.1, "lace": 11.0},
}

MTWELD = {
    "name": "Neodymium: Mt Weld-type monazite ore to Nd metal",
    "description": ("Weathered carbonatite ore at the Mt Weld reserve grade (6.4 % TREO) with monazite as the RE mineral; "
                    "same comminution, flotation and downstream design as the bastnaesite route, with sulphuric-acid "
                    "cracking and leaching, Th and Fe removal and MREC as the intermediate (the Kalgoorlie product)."),
    "sources": {
        "Ore grade": "Lynas Appendix 4E / Annual Report 2025, Table 1: Ore Reserves 32.0 Mt at 6.4 % TREO (p116); proved 7.4 %",
        "Monazite as RE mineral": "Lynas Annual Report 2024, Directors' Report (p14): monazite is the major RE mineral in the dolomite carbonatite",
        "RE distribution": "USGS Minerals Yearbook 2020, Rare earths, Table 2, Mt Weld Central Lanthanide row (p60.8)",
        "Gangue split, Th, concentrate grade, kinetics, SX factors": "ASSUMED (no figure in the Data Library)",
    },
    "value_mineral": "monazite",
    "gangue": ["hematite", "apatite", "quartz"],
    "extra_minerals": ["calcite"],
    "concentrate_name": "Monazite concentrate",
    "ore": {"tph": 100.0, "composition": {"monazite": 0.097, "hematite": 0.35, "apatite": 0.15, "quartz": 0.403}},
    "flotation": {"grind_p80": 75.0, "d50c": 50.0, "tau": [18.0, 12.0, 10.0, 10.0], "kinetics": [
        {"monazite": (0.60, 0.90), "hematite": (0.06, 0.40), "apatite": (0.20, 0.60), "quartz": (0.01, 0.10)},
        {"monazite": (0.65, 0.92), "hematite": (0.05, 0.40), "apatite": (0.18, 0.60), "quartz": (0.01, 0.10)},
        {"monazite": (0.65, 0.93), "hematite": (0.05, 0.40), "apatite": (0.17, 0.60), "quartz": (0.01, 0.10)},
        {"monazite": (0.65, 0.94), "hematite": (0.05, 0.40), "apatite": (0.16, 0.60), "quartz": (0.01, 0.10)}]},
    "leach": {"T": 95.0, "tau": 480.0, "extraction": {"monazite": 0.95, "apatite": 1.0, "hematite": 0.30}},
    "impurity_removal": {"Th4+ -> ThO2": 0.995, "Fe2+ -> hematite": 0.99},
    "sx": {"SX-1": dict(D_key=0.5, N_ext=12, N_scrub=12, OA=0.5, SA=0.3),
           "SX-2": dict(D_key=2.6, N_ext=30, N_scrub=30, OA=1.3, SA=1.9),
           "SX-3": dict(D_key=2.6, N_ext=70, N_scrub=70, OA=1.6, SA=3.0)},
    "reagents": {"leach_water": 250.0, "neutraliser": 5.0, "mrec": 13.0, "hcl_water": 250.0,
                 "ox_nd": 0.5, "ox_pr": 0.2, "ox_seg": 0.2, "lace": 11.0},
}

ENEABBA = {
    "name": "Neodymium: Eneabba-type monazite-xenotime concentrate to Nd metal",
    "description": ("Mineral-sands rare earth concentrate at 57 % REO fed straight to the refinery at 32 kt per year "
                    "(Iluka's worked example), through the same cracking, purification, SX and finishing steps."),
    "sources": {
        "Feed rate, grade, refinery recovery check": "Iluka Eneabba refinery presentation 2024, p25: 32 kt concentrate at ~57 % REO, 90 % simplified weighted average recovery",
        "RE distribution": "Iluka Eneabba refinery presentation 2024, p25 (La, Ce, Pr, Nd, Tb, Dy) and p27 (Sm, Gd, Y)",
        "SX order": "Iluka Eneabba refinery presentation 2024, p25: SX1 light/heavy, SX2 LaCe/NdPr, SX3 and SX5 Pr/Nd",
        "Gangue, Th, SX factors": "ASSUMED",
    },
    "value_mineral": "monazite-xenotime",
    "gangue": ["hematite", "quartz"],
    "extra_minerals": ["calcite"],
    "concentrate_name": "Monazite-xenotime concentrate",
    "concentrate": {"tph": 32000 / 8000, "p80": 150.0,
                    "composition": {"monazite-xenotime": 0.866, "hematite": 0.05, "quartz": 0.084}},
    "leach": {"T": 95.0, "tau": 480.0, "extraction": {"monazite-xenotime": 0.95, "hematite": 0.30}},
    "impurity_removal": {"Th4+ -> ThO2": 0.995, "Fe2+ -> hematite": 0.99},
    "sx": {"SX-1": dict(D_key=0.5, N_ext=12, N_scrub=12, OA=0.5, SA=0.3),
           "SX-2": dict(D_key=2.6, N_ext=30, N_scrub=30, OA=1.3, SA=1.9),
           "SX-3": dict(D_key=2.6, N_ext=70, N_scrub=70, OA=1.6, SA=3.0)},
    "reagents": {"leach_water": 40.0, "neutraliser": 1.0, "mrec": 3.0, "hcl_water": 40.0,
                 "ox_nd": 0.2, "ox_pr": 0.1, "ox_seg": 0.1, "lace": 3.0},
}

ROUTES = {"neodymium_bastnaesite": BASTNAESITE, "neodymium_monazite_mtweld": MTWELD, "neodymium_eneabba_concentrate": ENEABBA}

# Values produced by tools/tune_routes.py (2026-10-09): Mt Weld cleaner gangue Rmax scaled to give a 40 % REO
# concentrate (ASSUMED target); reagents sized at 1.15 x consumption; SX ratios are the lowest O/A + S/A meeting
# 99.5 % (SX-1, SX-2) and 99 % (SX-3) key recoveries at the fixed stage counts.
TUNED = {
 "neodymium_bastnaesite": {
  "flotation": {
   "grind_p80": 110.0,
   "d50c": 70.0,
   "tau": [
    16.0,
    12.0,
    10.0,
    10.0
   ],
   "kinetics": [
    {
     "bastnaesite": [
      0.75,
      0.93
     ],
     "calcite": [
      0.12,
      0.5
     ],
     "barite": [
      0.14,
      0.5
     ],
     "quartz": [
      0.01,
      0.1
     ]
    },
    {
     "bastnaesite": [
      0.8,
      0.95
     ],
     "calcite": [
      0.11,
      0.5
     ],
     "barite": [
      0.13,
      0.5
     ],
     "quartz": [
      0.01,
      0.1
     ]
    },
    {
     "bastnaesite": [
      0.8,
      0.96
     ],
     "calcite": [
      0.1,
      0.5
     ],
     "barite": [
      0.12,
      0.5
     ],
     "quartz": [
      0.01,
      0.1
     ]
    },
    {
     "bastnaesite": [
      0.8,
      0.97
     ],
     "calcite": [
      0.09,
      0.5
     ],
     "barite": [
      0.11,
      0.5
     ],
     "quartz": [
      0.01,
      0.1
     ]
    }
   ]
  },
  "sx": {
   "SX-1": {
    "D_key": 0.8,
    "N_ext": 12,
    "N_scrub": 12,
    "OA": 0.3,
    "SA": 0.3
   },
   "SX-2": {
    "D_key": 2.8,
    "N_ext": 30,
    "N_scrub": 30,
    "OA": 1.2,
    "SA": 1.9
   },
   "SX-3": {
    "D_key": 2.9,
    "N_ext": 70,
    "N_scrub": 70,
    "OA": 1.4,
    "SA": 2.9
   }
  },
  "reagents": {
   "leach_water": 250.0,
   "neutraliser": 0.05,
   "mrec": 14.27,
   "hcl_water": 250.0,
   "ox_nd": 0.43,
   "ox_pr": 0.2,
   "ox_seg": 0.05,
   "lace": 11.41
  }
 },
 "neodymium_monazite_mtweld": {
  "flotation": {
   "grind_p80": 75.0,
   "d50c": 50.0,
   "tau": [
    18.0,
    12.0,
    10.0,
    10.0
   ],
   "kinetics": [
    {
     "monazite": [
      0.6,
      0.9
     ],
     "hematite": [
      0.06,
      0.4824
     ],
     "apatite": [
      0.2,
      0.7236
     ],
     "quartz": [
      0.01,
      0.1206
     ]
    },
    {
     "monazite": [
      0.65,
      0.92
     ],
     "hematite": [
      0.05,
      0.4824
     ],
     "apatite": [
      0.18,
      0.7236
     ],
     "quartz": [
      0.01,
      0.1206
     ]
    },
    {
     "monazite": [
      0.65,
      0.93
     ],
     "hematite": [
      0.05,
      0.4824
     ],
     "apatite": [
      0.17,
      0.7236
     ],
     "quartz": [
      0.01,
      0.1206
     ]
    },
    {
     "monazite": [
      0.65,
      0.94
     ],
     "hematite": [
      0.05,
      0.4824
     ],
     "apatite": [
      0.16,
      0.7236
     ],
     "quartz": [
      0.01,
      0.1206
     ]
    }
   ]
  },
  "sx": {
   "SX-1": {
    "D_key": 0.8,
    "N_ext": 12,
    "N_scrub": 12,
    "OA": 0.3,
    "SA": 0.3
   },
   "SX-2": {
    "D_key": 2.8,
    "N_ext": 30,
    "N_scrub": 30,
    "OA": 1.2,
    "SA": 1.9
   },
   "SX-3": {
    "D_key": 2.9,
    "N_ext": 70,
    "N_scrub": 70,
    "OA": 1.4,
    "SA": 2.9
   }
  },
  "reagents": {
   "leach_water": 250.0,
   "neutraliser": 0.26,
   "mrec": 10.77,
   "hcl_water": 250.0,
   "ox_nd": 0.44,
   "ox_pr": 0.16,
   "ox_seg": 0.13,
   "lace": 6.52
  }
 },
 "neodymium_eneabba_concentrate": {
  "flotation": None,
  "sx": {
   "SX-1": {
    "D_key": 0.8,
    "N_ext": 12,
    "N_scrub": 12,
    "OA": 0.3,
    "SA": 0.3
   },
   "SX-2": {
    "D_key": 2.8,
    "N_ext": 30,
    "N_scrub": 30,
    "OA": 1.2,
    "SA": 1.9
   },
   "SX-3": {
    "D_key": 2.9,
    "N_ext": 70,
    "N_scrub": 70,
    "OA": 1.4,
    "SA": 2.9
   }
  },
  "reagents": {"leach_water": 120.0, "neutraliser": 0.23, "mrec": 5.93, "hcl_water": 110.0, "ox_nd": 0.26, "ox_pr": 0.09, "ox_seg": 0.21, "lace": 3.86}
 }
}


def _apply(cfg, tuned):
    if tuned.get("flotation"):
        cfg["flotation"]["kinetics"] = [{m: tuple(v) for m, v in st.items()} for st in tuned["flotation"]["kinetics"]]
    for k, v in tuned["sx"].items():
        cfg["sx"][k].update(v)
    cfg["reagents"].update(tuned["reagents"])


for _name, _tuned in TUNED.items():
    _apply(ROUTES[_name], _tuned)



# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
_ND = 144.242
_OX_PER_RE = {"La3+": 162.905, "Ce3+": 172.115, "Pr3+": 170.24, "Nd3+": 168.24, "Sm+HREE3+": 174.36}
_ORE_MINERALS = ("bastnaesite", "monazite", "monazite-xenotime")
_LOSS_STREAMS = ["Flotation tailings", "Leach residue", "Th residue", "Barren liquor", "SEG oxide cake", "LaCe carbonate",
                 "Pr oxide cake", "Nd barren liquor", "Process water return", "Calciner offgas", "Anode gas"]
_PRODUCTS = ["Nd oxide", "Pr oxide cake", "SEG oxide cake", "LaCe carbonate", "Th residue", "Nd metal", "Leach residue"]


def report(fs, cfg):
    """Nd balance by stage, losses, product basket, reagents and energy for a solved route."""
    ctx = fs.ctx

    def nd(sid):  # t/h-equivalent basis: kg/s of contained Nd
        s = fs.stream(sid)
        t = s.F[ctx.index("Nd3+")] * _ND / 1000
        for i, m in enumerate(ctx.mins):
            n = m.leach_products.get("Nd3+")
            if m.name == "Nd metal":
                t += s.S[i].sum()
            elif n:
                t += s.S[i].sum() / (m.mw / 1000) * n * _ND / 1000
        return float(t)

    def reo(sid):
        s = fs.stream(sid)
        r = 0.0
        for i, m in enumerate(ctx.mins):
            if m.name in _ORE_MINERALS:
                r += s.S[i].sum() / (m.mw / 1000) * sum(v * _OX_PER_RE[k] for k, v in m.leach_products.items() if k in _OX_PER_RE) / 1000
        return float(r)

    first = "Run-of-mine ore" if "ore" in cfg else cfg["concentrate_name"]
    nd0 = nd(first)
    seq = ([first] + (["Cyclone overflow", cfg["concentrate_name"]] if "ore" in cfg else []) +
           ["Pregnant leach solution", "MREC", "RE chloride feed", "PrNd strip liquor", "Nd strip liquor", "Nd oxide"])
    out = {"route": cfg["name"]}
    out["nd_balance"] = [(s, nd(s) * 3.6, 100 * nd(s) / nd0) for s in seq]
    metal = float(fs.stream("Nd metal").S[ctx.mindex("Nd metal")].sum())
    out["nd_balance"].append(("Nd metal (as metal)", metal * 3.6, 100 * metal / nd0))
    c = fs.stream(cfg["concentrate_name"])
    out["feed_REO"] = reo(first) / fs.stream(first).solids_mass
    out["conc_tph"] = c.solids_mass * 3.6
    out["conc_REO"] = reo(cfg["concentrate_name"]) / c.solids_mass
    out["concentrate_to_oxide"] = nd("Nd oxide") / nd(cfg["concentrate_name"])
    out["losses_pct"] = {s: 100 * nd(s) / nd0 for s in _LOSS_STREAMS if s in fs.values and nd(s) / nd0 > 1e-4}
    out["products_tph"] = {sid: {m: float(fs.stream(sid).S[i].sum() * 3.6) for i, m in enumerate(ctx.minerals)
                                 if fs.stream(sid).S[i].sum() > 1e-7} for sid in _PRODUCTS}
    ox = fs.stream("Nd oxide")
    out["nd_oxide_tph"] = ox.solids_mass * 3.6
    out["nd2o3_purity"] = float(ox.S[ctx.mindex("Nd2O3")].sum() / ox.solids_mass)
    out["metal_tph"] = metal * 3.6
    out["reagent_tph"] = {u: float(fs.units[u].results["Reagent consumed [t/h]"]) for u in fs.units
                          if "Reagent consumed [t/h]" in fs.units[u].results}
    out["power_kW"] = {u: float(fs.units[u].results.get("Power [kW]", fs.units[u].results.get("Electrical power [kW]", 0)))
                       for u in ["CR-100", "BM-100", "EC-100"] if u in fs.units}
    w = fs.stream("RE chloride feed")
    out["sx_feed_mol_per_L"] = sum(w.F[ctx.index(n)] for n in _OX_PER_RE) / (w.F[ctx.index("water")] * 0.018015)
    out["mass_closure_pct"] = fs.balances()["mass_closure_pct"]
    return out


def compare(names=None):
    """Solve routes and return one summary row per route (all ratios per tonne of Nd oxide)."""
    from .flowsheet import Flowsheet
    rows = {}
    for name in names or ROUTES:
        cfg = ROUTES[name]
        fs = Flowsheet(build(cfg)).solve()
        r = report(fs, cfg)
        p, ox = r["products_tph"], r["nd_oxide_tph"]
        rows[name] = {
            "feed_REO_pct": 100 * r["feed_REO"],
            "concentrate_REO_pct": 100 * r["conc_REO"],
            "nd_to_metal_pct": r["nd_balance"][-1][2],
            "concentrate_to_oxide_pct": 100 * r["concentrate_to_oxide"],
            "nd_oxide_tph": ox,
            "nd_oxide_purity_pct": 100 * r["nd2o3_purity"],
            "LaCe_per_t": sum(p["LaCe carbonate"].values()) / ox,
            "Pr_oxide_per_t": sum(p["Pr oxide cake"].values()) / ox,
            "SEG_per_t": sum(p["SEG oxide cake"].values()) / ox,
            "ThO2_per_t": p["Th residue"].get("ThO2", 0.0) / ox,
            "reagent_per_t": sum(r["reagent_tph"].values()) / ox,
        }
    return rows
