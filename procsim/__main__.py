"""procsim command line.

    procsim run FLOWSHEET.json [-o results.json] [--json]   solve a flowsheet file or a bundled example
    procsim examples [--copy DIR]                             list (or copy out) the bundled examples
    procsim units                                             list unit operations
    procsim editor [--out FILE] [--no-open]                   open the drag-and-drop editor in a browser
    procsim routes [NAME ...] [--json]                        compare the neodymium feedstock routes
"""
import argparse
import json
import os
import shutil
import sys
import tempfile
import webbrowser
from importlib import resources

from . import __version__
from .flowsheet import Flowsheet, _jsonable


def _data(*parts):
    return resources.files("procsim").joinpath("data", *parts)


def example_names():
    return sorted(p.name[:-5] for p in _data("examples").iterdir() if p.name.endswith(".json"))


def load_spec(ref):
    """A path to a JSON file, or the name of a bundled example."""
    if os.path.exists(ref):
        with open(ref, encoding="utf-8") as f:
            return json.load(f)
    name = ref[:-5] if ref.endswith(".json") else ref
    if name in example_names():
        return json.loads(_data("examples", name + ".json").read_text(encoding="utf-8"))
    raise SystemExit(f"procsim: '{ref}' is not a file or a bundled example. Bundled: {', '.join(example_names())}")


def _fmt(v):
    return f"{v:.4g}" if isinstance(v, float) else str(v)


def cmd_run(a):
    fs = Flowsheet(load_spec(a.flowsheet)).solve()
    res = _jsonable(fs.results())
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=1)
    if a.json:
        json.dump(res, sys.stdout, indent=1)
        print()
        return 0
    print(fs.table())
    print()
    for uid, u in fs.units.items():
        print(f"{uid} ({u.TYPE})")
        for k, v in u.results.items():
            print(f"   {k}: {_fmt(float(v) if hasattr(v, 'dtype') else v)}")
    for lp in fs.loops:
        state = "converged" if lp["converged"] else "NOT converged"
        print(f"Recycle {lp['units']} tears {lp['tears']}: {state} in {lp['iterations']} iterations")
    b = fs.balances()
    print(f"Mass closure {b['mass_closure_pct']:.2e} %, energy closure {b['energy_closure_kW']:.3g} kW")
    for w in fs.log:
        print("WARNING:", w)
    return 0 if all(lp["converged"] for lp in fs.loops) else 2


def cmd_examples(a):
    if a.copy:
        os.makedirs(a.copy, exist_ok=True)
        for n in example_names():
            dest = os.path.join(a.copy, n + ".json")
            with open(dest, "w", encoding="utf-8") as f:
                f.write(_data("examples", n + ".json").read_text(encoding="utf-8"))
            print(dest)
        return 0
    for n in example_names():
        spec = json.loads(_data("examples", n + ".json").read_text(encoding="utf-8"))
        print(f"{n:32s} {spec.get('name', '')}")
    return 0


def cmd_units(a):
    from .units import UNIT_TYPES
    for k, c in UNIT_TYPES.items():
        print(f"{k:20s} [{c.CATEGORY}] {c.DESCRIPTION}")
    return 0


def cmd_editor(a):
    out = a.out or os.path.join(tempfile.gettempdir(), "procsim-editor.html")
    with resources.as_file(_data("procsim.html")) as src:
        shutil.copyfile(src, out)
    print(out)
    if not a.no_open:
        webbrowser.open("file://" + os.path.abspath(out))
    return 0


def cmd_routes(a):
    from .ree import ROUTES, compare
    names = a.names or list(ROUTES)
    unknown = [n for n in names if n not in ROUTES]
    if unknown:
        raise SystemExit(f"procsim: unknown route(s) {unknown}. Routes: {', '.join(ROUTES)}")
    rows = compare(names)
    if a.json:
        json.dump(_jsonable(rows), sys.stdout, indent=1)
        print()
        return 0
    labels = [("feed_REO_pct", "Feed REO, %"), ("concentrate_REO_pct", "Concentrate REO, %"),
              ("nd_to_metal_pct", "Nd feed to metal, %"), ("concentrate_to_oxide_pct", "Concentrate to Nd oxide, %"),
              ("nd_oxide_tph", "Nd oxide, t/h"), ("nd_oxide_purity_pct", "Nd oxide purity (TREO), %"),
              ("LaCe_per_t", "LaCe carbonate per t Nd oxide"), ("Pr_oxide_per_t", "Pr oxide per t Nd oxide"),
              ("SEG_per_t", "SEG oxide per t Nd oxide"), ("ThO2_per_t", "ThO2 residue per t Nd oxide"),
              ("reagent_per_t", "Precipitant per t Nd oxide")]
    short = {n: n.replace("neodymium_", "") for n in names}
    w0 = max(len(l) for _, l in labels)
    print(" " * w0 + "".join(f"{short[n]:>24s}" for n in names))
    for key, label in labels:
        print(f"{label:{w0}s}" + "".join(f"{rows[n][key]:>24.4g}" for n in names))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="procsim", description="Steady-state process simulator for chemical, mineral and rare-earth flowsheets.")
    ap.add_argument("--version", action="version", version=f"procsim {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="Solve a flowsheet JSON file or a bundled example")
    r.add_argument("flowsheet", help="Path to a flowsheet JSON file, or a bundled example name")
    r.add_argument("-o", "--out", help="Write full results JSON here")
    r.add_argument("--json", action="store_true", help="Print results as JSON instead of tables")
    r.set_defaults(fn=cmd_run)
    e = sub.add_parser("examples", help="List the bundled example flowsheets")
    e.add_argument("--copy", metavar="DIR", help="Copy the example JSON files into DIR")
    e.set_defaults(fn=cmd_examples)
    sub.add_parser("units", help="List unit operations").set_defaults(fn=cmd_units)
    ed = sub.add_parser("editor", help="Open the flowsheet editor (a single offline HTML page) in a browser")
    ed.add_argument("--out", help="Where to write the HTML page (default: system temp folder)")
    ed.add_argument("--no-open", action="store_true", help="Write the page without opening a browser")
    ed.set_defaults(fn=cmd_editor)
    ro = sub.add_parser("routes", help="Compare the neodymium routes by feedstock")
    ro.add_argument("names", nargs="*", help="Route names (default: all)")
    ro.add_argument("--json", action="store_true", help="Print the comparison as JSON")
    ro.set_defaults(fn=cmd_routes)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
