"""Bundle web/app.html + web/engine.js + the example flowsheets into one self-contained page.

Writes procsim/data/procsim.html, which ships inside the package and is opened by `procsim editor`.
Run from the repository root after editing web/app.html, web/engine.js or the examples: python web/build.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXAMPLES = ("ngl_fractionation", "methanol_loop", "copper_concentrator", "neodymium_bastnaesite",
            "neodymium_monazite_mtweld", "neodymium_eneabba_concentrate")


def bundle():
    app = open(os.path.join(HERE, "app.html"), encoding="utf-8").read()
    engine = open(os.path.join(HERE, "engine.js"), encoding="utf-8").read()
    exs = [json.load(open(os.path.join(ROOT, "procsim", "data", "examples", n + ".json"), encoding="utf-8")) for n in EXAMPLES]
    return app.replace("/*ENGINE*/", engine).replace("/*EXAMPLES*/", json.dumps(exs))


if __name__ == "__main__":
    out = os.path.join(ROOT, "procsim", "data", "procsim.html")
    html = bundle()
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    print(out, len(html) // 1024, "KB")
