"""Write Python-core results for each example, for the JS parity check."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from procsim import Flowsheet
from procsim.flowsheet import _jsonable
here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(here, "reference"), exist_ok=True)
for name in ["ngl_fractionation", "methanol_loop", "copper_concentrator", "neodymium_bastnaesite", "neodymium_monazite_mtweld", "neodymium_eneabba_concentrate"]:
    spec = json.load(open(os.path.join(here, "..", "procsim", "data", "examples", name + ".json")))
    res = Flowsheet(spec).solve().results()
    json.dump(_jsonable(res), open(os.path.join(here, "reference", name + ".json"), "w"))
    print("reference", name)
