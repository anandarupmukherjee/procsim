"""Print the full Nd report for one or more routes: python tools/route_report.py neodymium_monazite_mtweld"""
import sys

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__))))  # checkout without install
from procsim import Flowsheet
from procsim.ree import ROUTES, build, report

if __name__ == "__main__":
    for name in sys.argv[1:] or list(ROUTES):
        cfg = ROUTES[name]
        r = report(Flowsheet(build(cfg)).solve(), cfg)
        print("=====", name)
        for k, v in r.items():
            print(k, v)
