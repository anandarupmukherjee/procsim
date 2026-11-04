"""Recalibrate the Nd routes: Mt Weld flotation to a target grade, reagent sizing, SX ratios.
Paste the printed values into TUNED in procsim/ree.py. Usage: python tools/tune_routes.py neodymium_bastnaesite ...
"""
import json, sys, copy, itertools
import numpy as np
import os
_H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__))))  # checkout without install
from procsim import Flowsheet
from procsim.units import SXCascade
from procsim.flowsheet import FlowsheetError
from procsim.ree import ROUTES, build, report
import re as _re

def solve(cfg):
    for _ in range(12):
        try:
            return Flowsheet(build(cfg)).solve()
        except FlowsheetError as e:
            m = _re.search(r'(PR-\d+) \(Precipitator\): Not enough reagent: needs ([\d.e+-]+) t/h', str(e))
            if not m: raise
            key = {'PR-100':'neutraliser','PR-200':'mrec','PR-300':'ox_nd','PR-400':'ox_pr','PR-500':'ox_seg','PR-600':'lace'}[m.group(1)]
            cfg['reagents'][key] = round(float(m.group(2)) / 3.6 / 0.060 * 1.25 + 0.05, 2)
    raise RuntimeError('reagent sizing did not settle')

def size_reagents(cfg):
    fs = solve(cfg)
    cons = {u: fs.units[u].results['Reagent consumed [t/h]'] for u in fs.units if 'Reagent consumed [t/h]' in fs.units[u].results}
    m = {'PR-100':'neutraliser','PR-200':'mrec','PR-300':'ox_nd','PR-400':'ox_pr','PR-500':'ox_seg','PR-600':'lace'}
    for u, k in m.items():
        cfg['reagents'][k] = round(max(cons[u] / 3.6 / 0.060 * 1.15, 0.05), 2)
    return solve(cfg)

def tune_sx(cfg, fs):
    ctx = fs.ctx
    cases = {'SX-1': ('RE chloride feed', 'Nd3+', 'Nd3+', 'Sm+HREE3+', 0.995),
             'SX-2': ('La-Nd raffinate', 'Pr3+', 'Ce3+', 'Pr3+', 0.995),
             'SX-3': ('PrNd strip liquor', 'Nd3+', 'Pr3+', 'Nd3+', 0.99)}
    for name, (sid, key, light, heavy, tgt) in cases.items():
        feed = fs.stream(sid); Ne, Ns = cfg['sx'][name]['N_ext'], cfg['sx'][name]['N_scrub']
        hi, li = ctx.index(heavy), ctx.index(light)
        best = None
        for D in np.round(np.arange(0.2, 3.01, 0.1), 2):
            for OA in np.round(np.arange(0.3, 4.01, 0.1), 2):
                for SA in np.round(np.arange(0.1, 4.01, 0.1), 2):
                    if best and OA + SA >= best[0]: continue
                    u = SXCascade('t', dict(key=key, D_key=D, N_ext=Ne, N_scrub=Ns, OA=OA, SA=SA))
                    o = u.run({'in': [feed]}, ctx)
                    rh = o['extract'].F[hi] / feed.F[hi]; rl = o['raffinate'].F[li] / feed.F[li]
                    if min(rh, rl) >= tgt:
                        best = (OA + SA, float(D), float(OA), float(SA), rh, rl)
        if best is None: raise RuntimeError(f'{name}: target not met')
        cfg['sx'][name].update(D_key=best[1], OA=best[2], SA=best[3])
        print(f'  {name}: D={best[1]} OA={best[2]} SA={best[3]}  heavy-key rec {best[4]:.4f}  light-key rec {best[5]:.4f}')

def calibrate_flotation(cfg, target_grade, lo=0.3, hi=3.0):
    """Scale gangue Rmax in the three cleaners to hit a concentrate REO grade."""
    base = copy.deepcopy(cfg['flotation']['kinetics'])
    vm = cfg['value_mineral']
    def grade(f):
        for st in range(4):
            for g, (k, r) in base[st].items():
                if g != vm: cfg['flotation']['kinetics'][st][g] = (k, round(min(r * f, 0.99), 4))
        fs = solve(cfg); return report(fs, cfg)['conc_REO'], fs
    a, b = lo, hi
    for _ in range(25):
        m = 0.5 * (a + b); gm, _ = grade(m)
        if gm > target_grade: a = m
        else: b = m
    return grade(0.5 * (a + b))

if __name__ == '__main__':
    out = {}
    for name in sys.argv[1:]:
        cfg = copy.deepcopy(ROUTES[name]); print('====', name)
        if name == 'neodymium_monazite_mtweld':
            g, fs = calibrate_flotation(cfg, 0.40)
            print('  flotation gangue Rmax scaled; concentrate REO', round(g, 4))
        fs = size_reagents(cfg)
        tune_sx(cfg, fs)
        fs = size_reagents(cfg)
        r = report(fs, cfg)
        print('  conc REO', round(r['conc_REO'],4), 'feed REO', round(r['feed_REO'],4), 'Nd to metal %', round(r['nd_balance'][-1][2],2), 'oxide purity', round(r['nd2o3_purity'],5))
        out[name] = {'flotation': cfg.get('flotation'), 'sx': cfg['sx'], 'reagents': cfg['reagents']}
    json.dump(out, open(os.path.join(_H, 'tuned_ree.json'), 'w'), indent=1)
