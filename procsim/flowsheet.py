"""Sequential-modular flowsheet solver with automatic tear-stream selection and
Wegstein-accelerated recycle convergence.

Flowsheet specification (JSON-compatible dict):
{
  "components": ["methane", ...],          # fluid species
  "minerals":   ["chalcopyrite", ...],     # optional particulate species
  "package":    "Peng-Robinson" | "Ideal",
  "units":   [{"id": "M1", "type": "Mixer", "params": {...}}],
  "streams": [{"id": "S1", "from": null | {"unit": "M1", "port": "out"},
               "to": null | {"unit": "F1", "port": "in"},
               "feed": {"T": K, "P": Pa, "flows": {...}, "solids": {...}}}]   # feed only when from is null
}
"""
import json
import math
import time
import numpy as np

from .components import MINERALS
from .stream import Context, Stream, NBINS
from .units import UNIT_TYPES, catalogue
from .numerics import ConvergenceError


class FlowsheetError(RuntimeError):
    pass


def _build_context(spec):
    comps = list(spec.get("components", []))
    mins = list(spec.get("minerals", []))
    if mins:
        for m in mins:
            for sp in MINERALS[m].leach_products:
                if sp not in comps:
                    comps.append(sp)
        if "dissolved other" not in comps:
            comps.append("dissolved other")
    return Context(comps, mins, spec.get("package", "Peng-Robinson"))


def _tarjan(nodes, edges):
    index, low, stack, on, out = {}, {}, [], set(), []
    counter = [0]

    def strong(v):
        index[v] = low[v] = counter[0]
        counter[0] += 1
        stack.append(v)
        on.add(v)
        for w in edges.get(v, ()):
            if w not in index:
                strong(w)
                low[v] = min(low[v], low[w])
            elif w in on:
                low[v] = min(low[v], index[w])
        if low[v] == index[v]:
            comp = []
            while True:
                w = stack.pop()
                on.discard(w)
                comp.append(w)
                if w == v:
                    break
            out.append(comp)

    for v in nodes:
        if v not in index:
            strong(v)
    return out[::-1]  # reverse topological order of Tarjan -> topological


class Flowsheet:
    def __init__(self, spec, tol=1e-6, max_iter=100, verbose=False):
        self.spec = spec
        self.tol, self.max_iter, self.verbose = tol, max_iter, verbose
        self.ctx = _build_context(spec)
        self.units = {}
        for u in spec.get("units", []):
            if u["type"] not in UNIT_TYPES:
                raise FlowsheetError(f"Unknown unit type '{u['type']}' for {u['id']}")
            self.units[u["id"]] = UNIT_TYPES[u["type"]](u["id"], u.get("params"))
        self.streams = {s["id"]: s for s in spec.get("streams", [])}
        self._validate()
        self.values = {}
        self.log = []

    # ---- topology ------------------------------------------------------
    def _validate(self):
        for sid, s in self.streams.items():
            for end in ("from", "to"):
                e = s.get(end)
                if e:
                    if e["unit"] not in self.units:
                        raise FlowsheetError(f"Stream {sid}: unknown unit {e['unit']}")
                    u = self.units[e["unit"]]
                    ports = u.OUTPORTS if end == "from" else u.INPORTS
                    if e["port"] not in ports:
                        raise FlowsheetError(f"Stream {sid}: unit {u.id} has no {end == 'from' and 'outlet' or 'inlet'} port '{e['port']}'")
            if not s.get("from") and not s.get("feed"):
                raise FlowsheetError(f"Stream {sid} has no source unit and no feed specification")
        for uid, u in self.units.items():
            if not any(s.get("to") and s["to"]["unit"] == uid for s in self.streams.values()):
                raise FlowsheetError(f"Unit {uid} has no inlet stream")

    def _edges(self):
        edges = {u: [] for u in self.units}
        for sid, s in self.streams.items():
            if s.get("from") and s.get("to"):
                edges[s["from"]["unit"]].append(s["to"]["unit"])
        return edges

    def _plan(self):
        """Return list of blocks: ('unit', uid) or ('loop', [uids in order], [tear stream ids])."""
        edges = self._edges()
        sccs = _tarjan(list(self.units), edges)
        plan = []
        for comp in sccs:
            cset = set(comp)
            internal = [sid for sid, s in self.streams.items()
                        if s.get("from") and s.get("to") and s["from"]["unit"] in cset and s["to"]["unit"] in cset]
            if len(comp) == 1 and not internal:
                plan.append(("unit", comp[0]))
                continue
            # start DFS from a unit that receives an external stream
            def has_external(u):
                return any(s.get("to") and s["to"]["unit"] == u and (not s.get("from") or s["from"]["unit"] not in cset)
                           for s in self.streams.values())
            starts = sorted(comp, key=lambda u: (not has_external(u), u))
            color, order, tears = {}, [], []
            sid_by_edge = {}
            for sid in internal:
                s = self.streams[sid]
                sid_by_edge.setdefault((s["from"]["unit"], s["to"]["unit"]), []).append(sid)

            def dfs(v):
                color[v] = 1
                for w in edges[v]:
                    if w not in cset:
                        continue
                    if color.get(w) == 1:
                        tears.extend(sid_by_edge[(v, w)])
                    elif w not in color:
                        dfs(w)
                color[v] = 2
                order.append(v)
            for st in starts:
                if st not in color:
                    dfs(st)
            # topological order of the loop with tear edges removed
            tear_set = set(tears)
            indeg = {u: 0 for u in comp}
            for sid in internal:
                if sid not in tear_set:
                    indeg[self.streams[sid]["to"]["unit"]] += 1
            ready = [u for u in starts if indeg[u] == 0]
            seq = []
            while ready:
                u = ready.pop(0)
                seq.append(u)
                for sid in internal:
                    s = self.streams[sid]
                    if sid not in tear_set and s["from"]["unit"] == u:
                        t = s["to"]["unit"]
                        indeg[t] -= 1
                        if indeg[t] == 0:
                            ready.append(t)
            plan.append(("loop", seq, sorted(set(tears))))
        return plan

    # ---- calculation ---------------------------------------------------
    def _inlets(self, uid):
        d = {}
        for sid, s in self.streams.items():
            if s.get("to") and s["to"]["unit"] == uid:
                st = self.values.get(sid)
                if st is not None:
                    d.setdefault(s["to"]["port"], []).append(st)
        return d

    def _run_unit(self, uid):
        u = self.units[uid]
        inl = self._inlets(uid)
        try:
            outs = u.run(inl, self.ctx)
        except Exception as e:  # report with context
            raise FlowsheetError(f"{uid} ({u.TYPE}): {e}") from e
        hin = sum(s.H for ss in inl.values() for s in ss)
        hout = sum(s.H for s in outs.values())
        u.results["Net heat input [kW]"] = (hout - hin) / 1e3
        u.net_heat = hout - hin
        for sid, s in self.streams.items():
            if s.get("from") and s["from"]["unit"] == uid:
                st = outs[s["from"]["port"]].copy(sid)
                self.values[sid] = st

    def _vec(self, sids):
        parts = []
        for sid in sids:
            s = self.values[sid]
            parts += [s.F, s.S.ravel(), [s.T / 100.0, s.P / 1e6]]
        return np.concatenate([np.asarray(p, float) for p in parts])

    def _set_vec(self, sids, v):
        k = 0
        nF, nS = len(self.ctx.names), len(self.ctx.minerals) * NBINS
        for sid in sids:
            s = self.values[sid]
            F = np.maximum(v[k:k + nF], 0.0); k += nF
            S = np.maximum(v[k:k + nS], 0.0).reshape(s.S.shape) if nS else s.S; k += nS
            T, P = v[k] * 100.0, v[k + 1] * 1e6; k += 2
            ns = Stream(self.ctx, T, P, F, S, sid)
            ns.flash_TP()
            self.values[sid] = ns

    def _init_tear(self, sid):
        s = self.streams[sid]
        if s.get("guess"):
            self.values[sid] = Stream.from_spec(self.ctx, s["guess"], sid)
        else:
            self.values[sid] = Stream(self.ctx, 298.15, 101325.0, name=sid).flash_TP()

    def solve(self):
        t0 = time.time()
        self.values = {}
        for sid, s in self.streams.items():
            if not s.get("from"):
                self.values[sid] = Stream.from_spec(self.ctx, s["feed"], sid)
        self.plan = self._plan()
        self.loops = []
        for block in self.plan:
            if block[0] == "unit":
                self._run_unit(block[1])
                continue
            _, seq, tears = block
            for sid in tears:
                if sid not in self.values:
                    self._init_tear(sid)
            x = self._vec(tears)
            xs_prev = gs_prev = None
            converged = False
            hist = []
            for it in range(1, self.max_iter + 1):
                for uid in seq:
                    self._run_unit(uid)
                g = self._vec(tears)  # g(x)
                scale = np.maximum(np.abs(g), 1e-6 * max(np.abs(g).max(), 1e-12) + 1e-12)
                err = float(np.max(np.abs(g - x) / scale))
                hist.append(err)
                if err < self.tol:
                    converged = True
                    break
                # Wegstein with bounded acceleration, after 2 direct substitutions
                if gs_prev is not None and it > 2:
                    dx = x - xs_prev
                    dg = g - gs_prev
                    with np.errstate(divide="ignore", invalid="ignore"):
                        sl = np.where(np.abs(dx) > 1e-14, dg / dx, 0.0)
                        q = np.where(np.abs(sl - 1) > 1e-12, sl / (sl - 1), 0.0)
                    q = np.clip(np.nan_to_num(q), -5.0, 0.0)
                    xn = q * x + (1 - q) * g
                else:
                    xn = g
                xs_prev, gs_prev = x, g
                x = xn
                self._set_vec(tears, x)
            self.loops.append({"units": seq, "tears": tears, "iterations": len(hist),
                               "converged": converged, "errors": hist})
            if not converged:
                self.log.append(f"Recycle loop {seq} did not converge in {self.max_iter} iterations (last error {hist[-1]:.2e})")
        self.elapsed = time.time() - t0
        return self

    # ---- results -------------------------------------------------------
    def balances(self):
        feeds = [sid for sid, s in self.streams.items() if not s.get("from")]
        prods = [sid for sid, s in self.streams.items() if not s.get("to")]
        m_in = sum(self.values[s].mass for s in feeds)
        m_out = sum(self.values[s].mass for s in prods)
        h_in = sum(self.values[s].H for s in feeds)
        h_out = sum(self.values[s].H for s in prods)
        q = sum(getattr(u, "net_heat", 0.0) for u in self.units.values())
        return {"mass_in_kg_s": m_in, "mass_out_kg_s": m_out,
                "mass_closure_pct": 100 * (m_out - m_in) / m_in if m_in else 0.0,
                "enthalpy_in_kW": h_in / 1e3, "enthalpy_out_kW": h_out / 1e3, "heat_added_kW": q / 1e3,
                "energy_closure_kW": (h_out - h_in - q) / 1e3}

    def results(self):
        return {
            "ok": True,
            "components": self.ctx.names,
            "minerals": self.ctx.minerals,
            "package": self.ctx.pkg.name if self.ctx.pkg else None,
            "streams": {sid: dict(self.values[sid].summary(), name=sid,
                                  source=(self.streams[sid].get("from") or {}).get("unit"),
                                  destination=(self.streams[sid].get("to") or {}).get("unit"))
                        for sid in self.streams if sid in self.values},
            "units": {uid: {"type": u.TYPE, "results": _clean(u.results)} for uid, u in self.units.items()},
            "loops": self.loops,
            "order": [b[1] if b[0] == "unit" else {"loop": b[1], "tears": b[2]} for b in self.plan],
            "balances": self.balances(),
            "warnings": self.log,
            "elapsed_s": self.elapsed,
        }

    def stream(self, sid):
        return self.values[sid]

    def table(self):
        """Plain-text stream table for terminals."""
        sids = list(self.values)
        rows = [("Stream", lambda s: s.name), ("T [C]", lambda s: f"{s.T - 273.15:.2f}"),
                ("P [bar]", lambda s: f"{s.P / 1e5:.3f}"), ("Vap frac", lambda s: f"{s.beta:.4f}"),
                ("Mass [kg/s]", lambda s: f"{s.mass:.4g}")]
        rows += [(f"{n} [mol/s]", (lambda i: lambda s: f"{s.F[i]:.4g}")(i)) for i, n in enumerate(self.ctx.names)]
        if self.ctx.minerals:
            rows += [(f"{m} [kg/s]", (lambda i: lambda s: f"{s.S[i].sum():.4g}")(i)) for i, m in enumerate(self.ctx.minerals)]
            rows += [("P80 [um]", lambda s: f"{s.psd_p80():.0f}" if s.solids_mass > 0 else "-")]
        w0 = max(len(r[0]) for r in rows)
        cols = [[f(self.values[sid]) for sid in sids] for _, f in rows]
        w = max(10, max(len(c) for col in cols for c in col))
        return "\n".join(r[0].ljust(w0) + " " + " ".join(c.rjust(w) for c in col) for r, col in zip(rows, cols))


def _clean(d):
    out = {}
    for k, v in d.items():
        if isinstance(v, (np.floating, np.integer)):
            v = v.item()
        if isinstance(v, float) and not math.isfinite(v):
            v = None
        out[k] = v
    return out


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        o = o.item()
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


def solve(spec, **kw):
    return Flowsheet(spec, **kw).solve()


def solve_json(text):
    """Entry point used by the web app: JSON spec in, JSON results out (never raises)."""
    try:
        fs = Flowsheet(json.loads(text)).solve()
        return json.dumps(_jsonable(fs.results()))
    except Exception as e:
        return json.dumps({"ok": False, "error": str(e)})


def catalogue_json():
    from .components import COMPONENTS, MINERALS
    return json.dumps(_jsonable({
        "units": catalogue(),
        "components": [{"name": c.name, "formula": c.formula, "mw": c.mw} for c in COMPONENTS.values() if not c.solute],
        "solutes": [c.name for c in COMPONENTS.values() if c.solute],
        "minerals": [{"name": m.name, "formula": m.formula, "density": m.density, "leach": m.leach_products} for m in MINERALS.values()],
        "packages": ["Peng-Robinson", "Ideal"],
    }))
