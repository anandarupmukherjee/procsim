// Parity check: the browser engine must reproduce the Python core on every example.
// Usage: python tests/make_reference.py && node tests/parity.js
const fs = require("fs"), path = require("path");
const PS = require("../web/engine.js");
const ex = path.join(__dirname, "..", "procsim", "data", "examples");
let worst = 0, fails = 0;
for (const name of ["ngl_fractionation", "methanol_loop", "copper_concentrator", "neodymium_bastnaesite", "neodymium_monazite_mtweld", "neodymium_eneabba_concentrate"]) {
  const spec = JSON.parse(fs.readFileSync(path.join(ex, name + ".json")));
  const ref = JSON.parse(fs.readFileSync(path.join(__dirname, "reference", name + ".json")));
  const t0 = Date.now(); const js = PS.solve(spec);
  if (!js.ok) { console.log(name, "JS ERROR", js.error); fails++; continue; }
  let w = 0, where = "";
  const cmp = (a, b, label, abs) => { const d = Math.abs(a - b) / Math.max(Math.abs(b), abs); if (d > w) { w = d; where = label; } };
  for (const [sid, s] of Object.entries(ref.streams)) {
    const j = js.streams[sid];
    cmp(j.T_K, s.T_K, sid + " T", 1); cmp(j.P_bar, s.P_bar, sid + " P", 1e-3);
    cmp(j.mass_flow_kg_s, s.mass_flow_kg_s, sid + " mass", 1e-6);
    cmp(j.enthalpy_kW, s.enthalpy_kW, sid + " H", 1);
    for (const [c, v] of Object.entries(s.flows_mol_s)) cmp(j.flows_mol_s[c], v, `${sid} ${c}`, 1e-4);
    if (s.solids_kg_s) for (const [m, v] of Object.entries(s.solids_kg_s)) cmp(j.solids_kg_s[m], v, `${sid} ${m}`, 1e-6);
  }
  for (const [uid, u] of Object.entries(ref.units)) for (const [k, v] of Object.entries(u.results))
    if (typeof v === "number") cmp(js.units[uid].results[k], v, `${uid} ${k}`, 1);
  console.log(`${name}: max relative difference ${w.toExponential(2)} (${where}), JS solve ${Date.now() - t0} ms`);
  worst = Math.max(worst, w);
  if (w > 1e-3) fails++;
}
console.log(fails ? `FAIL (${fails})` : `PASS, worst ${worst.toExponential(2)}`);
process.exit(fails ? 1 : 0);
