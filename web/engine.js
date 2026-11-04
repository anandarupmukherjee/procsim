/* procsim browser engine: a line-for-line JavaScript port of the Python core
 * (procsim/components.py, thermo.py, stream.py, units.py, flowsheet.py).
 * Same flowsheet JSON in, same results JSON out. Parity is checked by tests/parity.js. */
(function (root) {
  "use strict";
  const R = 8.314462618, T_REF = 298.15, K0 = 273.15, SQ2 = Math.SQRT2;

  // ------------------------------------------------------------------ numerics
  class ConvergenceError extends Error {}
  function brent(f, a, b, tol = 1e-9, maxiter = 200) {
    let fa = f(a), fb = f(b);
    if (fa === 0) return a;
    if (fb === 0) return b;
    if (fa * fb > 0) throw new ConvergenceError(`Root not bracketed: f(${a.toPrecision(4)})=${fa.toPrecision(4)}, f(${b.toPrecision(4)})=${fb.toPrecision(4)}`);
    let c = a, fc = fa, d = b - a, e = d;
    for (let i = 0; i < maxiter; i++) {
      if (fb * fc > 0) { c = a; fc = fa; d = e = b - a; }
      if (Math.abs(fc) < Math.abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
      const tol1 = 2 * 1e-15 * Math.abs(b) + 0.5 * tol, m = 0.5 * (c - b);
      if (Math.abs(m) <= tol1 || fb === 0) return b;
      if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
        const s = fb / fa; let p, q;
        if (a === c) { p = 2 * m * s; q = 1 - s; }
        else { const qq = fa / fc, r = fb / fc; p = s * (2 * m * qq * (qq - r) - (b - a) * (r - 1)); q = (qq - 1) * (r - 1) * (s - 1); }
        if (p > 0) q = -q; p = Math.abs(p);
        if (2 * p < Math.min(3 * m * q - Math.abs(tol1 * q), Math.abs(e * q))) { e = d; d = p / q; } else { d = e = m; }
      } else { d = e = m; }
      a = b; fa = fb;
      b = Math.abs(d) > tol1 ? b + d : b + (m >= 0 ? tol1 : -tol1);
      fb = f(b);
    }
    throw new ConvergenceError("Brent did not converge");
  }
  function bracketRoot(f, x0, step, lo, hi, maxiter = 60) {
    const f0 = f(x0); let a = x0, b = x0;
    for (let i = 0; i < maxiter; i++) {
      const a2 = Math.max(lo, a - step), b2 = Math.min(hi, b + step);
      if (f0 * f(a2) <= 0) return [a2, a];
      if (f0 * f(b2) <= 0) return [b, b2];
      a = a2; b = b2; step *= 1.6;
      if (a <= lo && b >= hi) break;
    }
    throw new ConvergenceError(`Could not bracket root in [${lo}, ${hi}]`);
  }
  const sum = (v) => v.reduce((s, x) => s + x, 0);
  const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };
  const zeros = (n) => new Array(n).fill(0);

  // ------------------------------------------------------------------ data
  function comp(name, formula, mw, tc, pcBar, omega, cp, hfKJ) {
    return { name, formula, mw, tc, pc: pcBar * 1e5, omega, cp, hf: hfKJ * 1e3, solute: false, cpLiq: 75 };
  }
  const COMPONENT_LIST = [
    comp("methane", "CH4", 16.043, 190.56, 45.99, 0.011, [19.25, 5.213e-2, 1.197e-5, -1.132e-8], -74.52),
    comp("ethane", "C2H6", 30.070, 305.32, 48.72, 0.099, [5.409, 1.781e-1, -6.938e-5, 8.713e-9], -83.82),
    comp("propane", "C3H8", 44.097, 369.83, 42.48, 0.152, [-4.224, 3.063e-1, -1.586e-4, 3.215e-8], -104.68),
    comp("n-butane", "C4H10", 58.123, 425.12, 37.96, 0.200, [9.487, 3.313e-1, -1.108e-4, -2.822e-9], -125.79),
    comp("n-pentane", "C5H12", 72.150, 469.70, 33.70, 0.252, [-3.626, 4.873e-1, -2.580e-4, 5.305e-8], -146.76),
    comp("n-hexane", "C6H14", 86.177, 507.60, 30.25, 0.301, [-4.413, 5.820e-1, -3.119e-4, 6.494e-8], -166.92),
    comp("n-heptane", "C7H16", 100.204, 540.20, 27.40, 0.350, [-5.146, 6.762e-1, -3.651e-4, 7.658e-8], -187.78),
    comp("benzene", "C6H6", 78.114, 562.05, 48.95, 0.210, [-33.92, 4.739e-1, -3.017e-4, 7.130e-8], 82.88),
    comp("toluene", "C7H8", 92.141, 591.75, 41.08, 0.264, [-24.35, 5.125e-1, -2.765e-4, 4.911e-8], 50.17),
    comp("water", "H2O", 18.015, 647.14, 220.64, 0.344, [32.24, 1.924e-3, 1.055e-5, -3.596e-9], -241.81),
    comp("nitrogen", "N2", 28.014, 126.20, 33.98, 0.037, [31.15, -1.357e-2, 2.680e-5, -1.168e-8], 0.0),
    comp("oxygen", "O2", 31.999, 154.58, 50.43, 0.022, [28.11, -3.680e-6, 1.746e-5, -1.065e-8], 0.0),
    comp("hydrogen", "H2", 2.016, 33.19, 13.13, -0.216, [27.14, 9.274e-3, -1.381e-5, 7.645e-9], 0.0),
    comp("carbon monoxide", "CO", 28.010, 132.85, 34.94, 0.045, [30.87, -1.285e-2, 2.789e-5, -1.272e-8], -110.53),
    comp("carbon dioxide", "CO2", 44.010, 304.21, 73.83, 0.224, [19.80, 7.344e-2, -5.602e-5, 1.715e-8], -393.51),
    comp("hydrogen sulfide", "H2S", 34.080, 373.40, 89.63, 0.090, [31.94, 1.436e-3, 2.432e-5, -1.176e-8], -20.63),
    comp("sulfur dioxide", "SO2", 64.064, 430.80, 78.84, 0.245, [23.85, 6.699e-2, -4.961e-5, 1.328e-8], -296.81),
    comp("ammonia", "NH3", 17.031, 405.40, 113.53, 0.257, [27.31, 2.383e-2, 1.707e-5, -1.185e-8], -45.90),
    comp("methanol", "CH3OH", 32.042, 512.64, 80.97, 0.565, [21.15, 7.092e-2, 2.587e-5, -2.852e-8], -200.94),
    comp("ethanol", "C2H5OH", 46.069, 513.92, 61.48, 0.649, [9.014, 2.141e-1, -8.390e-5, 1.373e-9], -234.95),
  ];
  for (const [n, f, mw] of [["Cu2+", "Cu2+", 63.546], ["Zn2+", "Zn2+", 65.38], ["Pb2+", "Pb2+", 207.2], ["Fe2+", "Fe2+", 55.845],
    ["Ni2+", "Ni2+", 58.693], ["Li+", "Li+", 6.94], ["dissolved other", "X", 100.0],
    ["La3+", "La3+", 138.905], ["Ce3+", "Ce3+", 140.116], ["Pr3+", "Pr3+", 140.908], ["Nd3+", "Nd3+", 144.242],
    ["Sm+HREE3+", "SmX3+", 150.36], ["Th4+", "Th4+", 232.038], ["Ca2+", "Ca2+", 40.078], ["reagent", "R", 60.0]])
    COMPONENT_LIST.push({ name: n, formula: f, mw, tc: 0, pc: 0, omega: 0, cp: [0, 0, 0, 0], hf: 0, solute: true, cpLiq: 75 });
  const COMPONENTS = Object.fromEntries(COMPONENT_LIST.map((c) => [c.name, c]));

  function hIg(c, T) {
    const [A, B, C, D] = c.cp;
    const F = (t) => A * t + B * t * t / 2 + C * t ** 3 / 3 + D * t ** 4 / 4;
    return c.hf + F(T) - F(T_REF);
  }
  const MINERAL_LIST = [
    ["chalcopyrite", "CuFeS2", 183.53, 4200, 540, { "Cu2+": 1 }, 10.0],
    ["pyrite", "FeS2", 119.98, 5010, 520, { "Fe2+": 1 }, 9.0],
    ["sphalerite", "ZnS", 97.44, 4050, 470, { "Zn2+": 1 }, 11.0],
    ["galena", "PbS", 239.27, 7500, 210, { "Pb2+": 1 }, 10.0],
    ["pentlandite", "Ni4.5Fe4.5S8", 771.9, 4800, 600, { "Ni2+": 4.5 }, 12.0],
    ["spodumene", "LiAlSi2O6", 186.09, 3150, 850, { "Li+": 1 }, 13.0],
    ["magnetite", "Fe3O4", 231.53, 5150, 650, { "Fe2+": 3 }, 10.0],
    ["hematite", "Fe2O3", 159.69, 5260, 650, { "Fe2+": 2 }, 12.0],
    ["quartz", "SiO2", 60.08, 2650, 740, {}, 16.0],
    ["bastnaesite", "(Ce,La,Nd,Pr)CO3F", 219.81, 4950, 600, { "La3+": 0.3514, "Ce3+": 0.4774, "Pr3+": 0.0415, "Nd3+": 0.1171, "Sm+HREE3+": 0.0126, "Th4+": 0.002 }, 12.0],
    ["monazite", "(Ce,La,Nd,Th)PO4", 241.59, 5150, 560, { "La3+": 0.2323, "Ce3+": 0.4378, "Pr3+": 0.048, "Nd3+": 0.1708, "Sm+HREE3+": 0.0511, "Th4+": 0.06 }, 13.0],
    ["monazite-xenotime", "(Ce,La,Nd,Y,Th)PO4", 242.34, 5000, 560, { "La3+": 0.2087, "Ce3+": 0.4078, "Pr3+": 0.0421, "Nd3+": 0.1539, "Sm+HREE3+": 0.1276, "Th4+": 0.06 }, 13.0],
    ["apatite", "Ca5(PO4)3F", 504.3, 3190, 700, { "Ca2+": 5 }, 9.0],
    ["calcite", "CaCO3", 100.09, 2710, 820, { "Ca2+": 1 }, 10.0],
    ["barite", "BaSO4", 233.39, 4480, 460, {}, 6.0],
    ["La2(CO3)3", "La2(CO3)3", 457.84, 2600, 600, { "La3+": 2 }, 5.0],
    ["Ce2(CO3)3", "Ce2(CO3)3", 460.26, 2600, 600, { "Ce3+": 2 }, 5.0],
    ["Pr2(CO3)3", "Pr2(CO3)3", 461.84, 2600, 600, { "Pr3+": 2 }, 5.0],
    ["Nd2(CO3)3", "Nd2(CO3)3", 468.51, 2600, 600, { "Nd3+": 2 }, 5.0],
    ["Sm2(CO3)3", "Sm2(CO3)3", 480.75, 2600, 600, { "Sm+HREE3+": 2 }, 5.0],
    ["La2O3", "La2O3", 325.81, 6510, 330, { "La3+": 2 }, 8.0],
    ["CeO2", "CeO2", 172.11, 7220, 360, { "Ce3+": 1 }, 8.0],
    ["Pr6O11", "Pr6O11", 1021.44, 6900, 330, { "Pr3+": 6 }, 8.0],
    ["Nd2O3", "Nd2O3", 336.48, 7240, 330, { "Nd3+": 2 }, 8.0],
    ["Sm2O3", "Sm2O3", 348.72, 8350, 330, { "Sm+HREE3+": 2 }, 8.0],
    ["ThO2", "ThO2", 264.04, 10000, 240, { "Th4+": 1 }, 8.0],
    ["Nd metal", "Nd", 144.242, 7010, 190, { "Nd3+": 1 }, 8.0],
  ].map(([name, formula, mw, density, cp, leach, wi]) => ({ name, formula, mw, density, cp, leach, workIndex: wi }));
  const MINERALS = Object.fromEntries(MINERAL_LIST.map((m) => [m.name, m]));

  const HC = ["methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "n-heptane", "benzene", "toluene"];
  const KIJ = {};
  const kkey = (a, b) => [a, b].sort().join("|");
  for (const h of HC) {
    KIJ[kkey("water", h)] = 0.48; KIJ[kkey("carbon dioxide", h)] = 0.12;
    KIJ[kkey("nitrogen", h)] = 0.08; KIJ[kkey("hydrogen sulfide", h)] = 0.07;
  }
  KIJ[kkey("water", "carbon dioxide")] = 0.12; KIJ[kkey("water", "nitrogen")] = 0.48;
  const kij = (a, b) => KIJ[kkey(a, b)] || 0;
  function component(n) { if (!COMPONENTS[n]) throw new Error(`Unknown component '${n}'`); return COMPONENTS[n]; }
  function mineral(n) { if (!MINERALS[n]) throw new Error(`Unknown mineral '${n}'`); return MINERALS[n]; }

  // ------------------------------------------------------------------ thermo
  function rachfordRice(z0, K0_) {
    const z = [], K = [];
    for (let i = 0; i < z0.length; i++) if (z0[i] > 1e-300) { z.push(z0[i]); K.push(K0_[i]); }
    const Km = Math.min(...K), Kx = Math.max(...K);
    if (Kx <= 1 + 1e-12 || Km >= 1 - 1e-12) return null;
    const lo = 1 / (1 - Kx), hi = 1 / (1 - Km);
    const g = (b) => { let s = 0; for (let i = 0; i < z.length; i++) s += z[i] * (K[i] - 1) / (1 + b * (K[i] - 1)); return s; };
    let eps = 1e-10 * (hi - lo), a = lo + eps, c = hi - eps;
    while (g(a) < 0) { eps *= 1e-3; const a2 = lo + eps; if (a2 === a || a2 <= lo) return a; a = a2; }
    eps = 1e-10 * (hi - lo);
    while (g(c) > 0) { eps *= 1e-3; const c2 = hi - eps; if (c2 === c || c2 >= hi) return c; c = c2; }
    return brent(g, a, c, 1e-14);
  }
  const normalize = (v) => { const s = sum(v); return v.map((x) => x / s); };

  class ThermoPackage {
    constructor(names) {
      this.names = names.slice(); this.comps = names.map(component); this.n = names.length;
      this.tc = this.comps.map((c) => c.tc); this.pc = this.comps.map((c) => c.pc); this.w = this.comps.map((c) => c.omega);
    }
    hIgMix(T, x) { let h = 0; for (let i = 0; i < this.n; i++) h += x[i] * hIg(this.comps[i], T); return h; }
    wilsonK(T, P) { return this.tc.map((tc, i) => this.pc[i] / P * Math.exp(5.373 * (1 + this.w[i]) * (1 - tc / T))); }
    flashTP(T, P, z) {
      z = normalize(z);
      if (this.n === 1) return this._flashPure(T, P, z);
      let K = this.wilsonK(T, P), beta = 0.5, x = z.slice(), y = z.slice();
      for (let it = 0; it < 300; it++) {
        const b = rachfordRice(z, K);
        if (b === null) {
          beta = Math.min(...K) >= 1 ? 1 : 0; x = z.slice(); y = z.slice();
          const Kn = this.K(T, P, z, z);
          if (Kn.every((k, i) => Math.abs(k - K[i]) <= 1e-8 * Math.abs(K[i]))) break;
          K = Kn;
          if (rachfordRice(z, K) === null) break;
          continue;
        }
        beta = b;
        x = z.map((zi, i) => zi / (1 + beta * (K[i] - 1)));
        y = x.map((xi, i) => K[i] * xi);
        x = normalize(x); y = normalize(y);
        const Kn = this.K(T, P, x, y);
        let err = 0; for (let i = 0; i < this.n; i++) err = Math.max(err, Math.abs(Math.log(Kn[i] / K[i])));
        K = Kn;
        if (err < 1e-10) break;
        if (Math.max(...K.map((k) => Math.abs(Math.log(k)))) < 1e-4) {
          beta = this._prefersVapour(T, P, z) ? 1 : 0; x = z.slice(); y = z.slice(); break;
        }
      }
      if (beta <= 0) return { T, P, beta: 0, x: z.slice(), y: z.slice(), phase: "L" };
      if (beta >= 1) return { T, P, beta: 1, x: z.slice(), y: z.slice(), phase: "V" };
      return { T, P, beta, x, y, phase: "VL" };
    }
    _prefersVapour(T, P, z) { return dot(z, this.tc) < T; }
    _flashPure(T, P, z) {
      let Tb = null;
      try { Tb = this.bubbleT(z, P); } catch (e) { Tb = null; }
      if (Tb === null || T > Tb) return { T, P, beta: 1, x: z, y: z, phase: "V" };
      return { T, P, beta: 0, x: z, y: z, phase: "L" };
    }
    enthalpy(fr) {
      if (fr.phase === "V") return this.h(fr.T, fr.P, fr.y, "V");
      if (fr.phase === "L") return this.h(fr.T, fr.P, fr.x, "L");
      return fr.beta * this.h(fr.T, fr.P, fr.y, "V") + (1 - fr.beta) * this.h(fr.T, fr.P, fr.x, "L");
    }
    _boundaryT(z, P, target, T0) {
      const f = (T) => this.flashTP(T, P, z).beta - target;
      const [a, b] = bracketRoot(f, T0, Math.max(2, 0.02 * T0), 20, 3000);
      return brent(f, a, b, 1e-5);
    }
    _wilsonT(z, P, bubble) {
      const gw = (T) => { const K = this.wilsonK(T, P); return bubble ? Math.log(dot(z, K)) : -Math.log(sum(z.map((zi, i) => zi / K[i]))); };
      try { const [a, b] = bracketRoot(gw, 300, 20, 20, 2500); return brent(gw, a, b, 1e-6); } catch (e) { return 300; }
    }
    bubbleT(x, P, T0) {
      x = normalize(x);
      if (x.filter((v) => v > 1e-12).length > 1) return this._boundaryT(x, P, 1e-7, T0 || this._wilsonT(x, P, true));
      const g = (T) => {
        let K = this.wilsonK(T, P), y;
        for (let k = 0; k < 50; k++) {
          const s = dot(K, x); y = x.map((xi, i) => K[i] * xi / s);
          const Kn = this.K(T, P, x, y);
          const done = Kn.every((v, i) => Math.abs(v - K[i]) < 1e-10); K = Kn; if (done) break;
        }
        return Math.log(dot(K, x));
      };
      return this._solveT(g, x, P, T0);
    }
    dewT(y, P, T0) {
      y = normalize(y);
      if (y.filter((v) => v > 1e-12).length > 1) return this._boundaryT(y, P, 1 - 1e-7, T0 || this._wilsonT(y, P, false));
      const g = (T) => {
        let K = this.wilsonK(T, P), x;
        for (let k = 0; k < 50; k++) {
          const s = sum(y.map((yi, i) => yi / K[i])); x = y.map((yi, i) => yi / K[i] / s);
          const Kn = this.K(T, P, x, y);
          const done = Kn.every((v, i) => Math.abs(v - K[i]) < 1e-10); K = Kn; if (done) break;
        }
        return -Math.log(sum(y.map((yi, i) => yi / K[i])));
      };
      return this._solveT(g, y, P, T0);
    }
    _solveT(g, z, P, T0) {
      if (!T0) {
        const gw = (T) => Math.log(dot(z, this.wilsonK(T, P)));
        try { const [a, b] = bracketRoot(gw, 300, 20, 20, 2500); T0 = brent(gw, a, b, 1e-6); } catch (e) { T0 = 300; }
      }
      const [a, b] = bracketRoot(g, T0, Math.max(2, 0.02 * T0), 20, 3000);
      return brent(g, a, b, 1e-7);
    }
  }

  class Ideal extends ThermoPackage {
    get name() { return "Ideal"; }
    psat(T) {
      return this.tc.map((tc, i) => {
        const Tr = T / tc;
        const f0 = 5.92714 - 6.09648 / Tr - 1.28862 * Math.log(Tr) + 0.169347 * Tr ** 6;
        const f1 = 15.2518 - 15.6875 / Tr - 13.4721 * Math.log(Tr) + 0.43577 * Tr ** 6;
        return this.pc[i] * Math.exp(Math.min(f0 + this.w[i] * f1, 50));
      });
    }
    K(T, P) { return this.psat(T).map((p) => p / P); }
    wilsonK(T, P) { return this.K(T, P); }
    dhvap(T) {
      return this.tc.map((tc, i) => { const Tr = Math.min(T / tc, 1); return R * tc * (7.08 * (1 - Tr) ** 0.354 + 10.95 * this.w[i] * (1 - Tr) ** 0.456); });
    }
    h(T, P, x, phase) { const hig = this.hIgMix(T, x); return phase === "V" ? hig : hig - dot(x, this.dhvap(T)); }
  }

  function cubicRoots(a2, a1, a0) {
    // real roots of Z^3 + a2 Z^2 + a1 Z + a0
    const q = (3 * a1 - a2 * a2) / 9, r = (9 * a2 * a1 - 27 * a0 - 2 * a2 ** 3) / 54;
    const disc = q ** 3 + r * r, off = -a2 / 3;
    if (disc > 0) {
      const s = Math.cbrt(r + Math.sqrt(disc)), t = Math.cbrt(r - Math.sqrt(disc));
      return [s + t + off];
    }
    const th = Math.acos(Math.max(-1, Math.min(1, r / Math.sqrt(-(q ** 3) || 1e-300)))), m = 2 * Math.sqrt(-q);
    return [m * Math.cos(th / 3) + off, m * Math.cos((th + 2 * Math.PI) / 3) + off, m * Math.cos((th + 4 * Math.PI) / 3) + off];
  }

  class PengRobinson extends ThermoPackage {
    get name() { return "Peng-Robinson"; }
    constructor(names) {
      super(names);
      this.kappa = this.w.map((w) => 0.37464 + 1.54226 * w - 0.26992 * w * w);
      this.ac = this.tc.map((tc, i) => 0.45724 * R * R * tc * tc / this.pc[i]);
      this.b = this.tc.map((tc, i) => 0.07780 * R * tc / this.pc[i]);
      this.k = this.names.map((a) => this.names.map((b) => kij(a, b)));
    }
    _a(T) {
      const n = this.n, ai = [], dai = [];
      for (let i = 0; i < n; i++) {
        const sTr = Math.sqrt(T / this.tc[i]), m = 1 + this.kappa[i] * (1 - sTr);
        ai.push(this.ac[i] * m * m);
        dai.push(-this.ac[i] * this.kappa[i] * m / Math.sqrt(T * this.tc[i]));
      }
      const sa = ai.map(Math.sqrt), aij = [], daij = [];
      for (let i = 0; i < n; i++) {
        aij.push([]); daij.push([]);
        for (let j = 0; j < n; j++) {
          const kk = 1 - this.k[i][j];
          aij[i].push(kk * sa[i] * sa[j]);
          daij[i].push(kk * 0.5 * (dai[i] / sa[i] * sa[j] + sa[i] * dai[j] / sa[j]));
        }
      }
      return [aij, daij];
    }
    _state(T, P, x, phase) {
      const [aij, daij] = this._a(T), n = this.n;
      let am = 0, dam = 0;
      const ax = zeros(n);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { ax[i] += aij[i][j] * x[j]; dam += x[i] * x[j] * daij[i][j]; }
      am = dot(x, ax);
      const bm = dot(x, this.b), A = am * P / (R * T) ** 2, B = bm * P / (R * T);
      let roots = cubicRoots(-(1 - B), A - 3 * B * B - 2 * B, -(A * B - B * B - B ** 3)).filter((z) => z > B);
      if (!roots.length) throw new ConvergenceError("No physical PR root");
      let Z;
      if (phase === "L") Z = Math.min(...roots);
      else if (phase === "V") Z = Math.max(...roots);
      else {
        const g = (Z) => Z - 1 - Math.log(Z - B) - A / (2 * SQ2 * B) * Math.log((Z + (1 + SQ2) * B) / (Z + (1 - SQ2) * B));
        Z = roots.reduce((best, z) => (g(z) < g(best) ? z : best), roots[0]);
      }
      return { ax, am, dam, bm, A, B, Z };
    }
    lnphi(T, P, x, phase) {
      const { ax, am, bm, A, B, Z } = this._state(T, P, x, phase);
      const L = Math.log((Z + (1 + SQ2) * B) / (Z + (1 - SQ2) * B));
      return this.b.map((bi0, i) => {
        const bi = bi0 / bm;
        return bi * (Z - 1) - Math.log(Z - B) - A / (2 * SQ2 * B) * (2 * ax[i] / am - bi) * L;
      });
    }
    K(T, P, x, y) { const l = this.lnphi(T, P, x, "L"), v = this.lnphi(T, P, y, "V"); return l.map((li, i) => Math.exp(li - v[i])); }
    h(T, P, x, phase) {
      const { am, dam, bm, B, Z } = this._state(T, P, x, phase);
      const L = Math.log((Z + (1 + SQ2) * B) / (Z + (1 - SQ2) * B));
      return this.hIgMix(T, x) + R * T * (Z - 1) + (T * dam - am) / (2 * SQ2 * bm) * L;
    }
    enthalpy(fr) { return fr.phase === "VL" ? super.enthalpy(fr) : this.h(fr.T, fr.P, fr.x, "auto"); }
    _prefersVapour(T, P, z) { return this._state(T, P, z, "auto").Z > 0.3; }
  }
  const PACKAGES = { "Peng-Robinson": PengRobinson, PR: PengRobinson, Ideal };

  // ------------------------------------------------------------------ PSD
  const SIEVES = Array.from({ length: 24 }, (_, i) => 50000 * 2 ** (-i / 2));
  const NBINS = SIEVES.length + 1;
  const MID = SIEVES.map((s) => s * 2 ** 0.25).concat([SIEVES[23] / 2]);
  function passingAtSieves(m) {
    const tot = sum(m); if (tot <= 0) return zeros(SIEVES.length);
    const cum = zeros(m.length); let s = 0;
    for (let i = m.length - 1; i >= 0; i--) { s += m[i]; cum[i] = s; }
    return cum.slice(1).map((c) => c / tot);
  }
  function pSize(m, frac = 0.8) {
    if (sum(m) <= 0) return NaN;
    const pas = passingAtSieves(m);
    if (pas[0] < frac) return SIEVES[0];
    for (let i = 0; i < pas.length - 1; i++) {
      const hi = pas[i], lo = pas[i + 1];
      if (hi >= frac && frac >= lo) {
        if (hi === lo) return SIEVES[i + 1];
        const t = (frac - lo) / (hi - lo);
        return Math.exp(Math.log(SIEVES[i + 1]) + t * (Math.log(SIEVES[i]) - Math.log(SIEVES[i + 1])));
      }
    }
    return SIEVES[SIEVES.length - 1];
  }
  function fracFromPassing(pas) {
    const f = zeros(NBINS);
    f[0] = 1 - pas[0];
    for (let i = 1; i < NBINS - 1; i++) f[i] = pas[i - 1] - pas[i];
    f[NBINS - 1] = pas[pas.length - 1];
    return f.map((v) => Math.max(v, 0));
  }
  function rrsbFractions(p80, spread = 1) {
    const d63 = p80 / Math.log(5) ** (1 / spread);
    const f = fracFromPassing(SIEVES.map((s) => 1 - Math.exp(-((s / d63) ** spread))));
    const t = sum(f); return f.map((v) => v / t);
  }
  const colSum = (S) => { const t = zeros(NBINS); for (const row of S) for (let j = 0; j < NBINS; j++) t[j] += row[j]; return t; };

  // ------------------------------------------------------------------ context & stream
  class Context {
    constructor(components, minerals = [], pkg = "Peng-Robinson") {
      this.names = components.slice(); this.comps = this.names.map(component);
      this.vol = []; this.sol = [];
      this.comps.forEach((c, i) => (c.solute ? this.sol : this.vol).push(i));
      this.mw = this.comps.map((c) => c.mw / 1000);
      this.minerals = minerals.slice(); this.mins = this.minerals.map(mineral);
      this.minCp = this.mins.map((m) => m.cp);
      if (!PACKAGES[pkg]) throw new Error(`Unknown property package '${pkg}'`);
      this.pkg = this.vol.length ? new PACKAGES[pkg](this.vol.map((i) => this.names[i])) : null;
    }
    index(n) { const i = this.names.indexOf(n); if (i < 0) throw new Error(`Species '${n}' is not in the component list`); return i; }
    mindex(n) { const i = this.minerals.indexOf(n); if (i < 0) throw new Error(`Mineral '${n}' is not in the mineral list`); return i; }
  }

  class Stream {
    constructor(ctx, T = 298.15, P = 101325, F = null, S = null, name = "") {
      this.ctx = ctx; this.name = name; this.T = T; this.P = P;
      this.F = F ? F.slice() : zeros(ctx.names.length);
      this.S = S ? S.map((r) => r.slice()) : ctx.minerals.map(() => zeros(NBINS));
      this.beta = 0; this.x = null; this.y = null; this.phase = "L"; this.H = 0;
    }
    static fromSpec(ctx, spec, name = "") {
      const s = new Stream(ctx, spec.T ?? 298.15, spec.P ?? 101325, null, null, name);
      for (const [k, v] of Object.entries(spec.flows || {})) s.F[ctx.index(k)] = +v;
      for (let [k, v] of Object.entries(spec.solids || {})) {
        if (typeof v === "number") v = { mass: v };
        const f = rrsbFractions(+(v.p80 ?? 10000), +(v.spread ?? 1));
        s.S[ctx.mindex(k)] = f.map((x) => x * +(v.mass || 0));
      }
      return s.flashTP();
    }
    copy(name) {
      const s = new Stream(this.ctx, this.T, this.P, this.F, this.S, name ?? this.name);
      Object.assign(s, { beta: this.beta, x: this.x, y: this.y, phase: this.phase, H: this.H });
      return s;
    }
    get nVol() { let s = 0; for (const i of this.ctx.vol) s += this.F[i]; return s; }
    get molarFlow() { return sum(this.F); }
    get fluidMass() { return dot(this.F, this.ctx.mw); }
    get solidsMass() { let s = 0; for (const r of this.S) s += sum(r); return s; }
    get mass() { return this.fluidMass + this.solidsMass; }
    _hNonvol(T) {
      let h = 0;
      for (const i of this.ctx.sol) { const c = this.ctx.comps[i]; h += this.F[i] * (c.hf + c.cpLiq * (T - T_REF)); }
      this.S.forEach((r, i) => { h += sum(r) * this.ctx.minCp[i] * (T - T_REF); });
      return h;
    }
    _flashVol(T) {
      const n = this.nVol, pkg = this.ctx.pkg;
      if (n <= 1e-15 || !pkg) return [null, 0];
      const fr = pkg.flashTP(T, this.P, this.ctx.vol.map((i) => this.F[i] / n));
      return [fr, n * pkg.enthalpy(fr)];
    }
    flashTP(T, P) {
      if (T !== undefined) this.T = T;
      if (P !== undefined) this.P = P;
      const [fr, hv] = this._flashVol(this.T);
      if (!fr) { this.beta = 0; this.x = this.y = null; this.phase = "L"; }
      else { this.beta = fr.beta; this.x = fr.x; this.y = fr.y; this.phase = fr.phase; }
      this.H = hv + this._hNonvol(this.T);
      return this;
    }
    enthalpyAt(T) { const [, hv] = this._flashVol(T); return hv + this._hNonvol(T); }
    flashPH(H, P, Tguess) {
      if (P !== undefined && P !== null) this.P = P;
      const T0 = Tguess || this.T || 300, scale = Math.max(Math.abs(H), 1);
      const f = (T) => (this.enthalpyAt(T) - H) / scale;
      try { const [a, b] = bracketRoot(f, T0, 5, 30, 3000); return this.flashTP(brent(f, a, b, 1e-6)); }
      catch (e) { throw new ConvergenceError(`PH flash failed for stream ${this.name}: ${e.message}`); }
    }
    splitPhases() {
      const ctx = this.ctx, V = new Stream(ctx, this.T, this.P), L = new Stream(ctx, this.T, this.P, this.F, this.S);
      if (this.beta > 0 && this.x) {
        const n = this.nVol;
        ctx.vol.forEach((i, k) => {
          if (this.beta >= 1) { V.F[i] = this.F[i]; L.F[i] = 0; }
          else { V.F[i] = this.beta * n * this.y[k]; L.F[i] = (1 - this.beta) * n * this.x[k]; }
        });
      }
      return [V.flashTP(), L.flashTP()];
    }
    p80() { const t = colSum(this.S); return sum(t) > 0 ? pSize(t) : null; }
    summary() {
      const ctx = this.ctx, mf = this.molarFlow;
      const d = {
        name: this.name, T_K: this.T, T_C: this.T - K0, P_bar: this.P / 1e5, vapour_fraction: this.beta, phase: this.phase,
        molar_flow_mol_s: mf, mass_flow_kg_s: this.mass, enthalpy_kW: this.H / 1e3,
        flows_mol_s: Object.fromEntries(ctx.names.map((n, i) => [n, this.F[i]])),
        mole_fractions: Object.fromEntries(ctx.names.map((n, i) => [n, mf > 0 ? this.F[i] / mf : 0])),
      };
      if (ctx.minerals.length) {
        d.solids_kg_s = Object.fromEntries(ctx.minerals.map((m, i) => [m, sum(this.S[i])]));
        d.solids_total_kg_s = this.solidsMass;
        d.solids_pct_w = this.mass > 0 ? 100 * this.solidsMass / this.mass : 0;
        d.P80_um = this.p80();
        d.psd_passing = this.solidsMass > 0 ? passingAtSieves(colSum(this.S)) : null;
      }
      return d;
    }
  }

  // ------------------------------------------------------------------ units
  const UNIT_TYPES = {};
  const P_ = (name, kind, def, units = "", help = "", options, when) => {
    const d = { name, kind, default: def, units, help };
    if (options) d.options = options; if (when) d.when = when; return d;
  };
  function mix(streams, ctx, P = null, name = "mix") {
    streams = streams.filter(Boolean);
    if (!streams.length) throw new Error("No inlet streams");
    if (streams.length === 1 && P === null) return streams[0].copy(name);
    const out = new Stream(ctx, 298.15, 101325, null, null, name);
    for (const s of streams) {
      s.F.forEach((v, i) => (out.F[i] += v));
      s.S.forEach((r, i) => r.forEach((v, j) => (out.S[i][j] += v)));
    }
    const H = sum(streams.map((s) => s.H)), Pm = P === null ? Math.min(...streams.map((s) => s.P)) : P;
    const Tg = sum(streams.map((s) => s.T * s.mass)) / Math.max(sum(streams.map((s) => s.mass)), 1e-30);
    out.P = Pm; out.T = Tg;
    if (out.molarFlow <= 0 && out.solidsMass <= 0) return out.flashTP();
    return out.flashPH(H, Pm, Tg);
  }
  function parseTable(text, keys) {
    const out = {};
    for (let line of (text || "").replace(/;/g, "\n").split("\n")) {
      line = line.trim(); if (!line || line.startsWith("#")) continue;
      let key, rest;
      if (line.includes(":")) { const k = line.indexOf(":"); key = line.slice(0, k); rest = line.slice(k + 1); }
      else { const p = line.split(/\s+/); key = p[0]; rest = p.slice(1).join(" "); }
      key = key.trim();
      const vals = {}, named = [...rest.matchAll(/([A-Za-z_]+)\s*=\s*([-+0-9.eE]+)/g)];
      if (named.length) named.forEach((m) => (vals[m[1]] = parseFloat(m[2])));
      else (rest.match(/[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?/g) || []).forEach((v, i) => { if (i < keys.length) vals[keys[i]] = parseFloat(v); });
      out[key] = vals;
    }
    return out;
  }
  const first = (o, d = 0) => { const v = Object.values(o || {}); return v.length ? v[0] : d; };

  class Unit {
    constructor(id, params) {
      this.id = id; this.p = {};
      for (const d of this.constructor.PARAMS) this.p[d.name] = d.default;
      Object.assign(this.p, params || {}); this.results = {};
    }
    static schema() {
      const c = this;
      return { type: c.TYPE, category: c.CATEGORY, label: c.LABEL, description: c.DESCRIPTION, inports: c.INPORTS, outports: c.OUTPORTS, params: c.PARAMS };
    }
    num(k) { const v = this.p[k]; return v === null || v === undefined || v === "" ? null : +v; }
    feed(inlets, ctx, port = "in") { return mix(inlets[port] || [], ctx, null, `${this.id}.${port}`); }
  }
  function reg(cls, TYPE, CATEGORY, LABEL, DESCRIPTION, INPORTS, OUTPORTS, PARAMS) {
    Object.assign(cls, { TYPE, CATEGORY, LABEL, DESCRIPTION, INPORTS, OUTPORTS, PARAMS }); UNIT_TYPES[TYPE] = cls;
  }
  const clamp = (v, a, b) => Math.min(Math.max(v, a), b);

  class Mixer extends Unit {
    run(inl, ctx) {
      const P = this.num("P_out"), out = mix(inl.in || [], ctx, P === null ? null : P * 1e5);
      this.results = { "Outlet T [C]": out.T - K0 }; return { out };
    }
  }
  reg(Mixer, "Mixer", "Fluid", "Mixer", "Adiabatic mixing of any number of streams.", ["in"], ["out"],
    [P_("P_out", "number", "", "bar", "Outlet pressure (blank = lowest inlet)")]);

  class Splitter extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), a = clamp(this.num("frac1"), 0, 1), o1 = f.copy(), o2 = f.copy();
      o1.F = f.F.map((v) => v * a); o1.S = f.S.map((r) => r.map((v) => v * a)); o1.H = f.H * a;
      o2.F = f.F.map((v) => v * (1 - a)); o2.S = f.S.map((r) => r.map((v) => v * (1 - a))); o2.H = f.H * (1 - a);
      this.results = { "Split to out1": a }; return { out1: o1, out2: o2 };
    }
  }
  reg(Splitter, "Splitter", "Fluid", "Splitter", "Splits a stream into two of identical composition and state.", ["in"], ["out1", "out2"],
    [P_("frac1", "number", 0.5, "-", "Fraction of the inlet sent to out1")]);

  class Heater extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), out = f.copy(`${this.id}.out`), P = f.P - (this.num("dP") || 0) * 1e5, mode = this.p.mode;
      if (mode === "temperature") out.flashTP(this.num("T_out") + K0, P);
      else if (mode === "duty") out.flashPH(f.H + this.num("duty") * 1e3, P);
      else {
        const vf = clamp(this.num("vf"), 0, 1); out.P = P;
        const n = out.nVol; if (n <= 0) throw new Error("Vapour fraction spec needs volatile components");
        const z = ctx.vol.map((i) => out.F[i] / n); let T;
        if (vf <= 0) T = ctx.pkg.bubbleT(z, P);
        else if (vf >= 1) T = ctx.pkg.dewT(z, P);
        else { const Tb = ctx.pkg.bubbleT(z, P), Td = ctx.pkg.dewT(z, P); T = brent((t) => ctx.pkg.flashTP(t, P, z).beta - vf, Tb + 1e-6, Td - 1e-6, 1e-6); }
        out.flashTP(T, P);
      }
      this.results = { "Duty [kW]": (out.H - f.H) / 1e3, "Outlet T [C]": out.T - K0, "Outlet vapour fraction": out.beta };
      return { out };
    }
  }
  reg(Heater, "Heater", "Fluid", "Heater / cooler", "Sets outlet temperature (or applies a duty) with an optional pressure drop.", ["in"], ["out"], [
    P_("mode", "select", "temperature", "", "Specification", ["temperature", "duty", "vapour fraction"]),
    P_("T_out", "number", 100.0, "C", "Outlet temperature", null, { mode: "temperature" }),
    P_("duty", "number", 0.0, "kW", "Heat added (negative = cooling)", null, { mode: "duty" }),
    P_("vf", "number", 1.0, "-", "Outlet vapour fraction (0 = bubble, 1 = dew)", null, { mode: "vapour fraction" }),
    P_("dP", "number", 0.0, "bar", "Pressure drop")]);

  class Valve extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), out = f.copy(`${this.id}.out`).flashPH(f.H, this.num("P_out") * 1e5);
      this.results = { "Outlet T [C]": out.T - K0, "Outlet vapour fraction": out.beta }; return { out };
    }
  }
  reg(Valve, "Valve", "Fluid", "Valve", "Isenthalpic pressure let-down (Joule-Thomson).", ["in"], ["out"],
    [P_("P_out", "number", 1.01325, "bar", "Outlet pressure")]);

  class Flash extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), s = f.copy(`${this.id}.drum`), P = this.num("P") * 1e5;
      if (this.p.mode === "T and P") s.flashTP(this.num("T") + K0, P);
      else if (this.p.mode === "adiabatic") s.flashPH(f.H, P);
      else s.flashPH(f.H + this.num("duty") * 1e3, P);
      const [vapour, liquid] = s.splitPhases();
      this.results = { "Duty [kW]": (s.H - f.H) / 1e3, "T [C]": s.T - K0, "Vapour fraction": s.beta };
      return { vapour, liquid };
    }
  }
  reg(Flash, "Flash", "Fluid", "Flash drum", "Two-phase vapour-liquid separator. Solids and dissolved species leave with the liquid.", ["in"], ["vapour", "liquid"], [
    P_("mode", "select", "T and P", "", "Specification", ["T and P", "adiabatic", "duty"]),
    P_("T", "number", 25.0, "C", "Drum temperature", null, { mode: "T and P" }),
    P_("P", "number", 1.01325, "bar", "Drum pressure"),
    P_("duty", "number", 0.0, "kW", "Heat added", null, { mode: "duty" })]);

  class ComponentSplitter extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), tab = {};
      for (const [k, v] of Object.entries(parseTable(this.p.splits, ["frac"]))) tab[k] = v.frac ?? v.f ?? first(v, 0);
      const d = this.num("default") || 0, o1 = new Stream(ctx, f.T, f.P), o2 = new Stream(ctx, f.T, f.P);
      ctx.names.forEach((n, i) => { const a = tab[n] ?? d; o1.F[i] = f.F[i] * a; o2.F[i] = f.F[i] * (1 - a); });
      ctx.minerals.forEach((m, i) => { const a = tab[m] ?? d; o1.S[i] = f.S[i].map((v) => v * a); o2.S[i] = f.S[i].map((v) => v * (1 - a)); });
      o1.flashTP(); o2.flashTP();
      this.results = { "Duty to keep T [kW]": (o1.H + o2.H - f.H) / 1e3 }; return { out1: o1, out2: o2 };
    }
  }
  reg(ComponentSplitter, "ComponentSplitter", "Fluid", "Component separator",
    "Black-box separator: specify the fraction of each component (or mineral) sent to out1.", ["in"], ["out1", "out2"], [
      P_("splits", "text", "water: 0.99", "", "Lines 'name: fraction to out1'; unlisted species go to out2"),
      P_("default", "number", 0.0, "-", "Fraction to out1 for unlisted species")]);

  function speciesIndex(ctx, tok) {
    tok = tok.trim().toLowerCase();
    const i = ctx.comps.findIndex((c) => c.name.toLowerCase() === tok || c.formula.toLowerCase() === tok);
    if (i < 0) throw new Error(`Species '${tok}' is not in the component list`); return i;
  }
  function parseReactions(text, ctx) {
    const out = [];
    for (let line of (text || "").split("\n")) {
      line = line.trim(); if (!line || line.startsWith("#")) continue;
      const parts = line.split(":").map((s) => s.trim());
      if (parts.length !== 3) throw new Error(`Reaction line needs 'equation : key : conversion' -> ${line}`);
      const [eq, key, conv] = parts;
      if (!eq.includes("->")) throw new Error(`Reaction needs '->': ${eq}`);
      const [lhs, rhs] = eq.split("->"), nu = zeros(ctx.names.length);
      for (const [side, sign] of [[lhs, -1], [rhs, 1]]) for (let term of side.split("+")) {
        term = term.trim(); const m = term.match(/^([0-9]*\.?[0-9]+)?\s*(.+)$/);
        nu[speciesIndex(ctx, m[2])] += sign * (m[1] ? parseFloat(m[1]) : 1);
      }
      const k = speciesIndex(ctx, key);
      if (nu[k] >= 0) throw new Error(`Key component ${key} must be a reactant`);
      out.push({ eq: eq.trim(), nu, key: k, conv: parseFloat(conv) });
    }
    return out;
  }
  class ConversionReactor extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx); let F = f.F.slice(); const ext = {};
      for (const r of parseReactions(this.p.reactions, ctx)) {
        const xi = r.conv * F[r.key] / -r.nu[r.key];
        F = F.map((v, i) => v + xi * r.nu[i]);
        if (F.some((v) => v < -1e-9)) throw new Error(`Reaction '${r.eq}' consumes more of a reactant than is available`);
        F = F.map((v) => Math.max(v, 0)); ext[r.eq] = xi;
      }
      const out = new Stream(ctx, f.T, f.P - (this.num("dP") || 0) * 1e5, F, f.S, `${this.id}.out`);
      if (this.p.mode === "isothermal") out.flashTP(f.T);
      else if (this.p.mode === "outlet temperature") out.flashTP(this.num("T_out") + K0);
      else out.flashPH(f.H);
      this.results = { "Duty [kW]": (out.H - f.H) / 1e3, "Outlet T [C]": out.T - K0 };
      for (const [k, v] of Object.entries(ext)) this.results[`Extent: ${k} [mol/s]`] = v;
      return { out };
    }
  }
  reg(ConversionReactor, "ConversionReactor", "Fluid", "Conversion reactor",
    "Stoichiometric reactor. Reactions are applied in sequence, each to a fractional conversion of its key reactant. Heat of reaction comes from formation enthalpies.",
    ["in"], ["out"], [
      P_("reactions", "text", "CO + 2 H2 -> CH3OH : CO : 0.3", "", "One per line: 'equation : key reactant : conversion'"),
      P_("mode", "select", "isothermal", "", "Energy specification", ["isothermal", "outlet temperature", "adiabatic"]),
      P_("T_out", "number", 250.0, "C", "Outlet temperature", null, { mode: "outlet temperature" }),
      P_("dP", "number", 0.0, "bar", "Pressure drop")]);

  class ShortcutColumn extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), pkg = ctx.pkg;
      if (!pkg) throw new Error("Column needs volatile components");
      const P = this.num("P") * 1e5, vol = ctx.vol, names = vol.map((i) => ctx.names[i]);
      const lk = names.indexOf(this.p.LK), hk = names.indexOf(this.p.HK);
      if (lk < 0 || hk < 0) throw new Error("Choose a light key and a heavy key from the component list");
      const Fv = vol.map((i) => f.F[i]), Ft = sum(Fv), z = Fv.map((v) => v / Ft);
      const rl = this.num("LK_rec"), rh = this.num("HK_rec");
      const Tbf = pkg.bubbleT(z, P), Tdf = pkg.dewT(z, P);
      const hL = pkg.h(Tbf, P, z, "L"), hV = pkg.h(Tdf, P, z, "V"), hF = (f.H - f._hNonvol(f.T)) / Ft;
      const q = (hV - hF) / (hV - hL);
      const Keq = (T, c) => { const fr = pkg.flashTP(T, P, c); return fr.phase === "VL" ? fr.y.map((y, i) => y / Math.max(fr.x[i], 1e-300)) : pkg.wilsonK(T, P); };
      let alpha = Keq(0.5 * (Tbf + Tdf), z); alpha = alpha.map((a) => a / alpha[hk]);
      let d, b, xD, xB, Ttop, Tbot, Nmin;
      for (let it = 0; it < 6; it++) {
        Nmin = Math.log((rl / (1 - rl)) * (rh / (1 - rh))) / Math.log(alpha[lk]);
        const dhb = (1 - rh) / rh, ratio = alpha.map((a) => dhb * a ** Nmin);
        ratio[lk] = rl / (1 - rl); ratio[hk] = dhb;
        d = Fv.map((v, i) => v * ratio[i] / (1 + ratio[i])); b = Fv.map((v, i) => v - d[i]);
        xD = normalize(d); xB = normalize(b);
        Ttop = pkg.dewT(xD, P); Tbot = pkg.bubbleT(xB, P);
        const Kt = Keq(Ttop - 0.01, xD), Kb = Keq(Tbot + 0.01, xB);
        const an = Kt.map((k, i) => Math.sqrt((k / Kt[hk]) * (Kb[i] / Kb[hk])));
        const conv = an.every((a, i) => Math.abs(a / alpha[i] - 1) < 1e-6); alpha = an; if (conv) break;
      }
      if (alpha[lk] <= 1) throw new Error("Light key is not more volatile than heavy key at column conditions");
      Nmin = Math.log((rl / (1 - rl)) * (rh / (1 - rh))) / Math.log(alpha[lk]);
      const theta = brent((th) => sum(alpha.map((a, i) => a * z[i] / (a - th))) - (1 - q), 1 + 1e-9, alpha[lk] - 1e-9);
      const Rmin = Math.max(sum(alpha.map((a, i) => a * xD[i] / (a - theta))) - 1, 1e-6), Rr = this.num("R_factor") * Rmin;
      const X = (Rr - Rmin) / (Rr + 1);
      const Y = X > 0 ? 1 - Math.exp((1 + 54.4 * X) / (11 + 117.2 * X) * (X - 1) / Math.sqrt(X)) : 0;
      const N = (Y + Nmin) / (1 - Y), D = sum(d), B = sum(b);
      const kirk = ((z[hk] / z[lk]) * (xB[lk] / xD[hk]) ** 2 * (B / D)) ** 0.206;
      const Tdb = pkg.bubbleT(xD, P);
      const dist = new Stream(ctx, Tdb, P, null, null, `${this.id}.distillate`); vol.forEach((i, k) => (dist.F[i] = d[k])); dist.flashTP(Tdb, P);
      const bott = new Stream(ctx, Tbot, P, f.F, f.S, `${this.id}.bottoms`); vol.forEach((i, k) => (bott.F[i] = b[k])); bott.flashTP(Tbot, P);
      const Qc = (Rr + 1) * D * (pkg.h(Ttop, P, xD, "V") - pkg.h(Tdb, P, xD, "L"));
      const Qr = dist.H + bott.H + Qc - f.H;
      this.results = {
        "Minimum stages (Fenske)": Nmin, "Minimum reflux (Underwood)": Rmin, "Reflux ratio": Rr, "Theoretical stages": N,
        "Feed stage from top (Kirkbride)": N / (1 + kirk), "Feed q": q, "Condenser duty [kW]": -Qc / 1e3, "Reboiler duty [kW]": Qr / 1e3,
        "Top T [C]": Ttop - K0, "Bottom T [C]": Tbot - K0, "LK/HK relative volatility": alpha[lk],
      };
      return { distillate: dist, bottoms: bott };
    }
  }
  reg(ShortcutColumn, "ShortcutColumn", "Fluid", "Distillation (shortcut)",
    "Fenske-Underwood-Gilliland column with Kirkbride feed location and a total condenser. Non-key components are distributed by Fenske; volatilities are re-evaluated at top and bottom conditions.",
    ["in"], ["distillate", "bottoms"], [
      P_("LK", "component", "", "", "Light key"), P_("HK", "component", "", "", "Heavy key"),
      P_("LK_rec", "number", 0.99, "-", "Light-key recovery in distillate"), P_("HK_rec", "number", 0.99, "-", "Heavy-key recovery in bottoms"),
      P_("R_factor", "number", 1.3, "-", "Reflux ratio / minimum reflux"), P_("P", "number", 1.01325, "bar", "Column pressure")]);

  // ---- minerals
  const splitS = (S, frac) => { const a = S.map((r, i) => r.map((v, j) => v * frac[i][j])); return [a, S.map((r, i) => r.map((v, j) => v - a[i][j]))]; };

  class Comminution extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), out = f.copy(`${this.id}.out`), tot = colSum(f.S);
      if (sum(tot) <= 0) throw new Error("Comminution feed has no solids");
      const F80 = pSize(tot), target = rrsbFractions(this.num("P80"), this.num("spread")), pasP = passingAtSieves(target);
      f.S.forEach((m, i) => {
        const ms = sum(m); if (ms <= 0) return;
        const pasF = passingAtSieves(m), pas = pasF.map((v, j) => Math.max(v, pasP[j]));
        out.S[i] = fracFromPassing(pas).map((v) => v * ms);
      });
      const P80 = pSize(colSum(out.S));
      let wi = this.num("work_index");
      if (wi === null) { const ms = f.S.map(sum); wi = sum(ms.map((m, i) => m * ctx.mins[i].workIndex)) / sum(ms); }
      const E = F80 > P80 ? 10 * wi * (1 / Math.sqrt(P80) - 1 / Math.sqrt(F80)) : 0, tph = sum(tot) * 3.6;
      out.flashTP();
      this.results = { "Feed F80 [um]": F80, "Product P80 [um]": P80, "Work index [kWh/t]": wi, "Specific energy [kWh/t]": E, "Solids [t/h]": tph, "Power [kW]": E * tph };
      return { out };
    }
  }
  reg(Comminution, "Comminution", "Minerals", "Crusher / mill",
    "Size reduction to a target P80. Product size follows a Rosin-Rammler curve (never coarser than the feed); specific energy from Bond's law.", ["in"], ["out"], [
      P_("P80", "number", 150.0, "um", "Target product P80"), P_("spread", "number", 1.0, "-", "Rosin-Rammler spread n of the product"),
      P_("work_index", "number", "", "kWh/t", "Bond work index (blank = mass-weighted mineral values)")]);

  class Classifier extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), a = this.num("sharpness"), d50 = this.num("d50c"), Rf = clamp(this.num("Rf"), 0, 1);
      const Ec = MID.map((d) => {
        const ax = Math.min(a * d / d50, 600); let E = (Math.exp(ax) - 1) / (Math.exp(ax) + Math.exp(Math.min(a, 600)) - 2);
        if (!isFinite(E)) E = 1; E = clamp(E, 0, 1); return Rf + (1 - Rf) * E;
      });
      const [sc, sf] = splitS(f.S, f.S.map(() => Ec));
      const c = new Stream(ctx, f.T, f.P, f.F.map((v) => v * Rf), sc, `${this.id}.coarse`).flashTP();
      const fi = new Stream(ctx, f.T, f.P, f.F.map((v) => v * (1 - Rf)), sf, `${this.id}.fine`).flashTP();
      this.results = { "Coarse solids [t/h]": c.solidsMass * 3.6, "Fine solids [t/h]": fi.solidsMass * 3.6, "Fine P80 [um]": fi.p80(), "Coarse P80 [um]": c.p80() };
      return { coarse: c, fine: fi };
    }
  }
  reg(Classifier, "Classifier", "Minerals", "Hydrocyclone / screen",
    "Size classification with a Whiten partition curve: coarse recovery = Rf + (1-Rf)E(d/d50c). Use a sharpness of 2-4 for cyclones, 8-15 for screens.",
    ["in"], ["coarse", "fine"], [
      P_("d50c", "number", 120.0, "um", "Corrected cut size"), P_("sharpness", "number", 3.0, "-", "Whiten sharpness alpha"),
      P_("Rf", "number", 0.3, "-", "Fraction of fluid (and fines bypass) to coarse product")]);

  const DEFAULT_KINETICS = "chalcopyrite: k=1.2, Rmax=0.95\nsphalerite: k=0.9, Rmax=0.92\ngalena: k=1.4, Rmax=0.95\npentlandite: k=0.8, Rmax=0.90\npyrite: k=0.25, Rmax=0.60\nquartz: k=0.02, Rmax=0.10";
  class FlotationBank extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), tab = parseTable(this.p.kinetics, ["k", "Rmax"]);
      const tau = this.num("tau"), N = Math.max(Math.trunc(this.num("cells")), 1), Rw = clamp(this.num("water_rec"), 0, 1), ent = this.num("entrainment");
      const sf = MID.map((d) => (this.p.size_effect === "on" ? Math.exp(-0.5 * (Math.log(d / 50) / 1.7) ** 2) : 1));
      const entf = MID.map((d) => (d < 40 ? ent * Rw * (1 - d / 40) : 0));
      const frac = ctx.minerals.map((m) => {
        const kv = tab[m] || { k: 0, Rmax: 0 }, k = kv.k ?? 0, Rmax = kv.Rmax ?? 1;
        return MID.map((_, j) => { const Rk = Rmax * (1 - (1 + k * tau / N) ** -N) * sf[j]; return clamp(Rk + (1 - Rk) * entf[j], 0, 1); });
      });
      const [sc, st] = splitS(f.S, frac);
      const conc = new Stream(ctx, f.T, f.P, f.F.map((v) => v * Rw), sc, `${this.id}.concentrate`).flashTP();
      const tail = new Stream(ctx, f.T, f.P, f.F.map((v) => v * (1 - Rw)), st, `${this.id}.tailings`).flashTP();
      const cm = conc.solidsMass;
      this.results = { "Concentrate solids [t/h]": cm * 3.6, "Mass pull [%]": 100 * cm / Math.max(f.solidsMass, 1e-30) };
      ctx.minerals.forEach((m, i) => {
        const fin = sum(f.S[i]); if (fin <= 0) return;
        this.results[`Recovery ${m} [%]`] = 100 * sum(sc[i]) / fin;
        this.results[`Grade ${m} in conc [% w]`] = cm > 0 ? 100 * sum(sc[i]) / cm : 0;
      });
      return { concentrate: conc, tailings: tail };
    }
  }
  reg(FlotationBank, "FlotationBank", "Minerals", "Flotation bank",
    "Bank of N perfectly mixed cells with first-order kinetics per mineral: R = Rmax[1-(1+k tau/N)^-N] times a size factor, plus entrainment of fines with water.",
    ["in"], ["concentrate", "tailings"], [
      P_("kinetics", "text", DEFAULT_KINETICS, "", "Lines 'mineral: k=rate [1/min], Rmax=fraction'"),
      P_("tau", "number", 12.0, "min", "Total residence time"), P_("cells", "number", 6, "-", "Number of cells in series"),
      P_("water_rec", "number", 0.15, "-", "Fraction of water (fluid) reporting to concentrate"),
      P_("entrainment", "number", 0.6, "-", "Entrainment factor for particles below 40 um"),
      P_("size_effect", "select", "on", "", "Reduce floatability of very fine and coarse particles", ["on", "off"])]);

  class LeachTank extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), tab = parseTable(this.p.extraction, ["x"]), out = f.copy(`${this.id}.out`), diss = {};
      ctx.minerals.forEach((m, i) => {
        const v = tab[m]; if (!v || !Object.keys(v).length) return;
        const x = clamp(first(v), 0, 1), mass = sum(f.S[i]) * x; if (mass <= 0) return;
        out.S[i] = f.S[i].map((s) => s * (1 - x));
        const mol = mass / (ctx.mins[i].mw / 1000); let metal = 0;
        for (const [sp, nu] of Object.entries(ctx.mins[i].leach)) { const j = ctx.index(sp); out.F[j] += mol * nu; metal += mol * nu * ctx.mw[j]; }
        const j = ctx.index("dissolved other"); out.F[j] += (mass - metal) / ctx.mw[j];
        diss[m] = 100 * x;
      });
      out.flashTP(this.num("T") + K0);
      this.results = { "Duty [kW]": (out.H - f.H) / 1e3, "Residence time [min]": this.num("tau") };
      for (const [m, x] of Object.entries(diss)) this.results[`Extraction ${m} [%]`] = x;
      for (const sp of ["Cu2+", "Zn2+", "Ni2+", "Li+", "Pb2+"]) if (ctx.names.includes(sp)) {
        const j = ctx.index(sp), dm = (out.F[j] - f.F[j]) * ctx.mw[j]; if (dm > 0) this.results[`${sp} dissolved [t/h]`] = dm * 3.6;
      }
      return { out };
    }
  }
  reg(LeachTank, "LeachTank", "Minerals", "Leach reactor",
    "Dissolves a specified fraction of each mineral into the liquor as metal ions. Non-metal mass is carried as 'dissolved other' so the mass balance closes. Operates at a set temperature; the duty is reported.",
    ["in"], ["out"], [
      P_("extraction", "text", "chalcopyrite: 0.90\npyrite: 0.10", "", "Lines 'mineral: fraction dissolved'"),
      P_("T", "number", 80.0, "C", "Operating temperature"), P_("tau", "number", 240.0, "min", "Residence time (reported only)")]);

  class Thickener extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), r = clamp(this.num("solids_rec"), 0, 1), c = clamp(this.num("uf_solids") / 100, 1e-3, 0.999);
      const Su = f.S.map((row) => row.map((v) => v * r)), need = sum(Su.map(sum)) * (1 - c) / c, fl = f.fluidMass;
      const a = fl > 0 ? Math.min(need / fl, 1) : 0;
      const u = new Stream(ctx, f.T, f.P, f.F.map((v) => v * a), Su, `${this.id}.underflow`).flashTP();
      const o = new Stream(ctx, f.T, f.P, f.F.map((v) => v * (1 - a)), f.S.map((row, i) => row.map((v, j) => v - Su[i][j])), `${this.id}.overflow`).flashTP();
      this.results = { "Underflow solids [% w]": ctx.minerals.length ? u.summary().solids_pct_w : 0, "Liquor to overflow [%]": 100 * (1 - a) };
      return { underflow: u, overflow: o };
    }
  }
  reg(Thickener, "Thickener", "Minerals", "Thickener / filter", "Solid-liquid separation to a target underflow solids content.", ["in"], ["underflow", "overflow"], [
    P_("solids_rec", "number", 0.995, "-", "Fraction of solids to underflow"), P_("uf_solids", "number", 60.0, "% w", "Underflow solids content")]);

  // ---- hydrometallurgy and refining
  function tridiagSolve(lo, di, up, rhs) {
    const n = di.length, c = zeros(n), d = zeros(n);
    c[0] = up[0] / di[0]; d[0] = rhs[0] / di[0];
    for (let i = 1; i < n; i++) {
      const m = di[i] - lo[i] * c[i - 1];
      c[i] = i < n - 1 ? up[i] / m : 0;
      d[i] = (rhs[i] - lo[i] * d[i - 1]) / m;
    }
    const x = zeros(n); x[n - 1] = d[n - 1];
    for (let i = n - 2; i >= 0; i--) x[i] = d[i] - c[i] * x[i + 1];
    return x;
  }
  const DEFAULT_BETA = "La3+: 1\nCe3+: 6\nPr3+: 12\nNd3+: 17.4\nSm+HREE3+: 150\nTh4+: 2000\nCa2+: 0.002";
  class SXCascade extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), beta = {};
      for (const [k, v] of Object.entries(parseTable(this.p.betas, ["beta"]))) beta[k] = first(v, 0);
      const key = this.p.key;
      if (!(key in beta) || beta[key] <= 0) throw new Error(`Key species '${key}' needs a separation factor in the list`);
      const Dk = this.num("D_key"), Ne = Math.max(Math.trunc(this.num("N_ext")), 1), Ns = Math.max(Math.trunc(this.num("N_scrub")), 0);
      const OA = this.num("OA"), SA = this.num("SA"), sw = clamp(this.num("strip_water"), 0, 1), N = Ne + Ns;
      const Lf = 1 + SA, Ls = SA > 0 ? SA : 1e-12;
      const ext = new Stream(ctx, f.T, f.P, null, null, `${this.id}.extract`);
      const raf = new Stream(ctx, f.T, f.P, null, f.S, `${this.id}.raffinate`);
      const rec = {};
      ctx.names.forEach((n, i) => {
        const Fi = f.F[i], c = ctx.comps[i];
        if (c.solute && n in beta && Fi > 0) {
          const D = Dk * beta[n] / beta[key];
          const e = Array.from({ length: N }, (_, k) => D * OA / (k < Ne ? Lf : Ls));
          const lo = [0].concat(e.slice(0, N - 1)), di = e.map((v) => -(1 + v));
          const up = Array(N - 1).fill(1).concat([0]), rhs = zeros(N); rhs[Ne - 1] = -Fi;
          const a = tridiagSolve(lo, di, up, rhs);
          const toExt = Math.min(Math.max(e[N - 1] * a[N - 1], 0), Fi);
          ext.F[i] = toExt; raf.F[i] = Fi - toExt; rec[n] = toExt / Fi;
        } else if (c.solute || !ctx.vol.includes(i)) raf.F[i] = Fi;
        else { ext.F[i] = Fi * sw; raf.F[i] = Fi * (1 - sw); }
      });
      ext.flashTP(); raf.flashTP();
      this.results = { "Stages (extraction + scrub)": N };
      const reIons = Object.keys(rec).filter((n) => n !== "Ca2+" && n !== "Th4+");
      const tot = sum(reIons.map((n) => ext.F[ctx.index(n)]));
      for (const [n, r] of Object.entries(rec)) this.results[`${n} to extract [%]`] = 100 * r;
      if (tot > 0) this.results[`${key} share of RE in extract [%]`] = 100 * ext.F[ctx.index(key)] / tot;
      return { extract: ext, raffinate: raf };
    }
  }
  reg(SXCascade, "SXCascade", "Hydrometallurgy", "Solvent extraction cascade",
    "Counter-current mixer-settler bank with an extraction section and a scrub section, solved stage by stage with constant distribution ratios D = D(key) x beta / beta(key). Typical P507 (PC88A) chloride-system separation factors are the default. Loaded organic is fully stripped; organic, scrub and strip acid circulate inside the unit and are not tracked.",
    ["in"], ["extract", "raffinate"], [
      P_("betas", "text", DEFAULT_BETA, "", "Separation factors 'species: beta' (unlisted species stay in the aqueous)"),
      P_("key", "solute", "Nd3+", "", "Key species for the distribution ratio"),
      P_("D_key", "number", 1.0, "-", "Distribution ratio of the key species (organic/aqueous)"),
      P_("N_ext", "number", 30, "-", "Extraction stages"), P_("N_scrub", "number", 30, "-", "Scrub stages"),
      P_("OA", "number", 1.0, "-", "Organic flow / aqueous feed flow (by volume)"),
      P_("SA", "number", 0.5, "-", "Scrub acid flow / aqueous feed flow (by volume)"),
      P_("strip_water", "number", 0.3, "-", "Strip liquor volume / aqueous feed volume")]);

  function parsePrecip(text, ctx) {
    const rx = [];
    for (let line of (text || "").split("\n")) {
      line = line.trim(); if (!line || line.startsWith("#")) continue;
      const parts = line.split(":").map((p) => p.trim());
      if (parts.length !== 2 || !parts[0].includes("->")) throw new Error(`Precipitation line needs 'ion -> solid : fraction' -> ${line}`);
      const [ion, solid] = parts[0].split("->").map((x) => x.trim());
      rx.push([ctx.index(ion), ctx.mindex(solid), parseFloat(parts[1])]);
    }
    return rx;
  }
  class Precipitator extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), out = f.copy(), jr = ctx.names.indexOf("reagent");
      let used = 0; const frac = rrsbFractions(this.num("p80"), 1.2), made = {};
      for (let [ion, sol, x] of parsePrecip(this.p.reactions, ctx)) {
        x = clamp(x, 0, 1);
        const m = ctx.mins[sol], nu = m.leach[ctx.names[ion]];
        if (!nu) throw new Error(`${m.name} does not contain ${ctx.names[ion]}`);
        const molIon = out.F[ion] * x, molSolid = molIon / nu, need = molSolid * m.mw / 1000 - molIon * ctx.mw[ion];
        out.F[ion] -= molIon;
        out.S[sol] = out.S[sol].map((v, j) => v + molSolid * m.mw / 1000 * frac[j]);
        used += need; made[m.name] = (made[m.name] || 0) + molSolid * m.mw / 1000;
      }
      if (used > 1e-12) {
        if (jr < 0) throw new Error("Add the 'reagent' species and a reagent feed stream to supply the precipitant");
        const have = out.F[jr] * ctx.mw[jr];
        if (used > have * (1 + 1e-9)) throw new Error(`Not enough reagent: needs ${(used * 3.6).toPrecision(4)} t/h, inlet carries ${(have * 3.6).toPrecision(4)} t/h`);
        out.F[jr] = Math.max(out.F[jr] - used / ctx.mw[jr], 0);
      }
      const mc = clamp(this.num("moisture") / 100, 0, 0.99), ms = out.solidsMass, fl = out.fluidMass;
      const a = fl > 0 ? Math.min(ms * mc / (1 - mc) / fl, 1) : 0, T = this.num("T") + K0;
      const cake = new Stream(ctx, T, f.P, out.F.map((v) => v * a), out.S, `${this.id}.solids`).flashTP();
      const liq = new Stream(ctx, T, f.P, out.F.map((v) => v * (1 - a)), null, `${this.id}.liquor`).flashTP();
      this.results = { "Duty [kW]": (cake.H + liq.H - f.H) / 1e3, "Reagent consumed [t/h]": used * 3.6 };
      for (const [k, v] of Object.entries(made)) this.results[`${k} precipitated [t/h]`] = v * 3.6;
      return { solids: cake, liquor: liq };
    }
  }
  reg(Precipitator, "Precipitator", "Hydrometallurgy", "Precipitation and filtration",
    "Precipitates a fraction of each dissolved ion as a named solid, then filters. The anion or oxygen mass of the solid is drawn from the 'reagent' species, so feed a reagent stream to the inlet. Oxalate precipitation followed by calcination can be represented directly as the oxide.",
    ["in"], ["solids", "liquor"], [
      P_("reactions", "text", "Nd3+ -> Nd2O3 : 0.995", "", "Lines 'ion -> solid : fraction precipitated'"),
      P_("T", "number", 60.0, "C", "Operating temperature"),
      P_("moisture", "number", 15.0, "% w", "Liquor content of the filter cake"),
      P_("p80", "number", 20.0, "um", "Precipitate particle size P80")]);

  class Electrolysis extends Unit {
    run(inl, ctx) {
      const f = this.feed(inl, ctx), out = f.copy();
      if (!ctx.names.includes("oxygen")) throw new Error("Add 'oxygen' to the components to carry the anode gas");
      const jo = ctx.index("oxygen"); let metalKg = 0;
      for (let line of (this.p.reactions || "").split("\n")) {
        line = line.trim(); if (!line || line.startsWith("#")) continue;
        let [eq, x] = line.split(":").map((p) => p.trim());
        const [ox, me] = eq.split("->").map((t) => ctx.mindex(t.trim()));
        x = clamp(parseFloat(x), 0, 1);
        const mo = ctx.mins[ox], mm = ctx.mins[me], ion = Object.keys(mm.leach)[0];
        const nOx = mo.leach[ion], nMe = mm.leach[ion];
        if (!nOx || !nMe) throw new Error(`${mo.name} and ${mm.name} must share a metal`);
        const massOx = out.S[ox].map((v) => v * x), tot = sum(massOx);
        const molMetal = tot / (mo.mw / 1000) * nOx, mMetal = molMetal / nMe * mm.mw / 1000;
        out.S[ox] = out.S[ox].map((v, j) => v - massOx[j]);
        if (tot > 0) out.S[me] = out.S[me].map((v, j) => v + massOx[j] / tot * mMetal);
        out.F[jo] += (tot - mMetal) / ctx.mw[jo];
        metalKg += mMetal;
      }
      out.flashTP(this.num("T") + K0);
      const [V, L] = out.splitPhases(); L.name = `${this.id}.metal`; V.name = `${this.id}.offgas`;
      this.results = { "Metal produced [t/h]": metalKg * 3.6, "Electrical power [kW]": this.num("energy") * metalKg * 3600,
        "Current efficiency": this.num("current_eff"), "Thermal duty [kW]": (out.H - f.H) / 1e3 };
      return { metal: L, offgas: V };
    }
  }
  reg(Electrolysis, "Electrolysis", "Hydrometallurgy", "Molten-salt electrolysis",
    "Reduces oxide feed to metal in a fluoride melt. Oxygen leaves as anode gas (reported as O2; with carbon anodes it is really CO and CO2). Electrical energy from a specific consumption.",
    ["in"], ["metal", "offgas"], [
      P_("reactions", "text", "Nd2O3 -> Nd metal : 0.98", "", "Lines 'oxide -> metal : fraction reduced'"),
      P_("T", "number", 1050.0, "C", "Cell temperature"),
      P_("energy", "number", 10.0, "kWh/kg", "Electrical energy per kg metal"),
      P_("current_eff", "number", 0.75, "-", "Current efficiency (reported)")]);

  // ------------------------------------------------------------------ flowsheet
  function buildContext(spec) {
    const comps = (spec.components || []).slice(), mins = (spec.minerals || []).slice();
    if (mins.length) {
      for (const m of mins) for (const sp of Object.keys(mineral(m).leach)) if (!comps.includes(sp)) comps.push(sp);
      if (!comps.includes("dissolved other")) comps.push("dissolved other");
    }
    return new Context(comps, mins, spec.package || "Peng-Robinson");
  }
  function tarjan(nodes, edges) {
    const index = {}, low = {}, stack = [], on = new Set(), out = []; let counter = 0;
    function strong(v) {
      index[v] = low[v] = counter++; stack.push(v); on.add(v);
      for (const w of edges[v] || []) {
        if (!(w in index)) { strong(w); low[v] = Math.min(low[v], low[w]); }
        else if (on.has(w)) low[v] = Math.min(low[v], index[w]);
      }
      if (low[v] === index[v]) { const c = []; let w; do { w = stack.pop(); on.delete(w); c.push(w); } while (w !== v); out.push(c); }
    }
    for (const v of nodes) if (!(v in index)) strong(v);
    return out.reverse();
  }

  class Flowsheet {
    constructor(spec, { tol = 1e-6, maxIter = 100 } = {}) {
      this.spec = spec; this.tol = tol; this.maxIter = maxIter; this.ctx = buildContext(spec);
      this.units = {}; this.unitOrder = [];
      for (const u of spec.units || []) {
        if (!UNIT_TYPES[u.type]) throw new Error(`Unknown unit type '${u.type}' for ${u.id}`);
        this.units[u.id] = new UNIT_TYPES[u.type](u.id, u.params); this.unitOrder.push(u.id);
      }
      this.streams = {}; this.streamOrder = [];
      for (const s of spec.streams || []) { this.streams[s.id] = s; this.streamOrder.push(s.id); }
      this.validate(); this.values = {}; this.log = [];
    }
    validate() {
      for (const sid of this.streamOrder) {
        const s = this.streams[sid];
        for (const end of ["from", "to"]) {
          const e = s[end]; if (!e) continue;
          const u = this.units[e.unit]; if (!u) throw new Error(`Stream ${sid}: unknown unit ${e.unit}`);
          const ports = end === "from" ? u.constructor.OUTPORTS : u.constructor.INPORTS;
          if (!ports.includes(e.port)) throw new Error(`Stream ${sid}: unit ${u.id} has no ${end === "from" ? "outlet" : "inlet"} port '${e.port}'`);
        }
        if (!s.from && !s.feed) throw new Error(`Stream ${sid} has no source unit and no feed specification`);
      }
      for (const uid of this.unitOrder)
        if (!this.streamOrder.some((sid) => this.streams[sid].to && this.streams[sid].to.unit === uid)) throw new Error(`Unit ${uid} has no inlet stream`);
    }
    edges() {
      const e = Object.fromEntries(this.unitOrder.map((u) => [u, []]));
      for (const sid of this.streamOrder) { const s = this.streams[sid]; if (s.from && s.to) e[s.from.unit].push(s.to.unit); }
      return e;
    }
    plan() {
      const edges = this.edges(), sccs = tarjan(this.unitOrder, edges), plan = [];
      for (const comp of sccs) {
        const cset = new Set(comp);
        const internal = this.streamOrder.filter((sid) => { const s = this.streams[sid]; return s.from && s.to && cset.has(s.from.unit) && cset.has(s.to.unit); });
        if (comp.length === 1 && !internal.length) { plan.push(["unit", comp[0]]); continue; }
        const hasExt = (u) => this.streamOrder.some((sid) => { const s = this.streams[sid]; return s.to && s.to.unit === u && (!s.from || !cset.has(s.from.unit)); });
        const starts = comp.slice().sort((a, b) => (hasExt(a) === hasExt(b) ? (a < b ? -1 : a > b ? 1 : 0) : hasExt(a) ? -1 : 1));
        const color = {}, tears = [], byEdge = {};
        for (const sid of internal) { const s = this.streams[sid], k = s.from.unit + "\u0000" + s.to.unit; (byEdge[k] = byEdge[k] || []).push(sid); }
        const dfs = (v) => {
          color[v] = 1;
          for (const w of edges[v]) {
            if (!cset.has(w)) continue;
            if (color[w] === 1) tears.push(...byEdge[v + "\u0000" + w]);
            else if (!(w in color)) dfs(w);
          }
          color[v] = 2;
        };
        for (const st of starts) if (!(st in color)) dfs(st);
        const tearSet = new Set(tears), indeg = Object.fromEntries(comp.map((u) => [u, 0]));
        for (const sid of internal) if (!tearSet.has(sid)) indeg[this.streams[sid].to.unit]++;
        const ready = starts.filter((u) => indeg[u] === 0), seq = [];
        while (ready.length) {
          const u = ready.shift(); seq.push(u);
          for (const sid of internal) {
            const s = this.streams[sid];
            if (!tearSet.has(sid) && s.from.unit === u) { if (--indeg[s.to.unit] === 0) ready.push(s.to.unit); }
          }
        }
        plan.push(["loop", seq, [...tearSet].sort()]);
      }
      return plan;
    }
    inlets(uid) {
      const d = {};
      for (const sid of this.streamOrder) {
        const s = this.streams[sid];
        if (s.to && s.to.unit === uid && this.values[sid]) (d[s.to.port] = d[s.to.port] || []).push(this.values[sid]);
      }
      return d;
    }
    runUnit(uid) {
      const u = this.units[uid], inl = this.inlets(uid); let outs;
      try { outs = u.run(inl, this.ctx); } catch (e) { const err = new Error(`${uid} (${u.constructor.TYPE}): ${e.message}`); err.unit = uid; throw err; }
      let hin = 0; for (const ss of Object.values(inl)) for (const s of ss) hin += s.H;
      const hout = sum(Object.values(outs).map((s) => s.H));
      u.results["Net heat input [kW]"] = (hout - hin) / 1e3; u.netHeat = hout - hin;
      for (const sid of this.streamOrder) { const s = this.streams[sid]; if (s.from && s.from.unit === uid) this.values[sid] = outs[s.from.port].copy(sid); }
    }
    vec(sids) {
      const v = [];
      for (const sid of sids) { const s = this.values[sid]; v.push(...s.F); for (const r of s.S) v.push(...r); v.push(s.T / 100, s.P / 1e6); }
      return v;
    }
    setVec(sids, v) {
      const nF = this.ctx.names.length; let k = 0;
      for (const sid of sids) {
        const F = v.slice(k, k + nF).map((x) => Math.max(x, 0)); k += nF;
        const S = this.ctx.minerals.map(() => { const r = v.slice(k, k + NBINS).map((x) => Math.max(x, 0)); k += NBINS; return r; });
        const T = v[k] * 100, P = v[k + 1] * 1e6; k += 2;
        this.values[sid] = new Stream(this.ctx, T, P, F, S, sid).flashTP();
      }
    }
    solve() {
      const t0 = Date.now(); this.values = {};
      for (const sid of this.streamOrder) { const s = this.streams[sid]; if (!s.from) this.values[sid] = Stream.fromSpec(this.ctx, s.feed, sid); }
      this._plan = this.plan(); this.loops = [];
      for (const block of this._plan) {
        if (block[0] === "unit") { this.runUnit(block[1]); continue; }
        const [, seq, tears] = block;
        for (const sid of tears) if (!this.values[sid]) {
          const g = this.streams[sid].guess;
          this.values[sid] = g ? Stream.fromSpec(this.ctx, g, sid) : new Stream(this.ctx, 298.15, 101325, null, null, sid).flashTP();
        }
        let x = this.vec(tears), xp = null, gp = null, converged = false; const hist = [];
        for (let it = 1; it <= this.maxIter; it++) {
          for (const uid of seq) this.runUnit(uid);
          const g = this.vec(tears), gmax = Math.max(...g.map(Math.abs));
          let err = 0;
          for (let i = 0; i < g.length; i++) { const sc = Math.max(Math.abs(g[i]), 1e-6 * Math.max(gmax, 1e-12) + 1e-12); err = Math.max(err, Math.abs(g[i] - x[i]) / sc); }
          hist.push(err);
          if (err < this.tol) { converged = true; break; }
          let xn;
          if (gp && it > 2) {
            xn = x.map((xi, i) => {
              const dx = xi - xp[i], dg = g[i] - gp[i], sl = Math.abs(dx) > 1e-14 ? dg / dx : 0;
              let q = Math.abs(sl - 1) > 1e-12 ? sl / (sl - 1) : 0;
              if (!isFinite(q)) q = 0; q = clamp(q, -5, 0);
              return q * xi + (1 - q) * g[i];
            });
          } else xn = g;
          xp = x; gp = g; x = xn; this.setVec(tears, x);
        }
        this.loops.push({ units: seq, tears, iterations: hist.length, converged, errors: hist });
        if (!converged) this.log.push(`Recycle loop ${seq.join(", ")} did not converge in ${this.maxIter} iterations (last error ${hist[hist.length - 1].toExponential(2)})`);
      }
      this.elapsed = (Date.now() - t0) / 1000;
      return this;
    }
    balances() {
      const feeds = this.streamOrder.filter((s) => !this.streams[s].from), prods = this.streamOrder.filter((s) => !this.streams[s].to);
      const m = (l) => sum(l.map((s) => this.values[s].mass)), h = (l) => sum(l.map((s) => this.values[s].H));
      const mi = m(feeds), mo = m(prods), hi = h(feeds), ho = h(prods), q = sum(Object.values(this.units).map((u) => u.netHeat || 0));
      return { mass_in_kg_s: mi, mass_out_kg_s: mo, mass_closure_pct: mi ? 100 * (mo - mi) / mi : 0,
        enthalpy_in_kW: hi / 1e3, enthalpy_out_kW: ho / 1e3, heat_added_kW: q / 1e3, energy_closure_kW: (ho - hi - q) / 1e3 };
    }
    results() {
      const streams = {};
      for (const sid of this.streamOrder) if (this.values[sid]) {
        const s = this.streams[sid];
        streams[sid] = Object.assign(this.values[sid].summary(), { name: sid, source: s.from ? s.from.unit : null, destination: s.to ? s.to.unit : null });
      }
      const units = {};
      for (const uid of this.unitOrder) units[uid] = { type: this.units[uid].constructor.TYPE, results: this.units[uid].results };
      return {
        ok: true, components: this.ctx.names, minerals: this.ctx.minerals, package: this.ctx.pkg ? this.ctx.pkg.name : null,
        streams, units, loops: this.loops, order: this._plan.map((b) => (b[0] === "unit" ? b[1] : { loop: b[1], tears: b[2] })),
        balances: this.balances(), warnings: this.log, elapsed_s: this.elapsed,
      };
    }
  }

  function solve(spec, opts) {
    try { return new Flowsheet(spec, opts).solve().results(); }
    catch (e) { return { ok: false, error: e.message, unit: e.unit || null }; }
  }
  function catalogue() {
    return {
      units: Object.values(UNIT_TYPES).map((c) => c.schema()),
      components: COMPONENT_LIST.filter((c) => !c.solute).map((c) => ({ name: c.name, formula: c.formula, mw: c.mw })),
      solutes: COMPONENT_LIST.filter((c) => c.solute).map((c) => c.name),
      minerals: MINERAL_LIST.map((m) => ({ name: m.name, formula: m.formula, density: m.density, leach: m.leach })),
      packages: ["Peng-Robinson", "Ideal"],
    };
  }

  const api = { solve, catalogue, Flowsheet, Stream, Context, PengRobinson, Ideal, UNIT_TYPES, COMPONENTS, MINERALS, SIEVES, MID, rrsbFractions, pSize, passingAtSieves };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.ProcSim = api;
})(typeof window !== "undefined" ? window : globalThis);
