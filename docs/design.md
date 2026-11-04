# procsim design notes

Steady-state flowsheet simulator for chemical, mineral and rare-earth processing. It has a Python engine for research and a browser editor. The editor ships in the package; open it with `procsim editor`.

The neodymium routes and their sources are described in [neodymium-routes.md](neodymium-routes.md).

## Architecture

1. A flowsheet JSON file is the contract. It holds components, minerals, the property package, units, streams and an optional `sources` map that the editor displays.
2. The Python core (`procsim/`, numpy only) solves it.
3. The browser engine (`web/engine.js`) is a line-for-line port. The artifact sandbox blocks the WebAssembly and package fetches that Pyodide needs, so the port replaces it.
4. `tests/parity.js` runs every example through both engines and compares streams and unit results. The worst relative difference is about 4e-11.
5. `web/build.py` bundles the editor, the engine and the examples into one page.
6. `procsim/ree.py` builds the three Nd routes from one parameterised design. `tools/tune_routes.py` calibrates flotation, sizes reagents and chooses SX ratios. `tools/route_report.py` prints the Nd balance and product basket.

## Solver

1. Units become a directed graph. Tarjan's algorithm finds strongly connected components, which are the recycle loops.
2. Inside each loop, a DFS from a unit with an external feed marks the back edges as tear streams.
3. The loop is iterated: two passes of direct substitution, then bounded Wegstein (q in [-5, 0]) on the tear vector. The vector holds flows, solids by size class, T and P.
4. Convergence is reached when the largest relative tear change is below 1e-6.

## Thermodynamics

- Peng-Robinson: vdW mixing, with a kij table for water, CO2, N2 and H2S against hydrocarbons.
- Ideal: Raoult's law, Lee-Kesler vapour pressure, Pitzer-type latent heat.
- Enthalpy is on an ideal-gas formation basis, so gas-phase heats of reaction come out of the balance. Minerals and dissolved ions carry sensible heat only.
- Checked against NIST normal boiling points (within 3 K) and benzene-toluene VLE (within 1 K).

## Unit library

| Group | Units |
|---|---|
| Chemical | Mixer, splitter, heater/cooler, valve, flash drum, component separator, conversion reactor, shortcut distillation (FUG, Kirkbride) |
| Mineral processing | Crusher/mill (Rosin-Rammler product, Bond energy), hydrocyclone/screen (Whiten curve), flotation bank (first-order kinetics, size factor, entrainment), leach reactor, thickener/filter |
| Hydrometallurgy and refining | SX cascade (extraction and scrub sections, constant-D stage-by-stage solve), precipitation and filtration (reagent-limited), molten-salt electrolysis |

Solids are tracked per mineral in 25 size classes (a √2 sieve series from 50 mm). Minerals are fully liberated.

## Species

- 20 fluid components.
- Dissolved ions, including La3+, Ce3+, Pr3+, Nd3+, Sm+HREE3+, Th4+, Ca2+ and Fe2+.
- A generic `reagent` species, which supplies precipitant mass.
- 27 minerals and products:
  - Ores: bastnäsite, monazite and monazite-xenotime, with RE splits from USGS MYB 2020 Table 2 and Iluka 2024.
  - Gangue: calcite, barite and apatite.
  - Precipitated products: RE carbonates, RE oxides, ThO2 and Nd metal.

## Examples

1. NGL recovery and fractionation.
2. Methanol synthesis loop with purge.
3. Copper concentrator with concentrate leach.
4. Neodymium from bastnäsite ore to metal.
5. Neodymium from Mt Weld-type monazite ore to metal.
6. Neodymium from Eneabba-type monazite-xenotime concentrate to metal.

## Checks

- 46 pytest tests, covering:
  - thermo against reference data
  - mass and energy closure on all six examples
  - an analytic single-stage SX case, and the tridiagonal solver against a dense solve
  - precipitator mass conservation and reagent limit
  - Nd product purity and recovery window per route
  - each mineral's RE fractions reproducing the cited oxide split
  - the command line, and the bundled editor page matching its sources
- The JS parity check passes on all six examples.

## Known limits

- Peng-Robinson with kij = 0 is poor for polar mixtures. Activity models are needed for aqueous chemistry.
- The shortcut column has a total condenser only.
- There are no pumps, compressors, rigorous columns, kinetic reactors or liberation classes.
- SX distribution ratios are constant, independent of loading and acidity. As a result, element recoveries do not depend on feed composition, and only purity does.
- Yttrium is lumped with the heavy group at the Sm molar mass.

## Roadmap

1. NRTL or UNIFAC activity models, and aqueous speciation for leach liquors.
2. Loading-dependent SX isotherms, so organic flow and stage count trade off realistically.
3. Split the heavy group (Sm, Gd, Dy, Tb, Y) and add the Dy/Tb SX circuits.
4. An ionic-clay route.
5. Population-balance mill model and liberation classes; fit the flotation size window to test work.
6. Equation-oriented mode for design specs and optimisation.
7. Feed flowsheet recoveries and co-product baskets into cost and price models.
