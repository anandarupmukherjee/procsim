# procsim

A steady-state process simulator for chemical, mineral and rare-earth processing flowsheets. It is in the spirit of Aspen Plus but small enough to read end to end.

One JSON flowsheet format runs in two engines:

- **Python core** (`procsim/`, numpy only): the research engine you script, extend and test.
- **Browser engine** (`web/engine.js`): a line-for-line port behind a drag-and-drop editor that works offline as a single HTML page. `tests/parity.js` checks it against the Python core; they agree to about 1e-11.

## Install

```bash
pip install git+https://github.com/anandarupmukherjee/procsim.git
```

Python 3.9 or later. The only dependency is numpy.

## Command line

```bash
procsim examples                          # list the bundled flowsheets
procsim run methanol_loop                 # solve a bundled example and print the stream table
procsim run my_flowsheet.json -o out.json # solve your own file, write full results as JSON
procsim run neodymium_bastnaesite --json  # machine-readable output
procsim units                             # list unit operations and their models
procsim editor                            # open the drag-and-drop editor in your browser
procsim routes                            # compare the neodymium routes by feedstock
procsim examples --copy ./flowsheets      # copy the example JSON files out to edit
```

`procsim run` exits with code 2 if a recycle loop fails to converge.

## Python

```python
from procsim import Flowsheet
from procsim.ree import ROUTES, build, report

fs = Flowsheet(build(ROUTES["neodymium_monazite_mtweld"])).solve()
print(fs.table())                                   # stream table
print(fs.units["SX-3"].results)                     # unit results
print(report(fs, ROUTES["neodymium_monazite_mtweld"])["nd_balance"])   # Nd balance by stage
```

A flowsheet is plain JSON:

```json
{"components": ["water"], "minerals": ["quartz"], "package": "Ideal",
 "units":   [{"id": "BM", "type": "Comminution", "params": {"P80": 150}}],
 "streams": [{"id": "Ore", "from": null, "to": {"unit": "BM", "port": "in"},
              "feed": {"T": 293.15, "P": 101325, "flows": {"water": 1000},
                       "solids": {"quartz": {"mass": 20, "p80": 9000}}}},
             {"id": "Product", "from": {"unit": "BM", "port": "out"}, "to": null}]}
```

Units in the file are SI (K, Pa, mol/s, kg/s). Unit parameters use engineering units (°C, bar, kW, µm, min). Any flowsheet saved from the editor's Import / export tab runs unchanged in Python.

## What is modelled

| Area | Model |
|---|---|
| Thermodynamics | Peng-Robinson EOS (vdW mixing, typical kij for water, CO2, N2, H2S with hydrocarbons) or Ideal (Raoult, Lee-Kesler vapour pressure, Pitzer latent heat). Enthalpy on a formation basis, so gas-phase heats of reaction are automatic. |
| Flash | TP (successive substitution, negative-flash Rachford-Rice), PH, bubble and dew points |
| Fluid units | Mixer, splitter, heater/cooler (T, duty or vapour fraction), valve, flash drum, component separator, conversion reactor (multiple reactions), shortcut distillation (Fenske-Underwood-Gilliland, Kirkbride) |
| Mineral units | Crusher/mill (Rosin-Rammler product, Bond energy), hydrocyclone/screen (Whiten partition curve), flotation bank (first-order kinetics, size factor, entrainment), leach reactor, thickener/filter |
| Hydrometallurgy | SX cascade (extraction and scrub sections, constant distribution ratios, P507 separation factors), precipitation and filtration (reagent-limited), molten-salt electrolysis |
| Solids | 27 minerals and products tracked by mass in 25 size classes (√2 sieve series, 50 mm down to 17 µm) |
| Solver | Sequential modular. Tarjan strongly connected components find recycle loops, DFS picks tear streams, and Wegstein accelerates convergence. |

## Bundled examples

| Name | What it shows |
|---|---|
| `ngl_fractionation` | Chiller, cold separator, de-ethaniser and de-propaniser (shortcut columns, Peng-Robinson) |
| `methanol_loop` | Two reactions, condensation, purge and recycle |
| `copper_concentrator` | Closed-circuit ball mill, rougher and cleaner flotation, concentrate leach |
| `neodymium_bastnaesite` | Mountain Pass-type ore to Nd metal |
| `neodymium_monazite_mtweld` | Mt Weld-type monazite ore to Nd metal |
| `neodymium_eneabba_concentrate` | Eneabba-type monazite-xenotime concentrate to Nd metal |

The three neodymium routes share one downstream design: acid digestion, Th and Fe removal, mixed RE carbonate, HCl redissolution, three SX splits (Nd/Sm, Ce/Pr, Pr/Nd), oxalate precipitation and calcination, then electrolysis. Feed figures are sourced from public reports; see [docs/neodymium-routes.md](docs/neodymium-routes.md) for the sources, the results and the assumptions still to check.

## Development

```bash
git clone https://github.com/anandarupmukherjee/procsim.git && cd procsim
pip install -e ".[test]"
pytest -q                                              # 46 tests
python tests/make_reference.py && node tests/parity.js # JS engine against the Python core
python examples/build_examples.py                      # regenerate the bundled example JSON
python web/build.py                                    # rebuild the bundled editor page
python tools/tune_routes.py neodymium_monazite_mtweld  # recalibrate a Nd route
```

The editor page in `procsim/data/procsim.html` is generated. After you edit `web/app.html`, `web/engine.js` or the examples, rebuild it with `python web/build.py`; a test fails if it is stale. Any change to a unit model must be made in both `procsim/units.py` and `web/engine.js`, and the parity check confirms they still agree.

Design notes are in [docs/design.md](docs/design.md).

## Known limits

- Peng-Robinson with kij = 0 is poor for polar mixtures. Water-alcohol shows a spurious heat of mixing. Activity models (NRTL, UNIFAC) are the next thermodynamics step.
- The shortcut column has a total condenser only.
- There are no pumps, compressors, rigorous columns, kinetic reactors or liberation classes; minerals are fully liberated.
- SX distribution ratios are constant, independent of loading and acidity. Element recoveries therefore do not depend on the rest of the feed; only purity does.
- Hydrometallurgical heats of reaction are not modelled; duties are sensible heat only. Precipitants are one generic `reagent` species.

## Licence

MIT; see [LICENSE](LICENSE).
