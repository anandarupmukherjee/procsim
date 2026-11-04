# Changelog

## 0.3.0 (2026-11-04)

- Packaged as an installable tool with a `procsim` command: `run`, `examples`, `units`, `editor`, `routes`.
- Example flowsheets and the browser editor ship inside the package (`procsim/data`).
- Neodymium route builder moved into the package as `procsim.ree`, with `report()` and `compare()`.
- CI on GitHub Actions: pytest on Python 3.10 to 3.13, JavaScript parity, editor bundle check.

## 0.2.0 (2026-10-09)

- Hydrometallurgy units: solvent-extraction cascade, precipitation and filtration, molten-salt electrolysis.
- Rare-earth species and minerals. Bastnäsite, monazite and monazite-xenotime splits come from USGS MYB 2020 Table 2 and Iluka (2024).
- Three neodymium routes (bastnäsite ore, Mt Weld-type monazite ore, Eneabba-type concentrate) on one downstream design.

## 0.1.0 (2026-10-09)

- Peng-Robinson and ideal thermodynamics, TP/PH flash, bubble and dew points.
- Fluid and mineral-processing units; sequential-modular solver with Tarjan loop detection and Wegstein acceleration.
- Browser engine (JavaScript port) with a drag-and-drop editor; parity test against the Python core.
