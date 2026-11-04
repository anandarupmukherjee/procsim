# Neodymium routes compared by feedstock

Three feedstocks run through one downstream design:

- **A. Bastnäsite** carbonatite ore (Mountain Pass type), `neodymium_bastnaesite`
- **B. Monazite** in weathered carbonatite ore (Mt Weld type), `neodymium_monazite_mtweld`
- **C. Monazite-xenotime** mineral-sands concentrate (Eneabba type), fed straight to the refinery, `neodymium_eneabba_concentrate`

Run `procsim routes` to reproduce the table below. The configurations are in `procsim/ree.py`. Ore routes run at 100 t/h of ore. The concentrate route runs at 32 kt per year (4.0 t/h at an assumed 8,000 operating hours).

## Sources

| Input | Value | Source |
|---|---|---|
| Mountain Pass RE split | La 34.00, Ce 48.80, Pr 4.20, Nd 11.70 % of TREO | USGS Minerals Yearbook 2020, Rare earths, Table 2 |
| Mountain Pass grades | Feed 7.7 to 8.4 % REO; concentrate 58.0 to 63 % REO | Rozelle, Tarka and Mamula, *Mineral processing of rare earth elements*, OSTI 2020, Table 2 |
| Mt Weld RE split | La 23.88, Ce 47.55, Pr 5.16, Nd 18.13, Sm 2.44, Y 0.76 % of TREO | USGS Minerals Yearbook 2020, Table 2 (Mt Weld Central Lanthanide) |
| Mt Weld ore grade | Ore Reserves 32.0 Mt at 6.4 % TREO | Lynas Rare Earths, Annual Report 2025, Ore Reserves Table 1 |
| Mt Weld RE mineral | Monazite is the major RE mineral in the dolomite carbonatite | Lynas Rare Earths, Annual Report 2024, Directors' Report |
| Eneabba concentrate | 32 kt at ~57 % REO; 90 % simplified weighted average recovery | Iluka Resources, Eneabba refinery presentation 2024, slide 25 |
| Eneabba RE split | La 21.8, Ce 45.0, Pr 4.6, Nd 16.6, Tb 0.2, Dy 0.9 % (slide 25); Sm 3, Gd 1, Y 6 % (slide 27) | Iluka Resources, Eneabba refinery presentation 2024 |

## Results

| | A. Bastnäsite | B. Mt Weld monazite | C. Eneabba concentrate |
|---|---|---|---|
| Feed | 8.04 % REO ore | 6.38 % REO ore | 57 % REO concentrate |
| Concentrate | 61.8 % REO | 40.0 % REO (assumed target) | (feed) |
| Feed to Nd metal | 58.3 % | 49.1 % | 88.7 % |
| Concentrate to Nd oxide | 91.9 % | 91.8 % | 90.5 % |
| Nd oxide purity (TREO basis) | 99.69 % | 99.72 % | 99.63 % |
| LaCe carbonate per t Nd oxide | 9.7 | 5.4 | 5.5 |
| SEG oxide per t Nd oxide | 0.12 | 0.32 | 0.87 |
| ThO2 residue per t Nd oxide | 0.03 | 0.56 | 0.63 |

Findings:

1. The concentrate route is checked against two Iluka figures without tuning. It gives 90.5 % from concentrate to oxide against Iluka's 90 %. Its Nd plus Pr oxide comes to 3.48 kt per year, against Iluka's SX3 capacity check of 3.5 kt NdPr.
2. The bastnäsite feed and concentrate grades sit inside the published Mountain Pass ranges.
3. Mt Weld ore carries more Nd per tonne than Mountain Pass despite its lower grade. Nd is 18.1 % of TREO at Mt Weld and 11.7 % at Mountain Pass. On the assumed flotation parameters, that advantage is lost to fines in the tailings.
4. Monazite feeds produce about half the LaCe per tonne of Nd oxide, and several times the heavy rare earths (SEG).
5. Thorium residue scales with monazite. The size of the gap follows from the assumed Th content, which is not sourced.

## Assumptions still to check

- Thorium: 6 mol % in both monazite feeds; 0.002 mol per mol in bastnäsite.
- Mt Weld concentrate target of 40 % REO, gangue split (hematite, apatite, quartz) and grind size.
- Flotation rate constants and Rmax for all ore routes. Recovery is set mainly by the assumed size and floatability window. In a sensitivity run, a 3x change in rate constant moved Mt Weld recovery by about one point, while grind size moved it by six.
- SX separation factors (typical P507 chloride-system values) and constant distribution ratios.
- Electrolysis at 10 kWh/kg and 75 % current efficiency; Bond work indices.
- Yttrium is lumped with Sm to Lu at the Sm molar mass.
