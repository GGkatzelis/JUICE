# Field examples

Three synthetic field cases, one per kind of run, to show how a case and its input table are built.
Every number is invented (plausible for Los Angeles in August); nothing here is a measurement.
Load one from the field window ("Load an example" on the Field case step) or from the start menu:
the case and its table are copied into your project folder, and the copy is yours to change.
`make_examples.py` writes all three; a test checks the shipped files are what it writes.

| folder | kind | what it shows |
|---|---|---|
| `site_day/` | Constrained site | A 24-h diurnal cycle: NO, NO2, O3, CO, isoprene and toluene held to the table (`C_` columns), one spin-up day, a 24-h lifetime for everything not held. HCHO is not held but observed (`O_HCHO`), so the Run step reports the model's skill against it and the viewer shows both. |
| `plume_ageing/` | Plume or transect | 3 h downwind of a source. The dilution rate comes from the decay of the tracer CO above its background (`dilution = tracer`), and CO itself is held to its column; `B_` columns give the air mixed in; the fresh plume is the case's initial values. |
| `la_parcel/` | Lagrangian air parcel | 4 h from Santa Monica Bay to Pasadena: lat/lon per row, boundary-layer height `h`, emissions in mol/s over an 8 × 8 km footprint (`E_` columns, `area_km2 = 64`), growth of the layer as dilution towards the `B_` background, deposition velocities for O3 and HNO3. |

## The two files of a case

`<name>.field.json` holds everything typed in the window: kind of run, site, chemistry (an MCM subset
for the listed precursors, or an existing mechanism include set), output step, emission unit and
footprint, dilution source and rates, spin-up, lifetime, initial values, deposition. The `table`
entry names the CSV, relative to the JSON when they sit together.

`<name>.csv` is the input table: one row per time, UTC. The comment lines at the top of each example
explain every column; in short

| column | meaning |
|---|---|
| `time` | `YYYY-MM-DD HH:MM` (UTC), or `HH:MM` / hours of day with the case's date, or seconds |
| `lat`, `lon` | position per row, or once in the case |
| `T`, `p`, `RH` or `H2O` | conditions (K, hPa, %, vol-%; °C and Pa are recognised) |
| `h` | boundary-layer height, m |
| `kdil` | first-order dilution rate, 1/s (when `dilution = table`) |
| `jNO2`, `jO1D`, `jHONO`, `jH2O2`, `jHCHO_r`, `jHCHO_m` | measured photolysis frequencies, 1/s |
| `E_<species> [unit]` | emission |
| `B_<species> [unit]` | background the air is diluted towards |
| `C_<species> [unit]` | measured value the species is held to |
| `O_<species> [unit]` | measured value only compared with the model (results viewer, skill table) |

Species names are those of the mechanism (MCM: `C5H8` isoprene, `TOLUENE`, `C2H5OH` ethanol, …).
Emissions are held from one row to the next; everything else is interpolated to the output step.
