# JUICE — Jülich Utility for Integrating Chemical Equations

A Python simulation framework for SAPHIR chamber chemistry, replacing the legacy MCM/EASY + FACSIMILE Fortran workflow with a modern Python/SciPy ecosystem. JUICE reads SAPHIR `.enz` environment files together with MCM `.easy` / `.fac` mechanism files, builds a stiff ODE system from user-selected reactions, solves it with `scipy.integrate.solve_ivp` (BDF/Radau), and provides an interactive viewer with process analysis.

## Requirements

- **Python ≥ 3.11** (tested on 3.14.3)
- Dependencies listed in [requirements.txt](requirements.txt): `numpy`, `scipy`, `pandas`, `matplotlib`, `PyQt6`

## Installation

```bash
git clone https://github.com/GGkatzelis/MCM_SAPHIR_Tool.git
cd MCM_SAPHIR_Tool
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS / Linux
source .venv/bin/activate

pip install -r requirements.txt
```

## Quick start

### GUI scenario builder
```bash
python saphir_gui.py
```
Opens the PyQt6 Scenario Builder. Workflow:
1. **Step 1** — Load chemistry: either an MCM `.fac` file (converted on the fly) or an existing EASY file, or build ENZ inputs from SAPHIR archive data.
2. **Step 2** — Select mechanism files and constrained variables (H2O, T, NO, O3, photolysis rates, …).
3. **Step 3** — Configure injections, custom reactions, ENZ interpolation mode (linear / piecewise).
4. **Run** — writes `scenario_config.json`, invokes `scenario_runner.py`, and launches `results_viewer.py` on completion.

### Headless run
```bash
python scenario_runner.py "MCM_tested_scenarios/IDL Input/ethanol/scenario_config.json"
```
Writes `simulation_results.csv`, `simulation_kinetics.npz`, `simulation_reactions.csv` next to the config, then launches the viewer.

### Benchmark suite (validates against FACSIMILE)
```bash
python run_benchmark.py
```
Runs all test scenarios (`scene_01`, `scene_02`, `ethanol`, `mcm_full`) in both linear and piecewise interpolation modes, compares every species against the `.cs` FACSIMILE reference, prints agreement tables.

## Repository layout

```
.
├── saphir_gui.py              # GUI scenario builder (main entry point)
├── scenario_runner.py         # Headless ODE executor
├── results_viewer.py          # Interactive viewer + process analysis (BudgetEngine)
├── easy_parser.py             # FACSIMILE .easy / .fac / .mech / .rate parser
├── idl_parser.py              # SAPHIR IDL-workflow .easy → JSON converter
├── mcm_converter.py           # Raw MCM .fac → EASY include files
├── mcm_kinetics.py            # KMT / Troe fall-off functions, COMPRC evaluator
├── load_saphir_environment.py # .enz file loader → scipy interpolators
├── mechanism.py               # Reaction / Mechanism dataclasses
├── mechanism_solver.py        # Legacy stub (migrated into scenario_runner.py)
├── run_benchmark.py           # Batch regression runner vs FACSIMILE
├── compare_vs_facsimile.py    # Single-scenario comparison CLI
├── diagnostic_comparison.py   # Side-by-side timeseries plot tool
├── mcm_species_db.py          # MCM species catalogue + formula lookup
├── enz_builder.py             # ENZ file builder dialog (GUI component)
├── architecture.md            # Internal design notes
├── mechanisms/                # base/<MCM download, unmodified, manifest+checksums>  addons/<EASY include sets, manifest>
├── MCM_database/              # legacy test subsets (the MCM v3.3.1 download moved to mechanisms/base/MCM_v3.3.1)
└── MCM_tested_scenarios/      # Reference scenarios + FACSIMILE .cs outputs
    └── IDL Input/
        ├── scene_01/  scene_02/  scene_03/   # Minimal validation cases
        ├── ethanol/                          # Single-VOC chamber day
        └── mcm_full/                         # Full-MCM (5842 species)
```

## Validation status

JUICE is continuously regression-tested against FACSIMILE. Current agreement on the shipped scenarios (median ratio within ±2% counts as "matched"):

| Scenario | Species within ±2% (linear / piecewise) | Notes |
|---|---|---|
| scene_01  | matched / matched | T=const, injection + dilution |
| scene_02  | matched / matched | T step, thermal rescaling |
| scene_03  | 15/15 / 15/15     | chemistry-source reactions |
| ethanol   | 9/38 / 12/38      | baseline HO2 −7.6% propagates to peroxides; see `architecture.md` |
| mcm_full  | ~20% / ~20%       | long-tail 5-30% from RO2 coupling + solver tolerance |

See [architecture.md](architecture.md) for the physics, validation deep-dive, and known limitations.

To run and compare your own scenario, follow [docs/validation_with_claude.md](docs/validation_with_claude.md) — a step-by-step playbook (also pasteable directly into Claude Code or claude.ai to have it run the comparison and return a markdown report).

## Citation

If you use JUICE in a publication, please cite the repository and contact the author for the corresponding methods reference.

## License

[MIT](LICENSE) — free to use, modify, and redistribute with attribution.

## Contact

- Georgios Gkatzelis, Forschungszentrum Jülich GmbH — [g.gkatzelis@fz-juelich.de](mailto:g.gkatzelis@fz-juelich.de)
- Franz Rohrer, Forschungszentrum Jülich GmbH — [f.rohrer@fz-juelich.de](mailto:f.rohrer@fz-juelich.de)

## Measured chamber conditions from your own file

JUICE needs the physical conditions of the chamber over time. At SAPHIR they are read from the
archive; anywhere else, choose **Measured — from a file you provide** on the *Where and when*
step and give one table (CSV or Excel, first sheet). `Template …` writes an example with the
recognised columns, units and time formats; `templates/measured_inputs_template.csv` is the same file.

| column | aliases | unit in the file | notes |
|---|---|---|---|
| time | utc, datetime, hour | `YYYY-MM-DD HH:MM` (UTC), `HH:MM` of the day, hours of day, unix s or s since 2000 | first column if none matches |
| T | temp, temperature, TK | K (°C recognised when < 150) | chamber air temperature |
| p | press, pressure | hPa (Pa when > 5000) | optional |
| H2O | water, H2O_percent | vol-% (ppm when > 100) | or give `RH` in % → converted with T and p |
| flow | dilute, total_flow | m³/h | dilution / replenishment flow |
| roof | hv, lights, open | 0 / 1 | optional; replaces the typed roof times |
| jNO2 jO1D jHONO jH2O2 jHCHO_r jHCHO_m | | 1/s | any subset; the rest is modelled |
| NO NO2 O3 | | ppb | only to hold a species to its measurement |

Unknown columns are ignored and listed. Missing conditions fall back to the typed constants, and
the review's *Needs attention* says so. The *What goes into the input files* list states for every
input whether it is measured, a constant, or modelled. Measurements for comparison in the results
viewer are a separate file, `<date>.observations.csv` next to the result (`hour_utc,key,value`).
