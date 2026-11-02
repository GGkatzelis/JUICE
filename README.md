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

### Desktop icons (Windows)
```bash
python make_shortcut.py              # JUICE + JUICE Results icons on the Desktop
python make_shortcut.py --startmenu  # also in the Start menu
```
**JUICE** opens the window for setting up and running experiments ([JUICE.bat](JUICE.bat)).
**JUICE Results** opens the results and plotting viewer ([JUICE Results.bat](JUICE%20Results.bat)).
Double-click it to pick a result, or drop a `.cs` file on it. Each click opens a new, independent
window. Shortcuts hold this PC's paths, so they are not in git. Run the script once after cloning.

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

Lamp chambers have no FACSIMILE reference; `python check_lamp_chamber.py` is their gate: a fixed
experiment in the generic UV-lamp chamber (NO₂ + H₂O₂ in the dark → lights on → ethanol → lights off)
is run headless and checked against physics that must hold (every J = 0 in the dark, dilution rate,
extended photostationary state jNO₂[NO₂] = [NO](k[O₃]+k[HO₂]+k[RO₂]), H₂O₂ loss = J3 + k[OH], Ox
production only while lit, N₂O₅ at night) and against the pinned result in
`MCM_tested_scenarios/UV lamps/` (every species within ±2 %; `--update` accepts a deliberate change).

See [architecture.md](architecture.md) for the physics, validation deep-dive, and known limitations.

To run and compare your own scenario, follow [docs/validation_with_claude.md](docs/validation_with_claude.md) — a step-by-step playbook (also pasteable directly into Claude Code or claude.ai to have it run the comparison and return a markdown report).

## Citation and collaboration

If you use JUICE in a publication, please cite it. GitHub's "Cite this repository" button gives the
reference from [CITATION.cff](CITATION.cff). Please also contact us for the corresponding methods reference.

JUICE is new research software. If you plan to use it for a publication, we strongly encourage you to
contact us early. We are happy to help check that the setup and the results are correct. Where we
contribute substantially to the analysis, co-authorship can be discussed in line with standard
authorship guidelines. This is a request, not a condition of the license.

## License

[MIT](LICENSE) — free to use, modify, and redistribute with attribution.

## Contact

- Georgios Gkatzelis, Forschungszentrum Jülich GmbH — [g.gkatzelis@fz-juelich.de](mailto:g.gkatzelis@fz-juelich.de)
- Franz Rohrer, Forschungszentrum Jülich GmbH — [f.rohrer@fz-juelich.de](mailto:f.rohrer@fz-juelich.de)

## Chambers other than SAPHIR

The chamber is a profile, `chambers/<id>.chamber.json`: volume, leak flow (the `DILUTE` expression is built
from them), wall material, light source (`sun` with latitude/longitude, `lamps` with their jNO2, or `dark`),
a data archive JUICE can read (SAPHIR only, so far), the default chamber effects, and per-species wall
losses. Choose the chamber on *Where and when*; picking another one applies that profile's switches and
`DILUTE` to the *Chamber* step, loading a saved experiment never does. `chambers/` ships SAPHIR (identical
to the built-in profile, checked by `python chambers.py --selftest`) and ~18 chambers from
[docs/chamber_survey.md](docs/chamber_survey.md) with public numbers, marked `verified: false`. *Edit …* next
to the chamber list changes a profile (volume, leak flow, light source, archive, default effects, wall losses,
injection standards), *New …* adds your own chamber pre-filled with a copy of the selected one; both write the
JSON file. SAPHIR opens read-only — *Save as new chamber …* derives another chamber from it.

**Lamp chambers** (profile light type `lamps`) do not use the sunlight parametrisation. jNO2 comes from the
measured column of your conditions file or from the profile's jNO2 at full power × the lights schedule
(a `lights` column may be 0–1 for dimmed lamps); every other photolysis frequency is
`J<n> = jNO2 × ∫σφF dλ / ∫σφF dλ (NO2)` with F the relative lamp spectrum — the method of Carter/SAPRC,
of the MCM's own chamber evaluations and of F0AM (see [docs/lamp_photolysis_practice.md](docs/lamp_photolysis_practice.md)).
σ·φ at 298 K for the MCM J numbers is `photolysis/sigma_phi_298K.csv`, built from the NCAR tuv-x data
(JPL/IUPAC) by `photolysis/build_sigma_phi.py`. The profile gives the spectrum as a two-column file
(`light.spectrum`), a synthetic stand-in (`light.spectrum_builtin`: `blacklight_350`, `blacklight_365` — a
Gaussian guess; j(O1D) and j(H2O2) hinge on the lamps' 300–320 nm tail, which only a measured spectrum
settles), `solar_like` for a Pyrex-filtered Xe/Ar arc (MCM clear-sky ratios at `effective_sza`), and
measured ratios `light.j_ratios` that override everything. The generated EASY states the ratios next to
the clear-sky ones. Measured J columns (jO1D, jHONO, …) always win. Literature values per lamp type:
[docs/lamp_photolysis_literature.md](docs/lamp_photolysis_literature.md). First-order wall losses of a
profile (`wall_loss_s1`) become `k[X --> ] = CONST(k)` lines, switchable on the *Chamber* step.
`python lamp_photolysis.py --selftest` prints the ratios of the stand-ins and of every lamp profile.
Indoor chambers have no SAPHIR-like wall HONO source and blacklights give little j(O¹D), so experiments
add an OH precursor: typically ~1 ppm H₂O₂ (UCR, Caltech, GTEC; J3 photolysis), HONO from NaNO₂ + H₂SO₄
(tens to hundreds of ppb, low-NOx) or CH₃ONO with NO (ppm, high-NOx; needs an add-on mechanism). The
UV-chamber template starts with 1 ppm H₂O₂ as an initial value.

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
the review's *Needs attention* says so. Photolysis: give `jNO2` whenever you have it — every J value
the file lacks is then modelled from it (mechanism clear-sky spectrum scaled by your jNO2); without
`jNO2` the run uses the clear-sky model. The template *SAPHIR — measured conditions from a file* on
the first step starts such an experiment. The *What goes into the input files* list states for every
input whether it is measured, a constant, or modelled. Measurements for comparison in the results
viewer are a separate file, `<date>.observations.csv` next to the result (`hour_utc,key,value`).

## Studies: many experiments at once, or a parameter sweep

File → "Study: several experiments at once, or a sweep …" (or `python study_window.py`). A study takes
your saved experiments and runs them all with one consistent change — a reaction rate (expression or
× factor), an extra add-on chemistry, a chamber parameter, scaled injections, initial values — and/or
a sweep over one or two axes (e.g. initial NO 0.1 … 20 ppb log × injected VOC × 0.5 … 4) to draw an
ozone isopleth or sensitivities. Every run gets its own folder under `<study>/runs/`, with provenance
lines in its `.easy`; the base experiments are never modified. "Run pilot, then estimate" runs the
first case alone, measures it and asks before starting the rest; finished runs are skipped when a
study is repeated or extended. Results: a table and CSV of metrics (O3 maximum, Ox gain while lit,
OH mean, NOx at the end, …), an isopleth plot for two axes, lines with normalised sensitivities for
one axis. Headless: `python study_runner.py <name>.study.json --workers 4 --yes`.
