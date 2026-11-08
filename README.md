# JUICE — Jülich Utility for Integrating Chemical Equations

Box-model simulations of atmospheric simulation-chamber experiments with the Master Chemical
Mechanism (MCM). JUICE replaces the IDL/EASY + FACSIMILE workflow used at the SAPHIR chamber
(Forschungszentrum Jülich) with Python:

- **Describe the experiment, not the files.** A step-by-step window takes the date and chamber
  conditions (from the SAPHIR archive, the SAPHIR wiki log, or your own measured-conditions file),
  the precursors, the injections and the chamber effects. From that one description JUICE writes the
  ENZ and EASY input files exactly as the IDL workflow did, extracts the MCM subset, and runs it.
- **FACSIMILE-equivalent solver.** A stiff BDF solver with an analytic Jacobian and a Numba kernel
  reproduces FACSIMILE on the SAPHIR reference experiments (see [Validation](#validation)).
- **Results you can work with.** Overview, time series, process analysis (production/destruction
  budgets), OH reactivity, chemical space (volatility, oxidation state), comparison with SAPHIR
  measurements, and a figure list that becomes a PowerPoint file.
- **Studies.** Many experiments with one consistent change, or 1–2 axis sweeps (ozone isopleths,
  NOx/VOC-sensitive regimes, sensitivities).
- **Other chambers and UV lamps.** Chamber profiles (volume, dilution, walls, light source) and lamp
  photolysis from the lamp spectrum.

## Requirements

- **Python ≥ 3.11** (tested on 3.14)
- Windows is the tested platform (the GUI and launchers); the solver and tools are plain Python
- Dependencies in [requirements.txt](requirements.txt): numpy, scipy, pandas, matplotlib, PyQt6,
  numba, netCDF4, requests, Pillow, python-pptx, imageio, imageio-ffmpeg

## Installation

Windows (PowerShell):
```powershell
git clone https://github.com/GGkatzelis/JUICE.git
cd JUICE
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python experiment_spec.py --selftest; if ($?) { python check_lamp_chamber.py }   # quick sanity check
```
macOS / Linux: the same, with `source .venv/bin/activate` and `&&` instead of `; if ($?) { … }`.
The desktop launchers look for the virtual environment first as `.venv` inside the repository, then
two folders up (`..\..\.venv`, a shared venv next to several projects), and otherwise use the
`pythonw.exe` on the PATH. Everything else (reference scenarios, chamber profiles, the σ·φ photolysis
table, the mechanism include sets) is in the repository; no other download is needed.

## Quick start

### Desktop icon (Windows)
```bash
python make_shortcut.py              # the JUICE icon on the Desktop
python make_shortcut.py --startmenu  # also in the Start menu
```
The **JUICE** icon ([JUICE.bat](JUICE.bat) → [juice_start.py](juice_start.py)) first asks
whether you model a chamber experiment or field observations. The chamber path opens this menu:

| Model an experiment | Look at results |
|---|---|
| New experiment (from scratch) | Results of one experiment (`.cs`) |
| Open a saved experiment (`.experiment.json`) | Results of a study (`.study.json`) |
| New study (sweep, or several experiments with one change) | |
| Import a FACSIMILE experiment (`.easy` + ENZ) | |

- **One window per choice.** The choice opens in its own window and the menu closes. Click the icon
  again for another window at any time, while others are still open.
- **Dropping a file on the icon** opens it straight away: a `.cs` in the results viewer, a
  `.study.json` in the study window, an `.experiment.json` in the experiment window, an `.easy` or
  `.elog` in the experiment window with the FACSIMILE import open on it.
- **Shortcuts are not in git.** They hold this PC's paths, so run the script once after cloning.
  It also removes the separate "JUICE Results" icon of earlier versions.
- **[JUICE Results.bat](JUICE%20Results.bat)** still opens the viewer directly. Started without a
  file, it first asks: one experiment, or a study?

### The experiment window
```powershell
python juice_gui.py
```
The window behind "New experiment" in the start menu (`python juice_gui.py <name>.experiment.json`
opens a saved one). Five steps down the left, a review panel on the right that
always shows what the experiment contains:

1. **Experiment**: name, folder, date and time window, chamber, conditions (SAPHIR archive, your own
   measured file or an ENZ, or constants), roof or lamp schedule; or start from a template, a saved
   experiment, an imported SAPHIR wiki log, or an imported FACSIMILE experiment.
2. **Chemistry**: precursors (the MCM subset is extracted for you) or an existing include set (from
   JUICE's include folder or any include folder you add), add-on chemistry, reaction-rate adjustments.
3. **Injections**: what, how much (ppb, µl by syringe, gas syringe, mass-flow controller) and when;
   mixtures; the OH reactivity each injection adds.
4. **Chamber**: wall deposition, wall sources, background reactivity, dilution, initial values.
5. **Run**: writes the input files, runs the solver and opens the results.

You only edit the experiment description (`<name>.experiment.json`); the ENZ, EASY and scenario JSON
are generated from it and can be opened from the review panel.

### Field observations
```powershell
python field_window.py                     # or: JUICE icon → Field observations → New field case
python field_scenario.py case.field.json --run   # the same without the window
```
The field path has the same five steps as the chamber window, with one input table in place of the
chamber's archive: one row per time with the conditions (T, p, H2O or RH, lat/lon, boundary-layer
height h, kdil), photolysis columns (jNO2 …) and species columns `E_<species>` (emission),
`B_<species>` (background the air is diluted towards) and `C_<species>` (held to the measurement),
each with an optional unit in brackets. Three kinds of run share the table:

| kind | what is held / added | typical question |
|---|---|---|
| Constrained site | measured species held; everything else runs free with a generic lifetime; spin-up days | radical budget, ozone production rate, closure of HCHO or PAN |
| Plume or transect | dilution towards background (kdil from the table, a constant, or the decay of a tracer) | ageing of a plume |
| Lagrangian air parcel | emissions E/(h·A), boundary-layer growth as dilution, background, deposition | does the inventory reproduce the site? NOx/VOC sensitivity |

1. **Field case**: kind of run, name, folder, site, input table (a template is written on request); the
   conditions the model will see.
2. **Chemistry**: the chamber window's chemistry step (MCM subset or an existing mechanism include set).
3. **Emissions and constraints**: what the table emits or holds, against time.
4. **Mixing and deposition**: dilution source, background, deposition velocities or lifetimes, initial values.
5. **Run**: input files, solver, results viewer.

Everything the field adds is written as ordinary EASY statements, so the generated files run in
FACSIMILE as well: `QE_X = INPUT(QE_X)` with `QE_X --> X`, dilution through the mechanism's
`X + DIL -->` reactions with `DILUTE = CONST(KDIL + …)`, background as `B_X --> X` with
`k = CONST(KDIL*1e-9*M)`, a held species as `X = CONST(C_X*1e-9*M)`, deposition as `X -->`.
The ENZ has one row per output step (60 s by default); emissions are held from one table row to the
next, everything else is interpolated. A species held to a measurement takes the value read at the
start of each row (FACSIMILE's convention); pressure may vary along the path (`P = INPUT(P)`).
`tests/test_field_scenario.py` checks the generated files and, on a dark night, the solver against
the analytic solutions for emission, dilution towards a background, a held ramp and deposition.

### Experiments set up for FACSIMILE
```powershell
python juice_gui.py path\to\<name>.easy      # the experiment window, import dialog open on the file
```
"Import a FACSIMILE experiment" (start menu, the Experiment step, or the File menu): an experiment
modelled with FACSIMILE earlier becomes an ordinary JUICE experiment. The dialog takes the `.easy`
(an `.elog` is accepted when its `.easy` lies next to it), finds the ENZ the `FILES[ENZ]` line names,
and reads both back into the experiment description ([facsimile_import.py](facsimile_import.py)):

- the ENZ gives the date, time window and step, the roof times, and the start and duration of every
  injection (its pulse columns); it is copied into the experiment folder as
  `<name>.measured.enz` and serves as the measured-inputs file (T, humidity, flow, solar zenith
  angle, J values, NO/NO2/O3 — the "from a file" source on the Experiment step accepts an ENZ);
- the `.easy` gives the mechanism include set, the chamber effects and their strengths, the
  J-value source, pressure, initial values, the injection channels and amounts, delayed-release
  chains, every `k[...] = CONST(...)` as a reaction adjustment, and whatever else as expert
  overrides (verbatim).

The dialog shows what it found and checks itself: the EASY regenerated from the description is
compared statement by statement with the original and the ENZ columns column by column; what
could not be read back exactly is listed and kept with the experiment under "Needs attention".
Then everything is on the steps as for any other experiment — look, change, prepare, run. The
include files stay where they are: the experiment folder and its `chemistry\` subfolder are added to
the experiment's **include folders** (Chemistry step; add more when the mechanism lives elsewhere),
and "Prepare" lists any include file still not found. When the FACSIMILE results (`<name>.cs`) lie
next to the inputs, the Run step offers **Compare with FACSIMILE** after the run (±2 %, as the
benchmark; the full report goes to the experiment folder).

JUICE reads the EASY version that made the run from the `.fac` / `.out` next to it and keeps it with
the experiment, because EASY 5.96b and 6.00 differ in one rule: a `FUNCTION` macro defined more than
once (as `PHTFKT` is in a SAPHIR scenario: clear-sky first, jNO2-scaled later) applies from its own
line on in 5.96b, while in 6.00 the last definition applies everywhere. The J values differ by about
40 % between the two, so the run follows the rule of its EASY. Without a `.fac` / `.out`, JUICE uses
the 6.00 rule, as for its own experiments.

### Results without the window
```powershell
python results_viewer.py path\to\<name>.cs
```
or drop a `.cs` file on the **JUICE** icon. Without a file, the viewer first asks whether you want
the results of one experiment or of a study (a study opens in the study window).

### Command line
```powershell
python scenario_runner.py path\to\<name>.json --toml FR_GUI.toml          # one prepared run
python study_runner.py path\to\<name>.study.json --workers 4 --yes        # a study
python scenario_runner.py path\to\<name>.json --solver precise            # the same run, precise profile
python run_benchmark.py --workers 4                                       # JUICE vs FACSIMILE, scene_01..09
python run_benchmark.py --workers 4 --solver precise                      # ... with the precise profile
python compare_juice_facsimile.py scene_04 --tol 0.02 --no-stop --juice-dir "MCM_tested_scenarios\JUICE archive\JUICE\benchmark\fast\scene_04"
```
`<name>.json` is the scenario the Run step writes next to the experiment's EASY and ENZ files. The
benchmark writes its runs to `MCM_tested_scenarios\JUICE archive\JUICE\benchmark\<profile>\` (not in
git); the last line compares one of them again and lists every column.

### Tests
```powershell
pip install -r requirements-dev.txt
python -m pytest                  # 62 tests, about 1.5 minutes
python -m pytest -m "not slow"    # input-side tests only, without solver runs
```
The tests check the MCM download against its checksums, the EASY round-trip and builder for every
reference scene, the experiment-description goldens (scene_05 and scene_09 rebuilt from a
description give the IDL-written files), the chamber profiles, the validation notes (dark
experiments, plot labels), the EASY-version rule for redefined `FUNCTION` macros, scene_04 against
FACSIMILE, and the UV-lamp chamber. They also run on GitHub for every push to `main`.

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

## Validation

`python run_benchmark.py` takes every reference scene in `MCM_tested_scenarios/IDL Input` that has
the IDL-written `.easy`, its `.enz` and the FACSIMILE result `.cs`, runs it through the full JUICE
pipeline with the validated solver settings, and compares **every column** of the result with
FACSIMILE: the maximum relative deviation over the whole run, not a median. Concentrations are
measured against max(max|FACSIMILE|, 100 molecules/cm³), the solver's absolute tolerance, so species
the solver does not resolve in either program are not reported as 100 % off.

Result of 3 November 2026 (tolerance ±2 %, 4 solver processes in parallel):

| scene | settings within ±2 % | species within ±2 % | worst species | species within ±0.1 % | run time |
|---|---|---|---|---|---|
| scene_04 | 61/61 | 30/30 | NO3 0.24 % | 26/30 | 28 s |
| scene_05 | 78/78 | 57/57 | NO3 0.12 % | 55/57 | 28 s |
| scene_06 | 107/107 | 747/747 | C531O2 0.35 % | 670/747 | 66 s |
| scene_07 | 92/92 | 1162/1165 | C924O2 2.6 % ¹ | 1098/1165 | 61 s |
| scene_08 | 95/95 | 612/612 | SA 0.31 % | 603/612 | 53 s |
| scene_09 | 110/110 | 3225/3228 | NC4PDCO2H 6.7 % ¹ | 1750/3228 | 3.8 min |
| scene_01 | 5/8 | 0/2 | Y 8.4 % ² | 0/2 | 10 s |
| scene_02 | 3/8 | 0/2 | Y 8.4 % ² | 0/2 | 10 s |
| scene_03 | 30/56 | 11/22 | O 91 % ² | 6/22 | 26 s |

Scenes 04–09 are SAPHIR experiments (one ENZ row per minute, as the archive and the step window write
them): every setting (J values, rate constants, T, M, dilution) agrees, and so does essentially every
species.
¹ The three species outside ±2 % in scene_07 and in scene_09 peak at a few molecules per cm³.
² scene_01–03 are small hand-written tests whose ENZ has only a few rows (lights on at 3600 s, the
next row an hour later) while results are written every 10 s. JUICE applies each new row one output
step (10 s) later than FACSIMILE does; that shift alone produces these deviations. It is the first
known limitation in [architecture.md](architecture.md), section A6.

The table is the **fast** solver profile, the default.

### Solver profiles: fast or precise

JUICE has two validated ways to solve the same chemistry. Choose one in the Run step of the
experiment window, in the study window (one profile per study), or with `--solver fast|precise` on
`scenario_runner.py` and `run_benchmark.py`. The default is `profile` in `FR_GUI.toml [solver]`.

| | **fast** (default) | **precise** |
|---|---|---|
| Numerics | rtol 1e-4, warm BDF continuation across ENZ rows | rtol 1e-6, cold restart at every ENZ row |
| Run time | 1× (scene_09: 3.8 min) | about 2× (scene_09: 7.2 min) |
| O3, NO, NO2, OH, HO2 vs FACSIMILE | within 0.25 % | within 0.03 % in scene_06–09 |
| Every species above 1 ppt | within 0.6 % | within 0.35 % (scene_09: 0.06 %) |
| Species below 0.01 ppt | at most 7 % off (all below 1e-6 ppt) | rare large *relative* excursions, up to 5× (≤ 5·10⁴ molecules/cm³, ≈ 0.002 ppt) |
| Use for | O3, SOA and product studies, sweeps, isopleths | validation against FACSIMILE, reference runs |

Benchmark of 4 November 2026, scenes 04–09; "within" = maximum relative deviation over the whole
run. Both profiles give the same answer for anything a chamber instrument can measure (~1 ppt
detection limit). In the small mechanisms (scene_04/05) fast is even slightly tighter. Under precise,
`run_benchmark.py` reports scene_04 as FAIL: sulfuric acid (SA), peaking far below 1 ppt, is 96 % off.
Compare runs only within one profile: the difference between the profiles (up to 0.5 %) is larger
than the difference between two runs with the same profile.

Lamp chambers have no FACSIMILE reference; `python check_lamp_chamber.py` is their gate: a fixed
experiment in the generic UV-lamp chamber (NO₂ + H₂O₂ in the dark → lights on → ethanol → lights off)
is run headless and checked against physics that must hold (every J = 0 in the dark, dilution rate,
extended photostationary state jNO₂[NO₂] = [NO](k[O₃]+k[HO₂]+k[RO₂]), H₂O₂ loss = J3 + k[OH], Ox
production only while lit, N₂O₅ at night) and against the pinned result in
`MCM_tested_scenarios/UV lamps/` (every species within ±2 %; `--update` accepts a deliberate change).

See [architecture.md](architecture.md) for how JUICE is built, the solver conventions, the known
limitations, and the development history.

## Repository layout

```
juice_start.py               start menu (the JUICE icon): opens one of the windows below
juice_gui.py                 experiment window; step pages in experiment_wizard.py
field_window.py              field observations window: five steps like the chamber window (field_spec.py: the case and its table;
                             field_scenario.py: case → ENZ + EASY → scenario, emissions / dilution / held species / deposition)
results_viewer.py            results viewer; tabs in overview_tab.py, timeseries_tab.py
study_window.py, study.py    studies; study_runner.py runs one from the command line
facsimile_import.py          reads an experiment set up for FACSIMILE (.easy + ENZ) into the experiment description
experiment_spec.py           experiment description → ENZ + EASY (easy_builder.py, easy_scene.py)
scenario_prep.py             EASY + ENZ → scenario JSON (elog_generator.py, idl_parser.py, easy_parser.py)
scenario_runner.py           the solver
chambers.py, chambers/       chamber profiles
lamp_photolysis.py, photolysis/   photolysis frequencies under lamps (σ·φ table and its build script)
mechanism_registry.py, mechanisms/   MCM v3.3.1 as downloaded (checksums) and add-on chemistry
mcm_converter.py             MCM subset → EASY include files
run_benchmark.py, compare_juice_facsimile.py, check_*.py, tests/   checks
templates/, mixtures/        experiment templates, injection mixtures
MCM_tested_scenarios/        reference scenes with their FACSIMILE results, include files,
                             UV-lamp reference; JUICE archive/ holds your own runs (not in git)
docs/                        chamber survey, lamp-photolysis literature and practice
architecture.md              how JUICE is built
FR_GUI.toml                  settings: the validated solver set; your own paths in FR_GUI.local.toml
```

## Citation and collaboration

If you use JUICE in a publication, please cite it. GitHub's "Cite this repository" button gives the
reference from [CITATION.cff](CITATION.cff). Please also contact us for the corresponding methods reference.

JUICE is new research software. If you plan to use it for a publication, we strongly encourage you to
contact us early. We are happy to help check that the setup and the results are correct. Where we
contribute substantially to the analysis, co-authorship can be discussed in line with standard
authorship guidelines. This is a request, not a condition of the license.

## License

Copyright (C) 2026 Georgios Gkatzelis and Franz Rohrer, Forschungszentrum Jülich GmbH.

JUICE is free software: you can redistribute it and/or modify it under the terms of the
[GNU General Public License](LICENSE) as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version. It is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Contact

- Georgios Gkatzelis, Forschungszentrum Jülich GmbH — [g.gkatzelis@fz-juelich.de](mailto:g.gkatzelis@fz-juelich.de)
- Franz Rohrer, Forschungszentrum Jülich GmbH — [f.rohrer@fz-juelich.de](mailto:f.rohrer@fz-juelich.de)
