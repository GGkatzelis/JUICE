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

New here? [docs/tutorial.md](docs/tutorial.md) walks through a chamber experiment, the same under UV
lamps, a field case and a study, step by step, without archive access.

Prefer code? [docs/juice_from_code.md](docs/juice_from_code.md) runs the same experiments from Python
(`import juice`), and every experiment the windows prepare also gets a `<name>_run.py` that rebuilds
and runs it. Examples: [examples/](examples).

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
`B_<species>` (background the air is diluted towards), `C_<species>` (held to the measurement) and
`O_<species>` (measured only for comparison), each with an optional unit in brackets. Three kinds of run share the table:

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
5. **Run**: input files, solver, results viewer. The `O_` columns are written next to the run in the viewer's
   observation format and the log reports bias, NMB, RMSE and R² against them (`<name>_skill.csv`).
   "Many tables as a study…" runs the case with a folder of input tables (one per trajectory or day),
   optionally over NOx and VOC emission factors, in the study window: a field study (`field_case` and
   `tables` in the `.study.json`; `python study_runner.py` runs it from the command line as well).

Three synthetic examples, one per kind of run, are in [field_examples/](field_examples/) with a commented
table each ("Load an example" copies one into your project folder; the README there explains the columns).

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
python run_benchmark.py --workers 4 --solver twostep                      # ... with Verwer's TWOSTEP (no Jacobian)
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

### GECKO-A mechanisms

JUICE reads the output of the GECKO-A mechanism generator (Aumont et al.; v1.0, 2026) and writes it as
EASY includes, so a GECKO-A mechanism runs like any other chemistry — in FACSIMILE as well as in JUICE.
The importer is Franz Rohrer's (`gecko_mechanism.py`, `gecko_photolysis.py`, `gecko_notation_to_smiles.py`,
`gecko_to_easy.py`, `make_soa_condensation.py`, table `gecko_data/solarlight.phot`), taken over from the
MCM_SAPHIR_Tool branch `easy-start-gui` on 9 Nov 2026 with his notes. From a generator folder
(`dictionary.out`, `reactionswithcom.dum`, `pvap.nan.dat`, `pero*.dat`, `scheme.log`) it writes:

- `anorg.geckoa.incl` — GECKO-A's inorganic core, and `geckoa.aliases.incl` — the GECKO-A names (HO, HNO2,
  HNO4, O3P, SULF) as `ALIAS` synonyms of the EASY/MCM names; a scenario loads this core *or* the MCM one;
- `<name>.geckoa.organic.incl` — the organic chemistry with GECKO-A's kinetics 1:1 (`A·Tⁿ·exp(−Ea/T)`,
  `TBODY`, `OXYGEN`, RO₂ counter classes as settings, Troe fall-off as `GTROE(...)`, the `EXTRA` forms
  written out, `ISOM` polynomials); multi-channel duplicates merged into summed rates; species whose fastest
  loss exceeds 1e6 s⁻¹ rescaled to the MCM `KDEC` convention with branching ratios unchanged;
- photolysis: GECKO-A tabulates J at 11 solar zenith angles per label and has no formula; each label is fitted
  to the MCM form l·cosᵐ(SZA)·exp(−n·sec SZA) and written as `GPHOT(...)`, the SZA-dependent ratio to the
  measured jNO₂ or jO¹D — the same construction as `j_clear_sky.incl`; the six SAPHIR-measured channels are
  used directly;
- `<name>.geckoa.names.incl` — MCM names where the structure is unambiguous (C2H5OH, HOCH2CH2O2, PAN, …),
  otherwise `G<GECKO code>`, every code accepted as an `ALIAS`; `<name>.dil.incl`; `<name>.SOA.include` —
  the Krechmer gas↔particle↔wall partitioning for the species GECKO-A itself partitions, with C* from its
  Nannoolal vapour pressures; `<name>.species.tsv` with SMILES (rdkit).

**In the app:** the Chemistry step has *Import GECKO-A output…* next to the base-mechanism chooser. Pick the
generator folder, name its primary species (the MCM name the Injections step will use, suggested where
known), optionally a generation cut, a C* limit and the wall scheme; the converter runs in the background
and files the result under `mechanisms/base/GECKOA_<name>/` with a manifest (`kind: gecko`, the generator
files' SHA-256, the options). From then on it is a base mechanism like the MCM download: choose
"GECKO-A: …", and the page leaves two decisions — the **inorganic core** (GECKO-A's own, with its photolysis
anchored on the measured jNO₂/jO¹D; or the MCM v3.3.1 core with the SAPHIR measured-J include, GECKO-A's
names aliased, for a like-for-like comparison with an MCM run) and whether the **partitioning** include is
loaded. Everything else — injections, chamber sources, dilution, photolysis, field emissions, studies —
is written as for any experiment; the run folder gets the includes under `chemistry/`, so FACSIMILE runs
it too, and the EASY header states the generator provenance. Two things an MCM include set brings in its
`.saphir.incl` and a GECKO-A mechanism has to bring itself are written by the import and the writer: the
physics lines (`CONST`, `DENSITY`, `M = DENSITY(P,T)`, `O2`, `N2` in `geckoa.physics.incl`, loaded first —
without them every termolecular rate, GECKO-A's `TBODY`/`OXYGEN` as much as the MCM `KMT` fall-offs, is
zero and OH runs away), and the six measured-J names GECKO-A's GPHOT ratios anchor on: the ones a run does
not measure are defined as the MCM clear-sky `PHTFKT(l,m,n)` (with the HV switch and the jNO₂ scaling of
`j_clear_sky.incl`) or, in a lamp chamber, from the lamp `J<n>`. Checked on a synthetic SAPHIR day with
the ethanol generator output: with either core O₃, OH, HONO and the ethanol decay follow the MCM-subset
run of the same experiment within a few per cent. The same import works from the command line:

```
python gecko_import.py <generator folder> --id GECKOA_isoprene --voc C5H8 [--maxgen N] [--cstar-max X] [--wall krechmer|gecko]
python gecko_to_easy.py <generator folder> --template <SAPHIR .easy> --name Isopren_GECKOA --mode pure --out <EASY root>
```

The second form is Franz's original converter: a complete scenario from a SAPHIR template EASY (`--mode
pure` = GECKO-A core, `--mode mcmcore` = the template's MCM core with only the organic includes replaced).
Every generated file carries the SHA of the generator files it came from. Imported mechanisms are
git-ignored (`mechanisms/base/GECKOA_*/`): regenerate them from the generator output rather than share
them. Verified 9 Nov 2026 on the isoprene mechanism (9 909 species, 44 204
generator reactions → 31 340 EASY reactions): regenerated from the raw output, the includes are byte for
byte those of Franz's run; JUICE parses them to the same 49 741 reactions; with the twostep profile the
SAPHIR day 18.07.2023 runs in 3 min 17 s (12 015 species with partitioning) and agrees with his BDF run to
0.035 % (median over species; O₃ 0.04 %, OH 0.34 %). The monoterpenes followed on 10 Nov 2026 (a
synthetic SAPHIR day, 20 ppb precursor + 5 ppb NO, MCM core, partitioning on, laptop): limonene (183 861
species / 936 583 reactions in the generator output; 213 110 species and 1 058 353 reactions in the ODE
system) imports in 1.3 min and integrates six chamber hours in 35 min, α-pinene (243 071 / 1 229 688;
284 234 species and 1 453 782 reactions in the ODE) imports in 1.6 min and takes 49 min — about 17 000
steps each, none rejected, 3.0 sweeps per step, the same step count as toluene, so the cost is linear in
the number of reactions. The toluene mechanism (118 165 species) runs a full 24-h day in 29.5 min. The
direct test of the solver against BDF on identical input is the 2-h toluene case without partitioning
(57 814 species, 348 108 reactions): TWOSTEP 2 min 26 s, BDF fast 9 h 1 min (1.64 M Jacobian terms); all
291 species above 1 ppt agree within 0.57 % of their run maximum (median 0.14 %; at the end O₃ 0.005 %,
OH 0.09 %, toluene 0.03 %). Output size is the next limit: the 24-h toluene `.cs` is 2.2 GB, the
α-pinene 6-h one 1.35 GB. GECKO-A
itself is open source (GPL-3.0-or-later on Zenodo, AGPL-3.0 in its GitLab repositories); JUICE reads its output and ships none of its code (the test fixture
`tests/fixtures/gecko_ethanol` is a small generator output folder).

### Solver profiles: fast, precise or twostep

JUICE has three validated ways to solve the same chemistry. Choose one in the Run step of the
experiment or field window, in the study window (one profile per study), or with
`--solver fast|precise|twostep` on `scenario_runner.py` and `run_benchmark.py`. The default is
`profile` in `FR_GUI.toml [solver]`. Fast and precise are scipy's BDF with the exact sparse Jacobian;
twostep is Verwer's TWOSTEP, which needs no Jacobian (see below).

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

**twostep** is the solver of the GECKO-A box model: Verwer's TWOSTEP (SIAM J. Sci. Comput. 15, 1243,
1994; Verwer and Simpson 1995), a variable-step BDF2 whose implicit stage is solved by fixed-point
sweeps y ← (w + gP)/(1 + gL) on the production/loss form instead of Newton iterations with an LU
factorisation — at least three sweeps, more (up to ten) until the iterate has converged to a tenth of
the tolerance, because an unconverged iterate pollutes the error estimate and forces needless small
steps (isoprene peroxy radicals did exactly that with a fixed count). Species that interconvert by
first-order reactions much faster than the step (CH3SO + O2 ⇌ CH3SOO2 in the DMS chemistry of the complete
MCM, the isoprene–OH adducts and their peroxy radicals) are found in the mechanism and solved as exact pairs
inside every sweep; without that, TWOSTEP sat at its minimum step on the complete MCM (since 10 Nov 2026:
complete MCM in 72 s against 555 s with fast, median deviation from FACSIMILE 0.027 %). No Jacobian is built, so the cost
per step is three to five rate evaluations and the run time grows with the number of reactions rather
than with the fill-in of the Jacobian factorisation
— the term that stops BDF in GECKO-A mechanisms (on this laptop, toluene GECKO-A without partitioning,
57 814 species / 348 108 reactions, two chamber hours: 2 min 26 s against 9 h 1 min with fast; on
F. Rohrer's machine 79 000 species: 73 min against 15 h, 118 000 species / 469 000 reactions: 2.5 h,
where BDF did not finish in 20 h; α-pinene GECKO-A, 284 234 species / 1.45 M reactions: six hours in
49 min). It restarts at
every ENZ row (implicit Euler), rtol 1e-3 with Verwer's max-norm error test, atol floor 100. Written
by Franz Rohrer ([twostep_solver.py](twostep_solver.py), an independent Python implementation of the
published algorithm, no code from boxmodel4gecko). Against FACSIMILE (9 Nov 2026, adaptive iteration):
every species above 1 ppt in scenes 04–09 within 0.72 % except ISOPBO2 in scene_08 at 1.6 %, 1.1 % with the pair solve (fast:
0.56 %), and in scenes 06, 07 and 09 *more* species within 0.1 % than fast; the median deviation from
fast over all species above 1e7 molecules/cm³ is 0.05–0.07 %.
Species below about 30 × the absolute-tolerance floor (3000 molecules/cm³ ≈ 1e-13 ppb at the default
`atol_floor` of 100) are numerical zero in its output: after every restart the two-step history of BDF2
carries the last real increment on (parasitic root 0.8 when the step doubles), by an amount the absolute
tolerance permits — so a species that should hold 3 molecules/cm³ can end at thousands (SA in scene_07/08),
and `run_benchmark.py` at its 100-molecule floor reports FAILs there, as it does for precise, whose cold
restarts show ≈ 1 × atol of the same effect. If sub-ppt species matter, lower `atol_floor`. Run time (two
benchmark workers, scenes 04–09 together): 191 s for the nine scenes against 465 s with fast; scene_08
15 s against 24 s, scene_09 75 s against 161 s on an idle machine; 3.1–4.6 sweeps per step. Its purpose
is the mechanism BDF cannot run. Per-ENZ-row mode only;
`JUICE_TS_CHECK=1` verifies at every row start
that P − L·y equals the BDF right-hand side, `JUICE_TS_ROWS=430,431` additionally integrates those rows
with BDF at rtol 1e-8 and prints the largest deviations.

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
field_examples/              three synthetic field cases (site, plume, Lagrangian parcel) with commented tables; make_examples.py writes them
field_window.py              field observations window: five steps like the chamber window (field_spec.py: the case and its table;
                             field_scenario.py: case → ENZ + EASY → scenario, emissions / dilution / held species / deposition)
results_viewer.py            results viewer; tabs in overview_tab.py, timeseries_tab.py
study_window.py, study.py    studies; study_runner.py runs one from the command line
facsimile_import.py          reads an experiment set up for FACSIMILE (.easy + ENZ) into the experiment description
experiment_spec.py           experiment description → ENZ + EASY (easy_builder.py, easy_scene.py)
scenario_prep.py             EASY + ENZ → scenario JSON (elog_generator.py, idl_parser.py, easy_parser.py)
twostep_solver.py            Verwer's TWOSTEP integrator (profile 'twostep'): BDF2 with fixed-point sweeps, no Jacobian
gecko_to_easy.py             GECKO-A generator output → EASY includes (gecko_mechanism.py reader, gecko_photolysis.py J fits,
                             gecko_notation_to_smiles.py, make_soa_condensation.py partitioning include; gecko_data/solarlight.phot)
gecko_import.py              a generator folder → base mechanism of the registry (kind 'gecko'); the Chemistry step's import dialog
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
