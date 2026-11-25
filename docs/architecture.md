# JUICE — Architecture

JUICE replaces the IDL/EASY + FACSIMILE workflow for SAPHIR (and now other) chamber experiments with
Python. Its contract is strict: for the same experiment JUICE writes the **same EASY and ENZ input
files** the IDL workflow writes, and its solver reproduces **FACSIMILE's results** from them.

This document describes JUICE as it is (November 2026); its sections are numbered A1 to A6. The bug hunts
and design decisions of March–April 2026 that led here are in
[development_history.md](development_history.md) (until November 2026 "Part B" of this document).

## A1. The pipeline

```
<name>.experiment.json                       what the user edits (ExperimentSpec), in juice_gui.py
        │  experiment_spec.write_experiment
        ▼
<name>.enz  +  <name>.easy                   IDL/FACSIMILE layout and naming (QNO, QO3, Q_<name>, q<NAME>)
        │  elog_generator.generate_both      #include expanded → <name>_juice.elog + canonical .elog
        │                                    (FUNCTION scope per EASY version, A4)
        │  idl_parser.parse_idl_easy         canonical elog + ENZ → scenario dict
        │  scenario_prep.prepare_scenario    + the validated [solver] keys of settings.toml
        ▼
<name>.json                                  scenario the runner reads
        │  scenario_runner.py
        ▼
<name>.cs  (+ .kv, .pd)                      FACSIMILE column format: concentrations, settings,
        │                                    rate constants and turnovers, production/destruction
        ▼
results_viewer.py                            Overview, Timeseries, Process analysis, Reactivity,
                                             Budgets, Chemical space, figure list → .pptx
```

[`docs/run_path.md`](run_path.md) follows one experiment through this pipeline function by function,
with the files each stage reads and writes and where a user can change something.

An existing IDL-written `.easy` enters at the second box and takes the same path, which is how the
reference scenes are checked. In the GUI it enters at the first box instead: `facsimile_import.py`
reads the `.easy` and its ENZ back into an `ExperimentSpec` (the ENZ becomes the measured-inputs
file, `measured_inputs.read_enz`), checks that the regenerated EASY and ENZ equal the originals, and
the experiment is then edited and run like any other. A study (`study.py`) is the same pipeline
repeated over a grid: it modifies the `ExperimentSpec`, never a generated file.

Two more entries since November 2026. The start window asks *Chamber experiment* or *Field
observations*: a field case (`field_window.py`, a `.field.json` through `field_spec.py`) enters at the
first box too, `field_scenario.py` writing its emissions, dilution, background and held species as EASY
statements plus an ENZ, so the rest of the path — and FACSIMILE — is unchanged. And a GECKO-A generator
folder enters through the Chemistry step: `gecko_import.py` files it as a base mechanism of kind `gecko`
under `chemistry/mechanisms/base/GECKOA_<name>/` (git-ignored, regenerated from the generator output), and
`experiment_spec.gecko_form` writes the physics include, the aliases, the chosen inorganic core and the
organic/partitioning includes into the `.easy` like any other include set (A4).

Two more since 10–11 November 2026. A mechanism of the F0AM chemistry library enters as a base of kind
`lumped` (`f0am_import.py`: the three F0AM files → physics, aliases, chemistry, photolysis-name and
dilution includes under `chemistry/mechanisms/base/<ID>/`, with the F0AM originals in `f0am_source/`; CRACMM2,
CRACMM1_aq, RACM2), used as a whole; `experiment_spec.lumped_form` writes its includes and defines the
J<n> (MCM numbers) and JN<k> (F0AM channels, `chemistry/photolysis/f0am_jn_clear_sky.csv`) its photolysis names refer
to, for sun, measured J or lamps. A modified MCM of the library enters as an add-on built from the
difference to F0AM's own MCM (`f0am_import.build_addon`; `chemistry/mechanisms/addons/NOAA_BB/`); add-ons may now
carry RO2 species (appended to the run's `.ro2.incl`) and JN channels. `chemistry/mechanisms/catalog.json` lists
every chemistry offered or planned; the Chemistry step shows the planned ones greyed out.

And every box of Fig. 1 is callable from Python (`import juice`, `juice/api.py`): `new_experiment` / `new_field_case`,
`write`, `prepare`, `solve`, `results`, `run`, plus `write_script`, which writes an experiment as
`<name>_run.py`. The chamber window's Run step, the wizard's Finish and the field window's Run step write
that script next to every prepared run; `tests/integration/test_juice_api.py` requires that the script's `build()`
returns the window's spec field by field and that its ENZ/EASY are byte-identical.

## A2. Modules

The Python code is one package, `juice/`, next to the data folders (`setups/chambers/`, `chemistry/mechanisms/`, `setups/mixtures/`,
`chemistry/photolysis/`, `setups/templates/`); `juice/_root.py` is the one place that knows where they are. The module names
below are file names: `juice/api.py` (what `import juice` gives), `juice/core/` (experiment and field case as
data, ENZ and EASY writers, parser, solver, runner, studies), `juice/chemistry/` (registry, species databases,
MCM converter, lamp photolysis), `juice/importers/` (F0AM, GECKO-A, FACSIMILE, SAPHIR wiki; F0AM export),
`juice/gui/` (the windows) and `juice/tools/` (checks, benchmark, shortcut). A module is started as
`python -m juice.<folder>.<module>`, e.g. `python -m juice.core.scenario_runner <name>.json`.

| Area | Modules |
|---|---|
| Step window | `juice_gui.py` (5 steps: Experiment, Chemistry, Injections, Chamber, Run), `experiment_wizard.py` (the step pages), `experiment_plots.py`, `juice_theme.py`, `loading_box.py` |
| Experiment → input files | `experiment_spec.py`, `easy_builder.py` (EASY form ↔ scene), `easy_scene.py` (lossless EASY document model and writer), `facsimile_import.py` (FACSIMILE `.easy` + ENZ → `ExperimentSpec`), `injection_calc.py` (µl / syringe / MFC → ppb, mixtures), `oh_reactivity.py` |
| Conditions | `enz_builder.py` (SAPHIR archive NetCDF → ENZ, solar zenith angle), `saphir_wiki.py` (action log import), `measured_inputs.py` (user's own conditions file), `load_saphir_environment.py` (ENZ → interpolators) |
| Chamber and light | `chambers.py` + `setups/chambers/<id>.chamber.json`, `lamp_photolysis.py` + `chemistry/photolysis/` (σ·φ at 298 K, lamp J ratios) |
| Python interface | `juice/api.py` (experiments and field cases in code, the four steps, `<name>_run.py` scripts), `examples/` |
| Chemistry | `mechanism_registry.py` + `chemistry/mechanisms/` (unmodified MCM download with checksums, add-ons, imported GECKO-A bases, converted F0AM mechanisms, `catalog.json`), `f0am_import.py` (F0AM `_AllRxns/_K/_J.m` → base of kind `lumped`; F0AM MCM variant → add-on), `mcm_converter.py` (MCM subset → EASY include files), `mcm_species_db.py`, `species_db.py`, `gecko_import.py` + Franz Rohrer's converter (`gecko_mechanism.py`, `gecko_photolysis.py`, `gecko_notation_to_smiles.py`, `gecko_to_easy.py`, `make_soa_condensation.py`, `chemistry/photolysis/`) for GECKO-A generator output |
| Parsing and evaluation | `elog_generator.py`, `idl_parser.py`, `easy_parser.py`, `evaluator.py`, `mechanism.py` |
| Solver | `scenario_prep.py`, `scenario_runner.py`, `twostep_solver.py` (Verwer's TWOSTEP, F. Rohrer) |
| Results | `results_viewer.py` (with `mcm_kinetics.py`: KRO2 and COMPRC as a fallback for its budgets), `overview_tab.py`, `timeseries_tab.py`, `plot_interaction.py`, `report_basket.py`, `report_dock.py`, `column_file_reader.py` |
| Studies | `study.py`, `study_runner.py` (command line), `study_window.py` |
| Field observations | `field_window.py` (five steps + run page), `field_spec.py` (`.field.json`: one table of E_/B_/C_/O_ columns, deposition, dilution), `field_scenario.py` (case → ENZ + EASY → scenario, observations file, skill table), `examples/field/` |
| Checks | `run_benchmark.py`, `compare_juice_facsimile.py`, `check_easy_roundtrip.py`, `check_lamp_chamber.py`, `tests/` (unit, integration, regression and window tests; [tests/README.md](../tests/README.md)) |

`mechanism_solver.py` holds only the species list and index the runner reads; the right-hand side and the
solver are in `scenario_runner.py`. The KMT rate functions once in `mcm_kinetics.py` were removed on 25 Nov 2026:
the solver never used them, and the mechanisms define their own KMT values.

## A3. Settings

`settings.toml` holds the **validated solver settings**; changing them changes results:
BDF with analytic Jacobian, Numba-compiled stoichiometry kernel, `atol_floor` 100 molec/cm³,
`piecewise` interpolation, a settings snapshot per output step (`perstep_snap`) with live settings on
top (below), Newton tolerance 1e-6 / 10 iterations. The step window never sets numerics; `scenario_prep.solver_from_toml` copies them into
every scenario JSON.

**Live settings** (`live_settings`, default on since 15 November 2026; F. Rohrer, MCM_SAPHIR_Tool branch
`easy-start-gui`, commits 2d19712 and 7f6c89e). A setting that depends on a species or on the model time,
directly or through other settings (the RO2 sum, COA, PERO*), and every rate constant that uses one, is
re-evaluated in every right-hand-side call (BDF) and every production/loss call (TWOSTEP). Everything else
stays in the per-row snapshot and the per-row rate-constant cache. This is what EASY 7.01 to 7.05 do. EASY up
to 7.00, and JUICE before this change, held these values for one ENZ row, so results depended on the row
length. `live_settings = "false"` restores that, bit-identical to the runner of 12 November in the BDF
profiles (all nine scenes, checked 15 November). Effect on scenes 04–09 at 60-s rows: O3, NOx, OH within
0.15 %, HO2 within 0.75 %, the RO2 sum within 1.2 %, single peroxy radicals up to 7 % for a few rows. The
FACSIMILE references of `validation/IDL Input` come from EASY 5.96b / 6.00 decks, so `run_benchmark` runs
with `live_settings = false` (like with like) unless `--live` is given; references for the live mode need
decks from EASY 7.05. If M itself is live, every species gets the density term c·(dM/dt)/M. The model time
is available to EASY expressions as `time` (seconds since the first ENZ line). A new way for a rate constant
to depend on the state must enter `_kdep_names` in `_solve`, or it is silently held per row. The analytic
Jacobian has no terms for the derivative of a live rate constant; Newton iterates on the true right-hand side.

**Speed for large mechanisms** (same commits): reaction kernels on `ode_threads` cores (default 4) from
`ode_threads_min_reactions` reactions on (20 000; `parallel_kernels.py`, bit-identical to one core); rate
constants of the form (row factor) × (live setting)^n updated as one array statement (`live_forms.py`); each
distinct rate text evaluated once per row; the parse of an elog of 5 MB or more kept in
`~/.juice/parse_cache` (`parse_cache.py`; the key includes the elog text and the parser sources).

**Not taken over from that branch:** its run-option table (`solver_options.py`: TOML above the scenario JSON,
rtol default by solver). JUICE keeps its profiles, where the JSON written at prepare time carries the
validated pair. `scenario_runner` accepts `--set key=value`, `--rtol`, `--threads`, `--hmax` for one run.

**Solver profiles.** `rtol` and the restart mode come as a validated pair, the profile
(`[solver] profile`, numbers in `[solver.profiles.*]`, reference values and descriptions in
`scenario_prep.SOLVER_PROFILES`):

- `fast`: rtol 1e-4 with warm BDF continuation across ENZ rows (`perstep_warm`). The default since
  2 November 2026.
- `precise`: rtol 1e-6 with a cold restart at every ENZ row, the set used before 2 November.
- `twostep`: Verwer's TWOSTEP (`solver_method = "TWOSTEP"`, rtol 1e-3, restart at every ENZ row;
  `twostep_numit` 3, `twostep_numit_max` 10, `twostep_conv` 0.1, `twostep_dtmin` 1e-4), since
  9 November 2026, for the mechanisms BDF cannot factorise. A profile may carry its own
  `solver_method`; the BDF profiles leave it unset. `jac_mode`, `perstep_warm` and `first_step` are
  ignored under TWOSTEP.

The JSON records `rtol`, `perstep_warm` and `solver_profile` (and, for twostep, `solver_method` and
the `twostep_*` keys). Choices:

- per run in the Run step, or with `scenario_runner.py --solver`
- per study, frozen in the `.study.json`; a study refuses to change it once runs have finished
- `run_benchmark.py --solver`

Rules learned on the way (benchmark of 4 November 2026, see the table in user_guide.md):

- **Never combine rtol 1e-6 with warm continuation.** It diverges in trace species: HCOCO3 1.4e5 %
  in scene_06, HO13C4O2 1.3e4 % in scene_09. `tests/unit/test_profiles.py` guards this.
- **Precise is tighter but not cleaner.** Its cold restarts re-ramp the step size at every row, and
  species below 0.01 ppt then show large relative excursions (BUTAL2NO3 5×, MGLYOO 1.9×, HO5C6CO3
  1.5×, SA ~1×; absolute size at most ~5·10⁴ molec/cm³, CH3OCHO at 0.002 ppt). The benchmark's ±2 %
  rule therefore fails scene_04 (SA) and scene_06/09 under precise, though no species above 1 ppt is
  affected. Fast has excursions too, but small: at most 7 %, all in species below 1e-6 ppt.
- **Fast costs accuracy only in large mechanisms.** In scene_06–09, main species are 0.1–0.25 % off
  under fast vs 0.02–0.03 % under precise. In scene_04/05, fast is slightly tighter.
- **Results from the two profiles differ by up to ~0.5 %.** Compare runs within one profile.
- **TWOSTEP against the others:** species above 1 ppt within 0.72 % of FACSIMILE in scenes 04–09
  (ISOPBO2 in scene_08 1.6 %) and within 0.57 % of the BDF solution on a GECKO-A toluene case
  (A5); species below ~30 × `atol_floor` are numerical zero in its output (A6).
- **A pre-profile TOML** (with `rtol_default` and no `profile` key) is used exactly as written
  (profile `custom`). Example: the 2 November file `f91383b^:settings.toml` is rtol 1e-6, cold.
- **A pre-profile scenario JSON** (with `rtol` and no `perstep_warm`) gets the TOML's profile when its
  rtol matches. Otherwise it gets a safe pairing: warm only if rtol > 1e-5.

`[paths]` is absent on purpose: JUICE then uses the folders inside the repository
(`validation/IDL Input` for references and include files, `runs/`
for new experiments). A machine-specific `settings.local.toml` (never committed) overrides
any key.

## A4. Solver conventions worth knowing

- **ENZ values are step functions** (`interpolation = "piecewise"`, as in FACSIMILE). Every ENZ row
  boundary and every output step starts a new segment; all settings (J values, rate constants, T, M,
  switches) are evaluated once at the segment start and held (`perstep_snap`).
- **Boundary convention.** An ENZ line holds from its own time to the next line; a value looked up
  *at* a row time is that row's value. (Until 17 Nov 2026 the rows were shifted by +1 ns, which only
  mattered for ENZ times near 0; A6 item 1.) This reproduces FACSIMILE in scene_01 … scene_09.
- **Thermal expansion.** At a row boundary where M changes, concentrations (and the BDF history) are
  rescaled by M_new/M_old, as FACSIMILE does.
- **O1D** is not integrated; it is written as its pseudo-steady state J1·[O3]/k_quench.
- **Last definition wins** for parameters and `k[...] = ...` overrides; reactions are never
  deduplicated, so an include set must not repeat a reaction ([development history](development_history.md), section 8).
- **A `FUNCTION` macro defined more than once** follows the EASY version that made the run
  (`elog_generator.easy_version`, read from the `.fac`/`.out` next to a FACSIMILE run;
  `function_scope_for`): EASY 6.00 applies the last definition everywhere (the canonical elog keeps
  only that one), EASY 5.96b applies a definition from its own line on (the canonical elog keeps
  every definition and `easy_parser` versions them as `NAME`, `NAME__2`, … so each use binds to the
  definition in effect at its line). `PHTFKT` is the case that matters: clear-sky in
  `<mech>.fac.saphir.incl` where `J1..J61` are defined, jNO2-scaled in `j_clear_sky.incl` and
  `saphir.jmeasured_02.incl` — about 40 % apart. scene_01/02 and campaign runs made before 2026 are
  5.96b, scene_03–09 are 6.00; JUICE's own experiments (no `.fac`) use the 6.00 rule.
- **`#include` names may contain spaces**; the name runs to the end of the line or a `;` comment.
- **TWOSTEP restarts at every ENZ row** with an implicit Euler step and integrates the row relative to
  its start; held species are frozen at the row start, and rate constants follow the same rules as in
  the BDF per-row mode (live settings; with `live_settings = false`, O1D, rate constants that name a
  species and the sources are still refreshed per call since 15 November, which moved the frozen TWOSTEP
  results closer to FACSIMILE: scene_09 3227 of 3228 species within 0.1 %).
  Debug switches (`JUICE_TS_CHECK=1` verifies at every row that P − L·y equals the BDF right-hand side,
  `JUICE_TS_ROWS=a,b` integrates those rows with BDF at rtol 1e-8 for comparison). The error test is
  Verwer's max-norm at rtol 1e-3 with the same `atol_floor`; the fixed-point iteration runs at least
  `twostep_numit` sweeps and up to `twostep_numit_max` until the change is below `twostep_conv` ×
  tolerance — an unconverged iterate pollutes the error estimate, and with a fixed count the isoprene
  peroxy radicals forced sub-second steps for a whole day. The BDF path is bit-identical with and
  without the TWOSTEP code (checked on all nine scenes after every runner change).
- **Held species** (field mode): `X = CONST(<expression of an INPUT>)` makes X a held species
  (`idl_parser.held_species`). The runner sets it per ENZ row, zeroes its derivative and drops its
  Jacobian rows (BDF) or its P/L entries (TWOSTEP), and restores it exactly after each row.
- **A GECKO-A base brings its own physics.** An MCM include set defines `CONST`, `DENSITY`, `M`, `O2`
  and `N2` in its `.saphir.incl`; a GECKO-A mechanism has no such file, so the import writes
  `geckoa.physics.incl` and the writer loads it first. Without it every termolecular rate (the MCM
  `KMT` fall-offs as much as GECKO-A's `TBODY`/`OXYGEN`) is zero and OH has no sink — the symptom to
  recognise is OH at 1e8–1e9 and HONO near zero.
- **Photolysis anchors.** GECKO-A's `GPHOT` ratios anchor on the six SAPHIR-measured J names (jo1d,
  jno2, jh2o2, jhono, jhcho_r, jhcho_m). A run that does not measure them gets them as the MCM
  clear-sky `PHTFKT(l,m,n)` under the sun or from the lamp `J<n>` in a lamp chamber — never as
  `CONST(J<n>)` under the sun, because `J1…J61` exist only with the MCM photolysis include, and an
  undefined name in a rate expression evaluates as 1 s⁻¹, not 0 (A6).
- **Two FACSIMILE habits JUICE reproduces with a warning.** A plain name in a parameter expression
  that is in no reaction (a truncated include: `COHM2CO` for `COHM2CO3` in an RO2 sum) is kept as a
  species without chemistry holding `INITIAL[*]`, as FACSIMILE does; only a missed `NAME(` call
  still stops the run (`scenario_runner._finalize_mechanism`). And FACSIMILE is case-sensitive while
  JUICE is not: a name written `XNOY` and `XNOy` is two species there and one here —
  `elog_generator.case_variants` reports such names when the scenario is prepared.

## A5. Validation

`python -m juice.tools.run_benchmark` runs every reference scene that has `.easy`, `.enz` and a FACSIMILE `.cs`
(scene_01 … scene_09) through the full pipeline and compares every column of the result with
FACSIMILE: maximum relative deviation over the whole run, concentrations measured against
max(max|FAC|, 100 molec/cm³), the solver's absolute tolerance. The current table is in [user_guide.md](user_guide.md#validation).
`python -m pytest` (about 11 min, also on GitHub Actions, with a coverage table) runs the test suite in four
layers ([tests/README.md](../tests/README.md)): unit tests that call one function with a hand-made input; the
input-side self-tests for every scene and the importers; regression tests against a reference — scene_04
against FACSIMILE, the analytic scenes `validation/scene_00_clock` and `validation/scene_T1`, a small experiment
whose seven pieces of chemistry all have closed forms (`tests/exact_small.py`: decay, equilibrium, second order,
light switch, injection, dilution, temperature step; every solver profile), and the UV-lamp chamber check, which
has no FACSIMILE reference and is checked against physics instead (photostationary state, H2O2 loss, dark
chemistry) and against a pinned run; and the windows opened without a screen. The budget engine of the results
viewer is checked on the run of that small experiment, term by term.

TWOSTEP has the same gates: `run_benchmark.py --solver twostep` (scenes 04–09, 9 Nov 2026: every
species above 1 ppt within 0.72 % of FACSIMILE, ISOPBO2 in scene_08 1.6 %; the 100-molecule floor FAILs
on floor species, A6), `tests/regression/test_twostep.py`, and the BDF path checked bit-identical. At GECKO-A
sizes the reference is the BDF solution of the same scenario: the 2-h toluene case without partitioning
(57 814 species, 348 108 reactions, 10 Nov 2026) agrees within 0.57 % of each species' run maximum for
all 291 species above 1 ppt (median 0.14 %) — 2 min 26 s against 9 h 1 min. The field mode is verified
against analytic solutions (`tests/regression/test_field_scenario.py`, slow) and against F0AM's Lagrangian core on
a dark night; the GECKO-A path by `tests/integration/test_gecko_import.py` and `tests/integration/test_gecko_front_door.py` on
the ethanol generator output in `tests/fixtures/gecko_ethanol`, including a run of the GECKO-A
mechanism on the MCM core against the MCM subset of the same experiment (within 5–10 %).

The converted F0AM mechanisms are checked by `tests/integration/test_f0am_import.py` (converter constructs on a
synthetic file, installed counts, no undefined name in any rate, the corrections of the library files,
one run each for CRACMM2 and NOAA_BB) and, outside the suite, against F0AM running its own files with
identical constant J (dark / lit, every species above 1 ppt: CRACMM2 0.39 / 0.23 %, CRACMM1 0.30 / 0.26 %
without F0AM's marine O3 sink, RACM2 0.40 / 1.0 %; RACM2's 16 K-file constants equal MATLAB's to six
digits). They have not been run in FACSIMILE yet; CRACMM writes negative yields for its DELTA_C/DELTA_N
counters, which FACSIMILE has to accept.

## A6. Known limitations (November 2026)

1. *(Fixed 17 Nov 2026.)* **One-output-step lag with ENZ files that start near 0.** scene_01–03
   (ENZ times from 0, few rows) were off by up to 8 % at injection ramps and more at light and
   temperature steps: the +1 ns shift of the ENZ breakpoints (A4) is real at small times and rounds
   away at SAPHIR's absolute times (~7.7·10⁸ s), so JUICE read each row from the next line on.
   F. Rohrer removed the shift (MCM_SAPHIR_Tool 9e4fc04, test scene `validation/scene_T1`); scene_01–03
   now pass (worst 0.91 %) and scene_04–09 are unchanged bit for bit.
2. **Trace species below ~10 molec/cm³.** In scene_07 and scene_09, three of ~1 200 and ~3 200
   species differ by 2.6–6.7 %; all peak at a few molecules per cm³, near the solver tolerance.
3. **mcm_full** (5 842 species) is not in the benchmark: its files are named `full_MCM.*`, not after
   the folder, and its 100 MB `.cs` would dominate the run time.
4. **Chamber profiles other than SAPHIR** come from public numbers (`docs/chamber_survey.md`) and are
   marked `verified: false`.
5. **The printed RO2 sum differs from FACSIMILE's by up to 1.4 % for a few minutes while it rises
   steeply** (after the roof opens or after an injection; 20–70 of ~1 000 rows per scene, median
   below 0.02 %; the same in both solver profiles and in the runs of the earlier repository). It is
   FACSIMILE's printed value, not the chemistry: in FACSIMILE's own `.cs` the RO2 column is that
   0.8–1.4 % away from the sum of its own member radicals (scene_04, 06, 09), while JUICE prints
   exactly the sum of its members at the output time, and the member sums of the two programs
   agree within 0.25 %. Nothing to correct in JUICE; the RO2 column stays inside the ±2 % tolerance.
6. **Floor species under TWOSTEP.** Species below about 30 × `atol_floor` (3 000 molec/cm³) are
   numerical zero in its output: after every restart the BDF2 history carries the last increment on by
   an amount the absolute tolerance permits (SA 3 → 2 458 molec/cm³ in scene_08, DHPMPAL in scene_07).
   `run_benchmark.py` at its 100-molecule floor therefore FAILs scenes 04–09 for twostep; at the 1-ppt
   floor all pass. Lower `atol_floor` if sub-ppt species matter.
7. **Undefined names.** In a *parameter* (a setting such as `jo1d = CONST(J1)`), an undefined plain name is
   kept as a species without chemistry holding `INITIAL[*]`, as FACSIMILE does, with a warning; with the SAPHIR
   `INITIAL[*] = 1` it reads as 1 (met as O3 stuck at 1 ppb on the GECKO-A inorganic core; the writer now avoids
   it). In a *reaction rate* it reads as 0 and the reaction does not run; since 25 Nov 2026 the runner names
   such a name in a warning (`tests/regression/test_failed_runs.py`); no SAPHIR scene has one.
10. **Rate sanity check** (25 Nov 2026). Once per row the runner checks every rate constant before it is used: a
    value that is not finite (it is set to 0, as before), or one beyond any physical value for its order (first
    order above 1e13 s⁻¹, second order above 1e-8 cm³ s⁻¹), is named in a warning and repeated at the end of the
    run; the run itself is unchanged. No warning in CRACMM 2, CRACMM 1, RACM2, NOAA_BB, GECKO-A ethanol and
    toluene, the MCM or the field examples. It does fire in scene_08 and scene_09: the photolysis of NO3CH2CHO
    in `saphir.jmeasured_02.incl` (PHOTMMSC, exp(+1.389/cos)) reaches 6e15 and 2e33 s⁻¹ at twilight, where the
    cosine is held at its 0.01 floor while jNO2 is still above 0. FACSIMILE computes the same, so the two still
    agree; the include file needs a fix, and every SAPHIR run through sunrise or sunset with it shows the warning.
8. **The runner builds the analytic Jacobian structure under TWOSTEP** (1.86 M terms for toluene)
   although the solver never uses it — wasted set-up time at GECKO-A sizes; skip it when
   `_use_twostep`.
9. **Output size above ~10⁵ species.** The `.cs` of the 24-h toluene run is 2.2 GB, the α-pinene 6-h one
   1.35 GB (more with the `.kv` export); the output and the analysis tools need a species selection at
   that size. Mechanisms of 10⁶ reactions run (α-pinene: 1.45 M reactions in the ODE system, 49 min for
   six chamber hours), but the solver is validated against BDF only up to 3.5·10⁵ reactions and against
   FACSIMILE up to the full MCM; beyond that there is no reference yet.
