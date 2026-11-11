# One run, file by file

[`architecture.md`](architecture.md) lists what each module is. This page follows **one chamber
experiment** from the moment it is defined to the moment its results are on screen: which file is called,
what it does, what it hands on, and where you can step in.

The example is [`examples/chamber/saphir_isoprene.py`](../examples/chamber/saphir_isoprene.py), the base case. The window
takes the same path: its Run step calls the same functions (Sect. 7).

## 1. The tree

```
examples/chamber/saphir_isoprene.py   (or <name>_run.py, or the window)
│
├─ build()                                   STEPS 1–4: define the experiment (nothing is written yet)
│   └─ juice.new_experiment(name, chamber)
│       └─ experiment_spec.ExperimentSpec          the experiment: conditions, chemistry, injections, chamber
│           └─ experiment_spec.apply_chamber
│               └─ chambers.chamber(<id>)          setups/chambers/<id>.chamber.json: volume, place, light, wall defaults
│
└─ main()                                    STEP 5: run
    ├─ juice.check(exp)                            the review panel's checks
    │   └─ experiment_spec.validate
    │
    ├─ ① juice.write(exp, OUT)  ──────────────────────────────────────── writes the input files
    │   └─ experiment_spec.write_experiment
    │       ├─ ensure_mechanism                    the chemistry → OUT/chemistry/
    │       │   ├─ mechanism_registry               which base, which add-ons (chemistry/mechanisms/)
    │       │   └─ mcm_converter.convert_mcm_file   cuts the subset of exp.precursors out of the full MCM
    │       ├─ enz_columns → write_enz_file         the time grid and its columns → OUT/<name>.enz
    │       │   ├─ chambers / enz_builder           solar zenith angle, roof or lights
    │       │   ├─ measured_inputs                  (source 'file') your conditions table
    │       │   └─ enz_builder                      (source 'archive') the SAPHIR data archive
    │       ├─ render_easy → easy_form               the experiment as EASY statements → OUT/<name>.easy
    │       │   ├─ easy_builder.render_form
    │       │   │   └─ easy_scene.write_scene       the FACSIMILE layout of the file
    │       │   └─ lamp_photolysis                  (lamp chambers) J values of the light source
    │       └─ save_spec                            → OUT/<name>.experiment.json
    │
    ├─ ② juice.prepare(exp, OUT, solver, files=…)  ────────────────────── makes the solver input
    │   └─ scenario_prep.prepare_scenario
    │       ├─ elog_generator.generate_both         EASY with every #include expanded
    │       │                                       → OUT/<name>_juice.elog, <name>_juice.canonical.elog
    │       ├─ idl_parser.parse_idl_easy            canonical elog + ENZ → one scenario dictionary
    │       ├─ scenario_prep.solver_from_toml       the solver profile (settings.toml [solver])
    │       └─                                      → OUT/<name>.json
    │
    ├─ ③ juice.solve(json)  ──────────────────────────────────────────── integrates
    │   └─ python -m juice.core.scenario_runner <name>.json   (its own process; output → OUT/<name>_run.log)
    │       └─ ScenarioRunner.run
    │           ├─ _setup_environment               ENZ → time series (load_saphir_environment)
    │           ├─ _build_mechanism                 elog → reactions, rate expressions (easy_parser,
    │           │                                   evaluator, mcm_kinetics, mechanism)
    │           ├─ _solve                           scipy BDF (fast, precise) or twostep_solver (twostep)
    │           ├─ _export_results                  → OUT/<name>.cs
    │           └─ _export_kinetics                 (kv_pd) → <name>.kv, <name>.pd, _reactions.csv,
    │                                               _kinetics.npz
    │
    └─ ④ juice.results(cs, exp)  ─────────────────────────────────────── reads the results
        └─ juice.Results                            run_benchmark.parse_cs_file → a table; .ppb(), .hour
            (or: python -m juice.gui.results_viewer <name>.cs)
```

Each of ①–④ writes files and the next one reads only those files. That is why you can stop after any of
them, look at or change what it wrote, and go on.

## 2. Defining the experiment (Steps 1–4)

| | |
|---|---|
| **Files** | `juice/api.py`, `experiment_spec.py`, `chambers.py` + `setups/chambers/<id>.chamber.json` |
| **Reads** | the chamber profile (`setups/chambers/SAPHIR.chamber.json`), or a template from `setups/templates/` |
| **Does** | builds an `ExperimentSpec` in memory. `new_experiment` sets the chamber's defaults (volume, location, light type, wall sources, dilution); every line after it changes one field |
| **Writes** | nothing |
| **Step in** | any field. The four parts are `exp.conditions` (Step 1), `exp.base_mechanism` / `precursors` / `addons` / `initials` (Step 2), `exp.injections` (Step 3), `exp.chamber` (Step 4) |

The experiment is plain data. A field you don't set keeps the chamber's default, and the saved
`<name>.experiment.json` holds every field.

## 3. ① Writing the input files

| | |
|---|---|
| **Files** | `experiment_spec.py` (`write_experiment`), `mcm_converter.py`, `mechanism_registry.py`, `enz_builder.py`, `measured_inputs.py`, `easy_builder.py`, `easy_scene.py`, `lamp_photolysis.py` |
| **Reads** | the spec; the full MCM in `chemistry/mechanisms/base/MCM_v3.3.1/`; your conditions file or the archive, if chosen |
| **Does** | three things, in this order: the chemistry, the ENZ, the EASY |
| **Writes** | `chemistry/`, `<name>.enz`, `<name>.easy`, `<name>.experiment.json` |
| **Step in** | edit `<name>.easy` or `<name>.enz` by hand before ② (example, off-road C) |

- **Chemistry** (`ensure_mechanism`). For the MCM, `mcm_converter` keeps only the reactions reachable from
  `exp.precursors` and writes them as EASY include files into `chemistry/`. These are the `.init.mech`
  and `.init.rate` files plus the RO2 sum (`.ro2.incl`). Add-on files are copied next to them. A lumped
  base (CRACMM2, RACM2, …) or a GECKO-A base is copied whole instead.
- **ENZ** (`enz_columns`). It contains one row per output step. The columns are the solar zenith angle,
  roof open or shut (`hv`), dilution, temperature, water, one switch column per injection, and any
  measured species or J values. With source `manual` the physical columns are your constants. With `file`
  they come from your table, and with `archive` from the SAPHIR archive.
- **EASY** (`easy_form` → `easy_builder` → `easy_scene`). Each part of the spec becomes FACSIMILE statements:
  - the includes of the chemistry;
  - the chamber processes (HONO and HCHO wall sources, background reactivity, dilution, wall
    deposition);
  - the photolysis block;
  - initial values;
  - the injections;
  - your rate changes and verbatim statements.

  `easy_scene` lays the file out the way the IDL tools did, so FACSIMILE reads it unchanged.

## 4. ② Making the solver input

| | |
|---|---|
| **Files** | `scenario_prep.py`, `elog_generator.py`, `idl_parser.py`, `settings.toml` |
| **Reads** | `<name>.easy`, `<name>.enz`, the include files (run folder, `chemistry/`, `exp.include_dirs`, the settings file's `include_dir`, in that order) |
| **Does** | expands every `#include` into one flat file, the elog, as FACSIMILE's preprocessor does. Then parses the elog and the ENZ into one dictionary: species, reactions, rate expressions, parameters, the time series and the solver settings |
| **Writes** | `<name>_juice.elog` (expanded), `<name>_juice.canonical.elog` (normalised names), `<name>.json` |
| **Step in** | the solver profile (`solver=` here); the include search order; the elog shows exactly what the solver will see, so look there first when a reaction seems missing |

An include that is not found is marked `! NOT FOUND` in `<name>_juice.elog`, and the window's Run step
lists it.

## 5. ③ The solver

| | |
|---|---|
| **Files** | `scenario_runner.py` (`ScenarioRunner`), `twostep_solver.py`, `easy_parser.py`, `evaluator.py`, `mcm_kinetics.py`, `mechanism.py`, `load_saphir_environment.py` |
| **Reads** | `<name>.json` (it points to the elog and the ENZ) |
| **Does** | builds the reactions and their rate expressions, integrates the equations row by row of the ENZ, and writes the results |
| **Writes** | `<name>.cs` (concentrations), `<name>_run.log`; with `kv_pd`: `<name>.kv` (rate constants), `<name>.pd` (production and destruction), `<name>_reactions.csv`, `<name>_kinetics.npz` |
| **Step in** | the solver profile (below); `kv_pd=False` makes the run faster but leaves Process Analysis empty |

The run is five calls in order: `_setup_environment` → `_build_mechanism` → `_solve` →
`_export_results` → `_export_kinetics`.

The three profiles are set in `settings.toml` `[solver.profiles]` and chosen with `solver=`:

| profile | method | settings |
|---|---|---|
| `fast` (default) | scipy BDF with the exact Jacobian | rtol 1e-4, warm continuation across ENZ rows |
| `precise` | the same | rtol 1e-6, cold restart at every ENZ row |
| `twostep` | Verwer's TWOSTEP, no Jacobian (`twostep_solver.py`) | rtol 1e-3 |

`juice.solve` starts `scenario_runner.py` as its own process, as the window does. It can also be run by
hand: `python -m juice.core.scenario_runner <name>.json --solver precise`.

## 6. ④ Results

| | |
|---|---|
| **Files** | `juice/api.py` (`Results`), `run_benchmark.py` (`parse_cs_file`), `results_viewer.py` |
| **Reads** | `<name>.cs` (and `.kv`, `.pd` for the viewer's Process Analysis) |
| **Does** | reads the columns into a table. `ppb()` converts with the model's air density (column `M`) |
| **Step in** | anything you do with the numbers. `r.df` is a pandas DataFrame |

## 7. The same path from the window

The windows don't have a separate engine. The chamber window's Run step (`juice_gui.py`, `RunPage`) does
the same work as the script:

| window | in the script | function behind it |
|---|---|---|
| Steps 1–4 | `build()` | the fields of `ExperimentSpec` (`experiment_wizard.py` holds the pages) |
| review panel | `juice.check` | `experiment_spec.validate` |
| *1. Prepare input files* | `juice.write` + `juice.prepare` | `write_experiment` + `scenario_prep.prepare_scenario`; also writes `<name>_run.py` |
| *2. Run simulation* | `juice.solve` | `scenario_runner.py` in its own process (`QProcess`) |
| *3. Open results* | `juice.results` | `results_viewer.py <name>.cs` |

A **field case** (`field_window.py`, `juice.new_field_case`) differs only in ①. `field_spec.check` reads
its one table, and `field_scenario.write_case` writes the ENZ and EASY: emissions, background, held
species, dilution, deposition. Then ②–④ are the same files and functions. After ④, `juice.skill`
(`field_scenario.skill_table`) compares the model with the table's `O_` columns → `<name>_skill.csv`.
The walk-through is [`examples/field/field_site.py`](../examples/field/field_site.py); the lamp chamber's is
[`examples/chamber/uv_lamp_chamber.py`](../examples/chamber/uv_lamp_chamber.py).

## 8. I want to change …

| … | where | how |
|---|---|---|
| a condition, an injection, a chamber process | the spec, Steps 1–4 | a field in `build()` |
| a rate constant | the spec | `exp.reaction_adjustments` (the example, Step 2) |
| reactions of my own | an include file | `exp.include_dirs` + `exp.extra_includes` (off-road B) |
| a single EASY statement | the EASY | `exp.overrides`, written after the mechanism (off-road B) |
| the conditions over time | the ENZ | `c.source = "file"` with your table, or edit `<name>.enz` between ① and ② |
| holding NO2 or O3 to a measurement | ENZ + EASY | `c.constrain_species`, with the column in your table |
| the solver or its tolerance | `<name>.json` | `solver=` in ②; the numbers of each profile in `settings.toml` |
| a new chamber | the chamber profile | the chamber editor, or a `setups/chambers/<id>.chamber.json`; then `new_experiment(chamber="<id>")` |
| the whole design (loops, scans) | your script | call `build()` with arguments and loop (off-road A); for large grids `study_runner.py` |

The first rows change the spec. A change there is saved in `<name>.experiment.json` and repeated by
`<name>_run.py`. Edits to the written files (`.easy`, `.enz`) after ① are not. They live only in that run
folder, and running ① again overwrites them.
