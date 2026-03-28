🏗️ SAPHIR-MCM Tool: Architecture Overview

The tool is designed as a full-stack Simulation Framework that replaces the legacy MCM/EASY Fortran-based workflow with a modern Python/ODE-solver ecosystem. It bridges the gap between raw chamber measurement data (.enz files) and the Master Chemical Mechanism (MCM) kinetics, incorporating robust GUI generation, mathematical stiffness solving, and interactive visualization.
1. Data-Centric Design Pattern

The architecture follows a Decoupled Workflow:

    Front-End (saphir_gui.py): A PyQt6-based "Scenario Builder." Its sole purpose is to generate a scenario_config.json. This JSON acts as the "Chemical Matrix," containing user-selected mechanism files, coordinate mappings for environmental data, and chamber-specific artifact constants.

    The Matrix (scenario_config.json): The single source of truth. It defines what is "Measured" (dynamic interpolation) vs. what is "Fixed" (constants) and stores the "Scaling Factors" for chamber emissions.

    Back-End Engine (scenario_runner.py): A headless simulation script that consumes the JSON. It uses scipy.integrate.solve_ivp (BDF method) to solve the stiff chemical ODEs.

    Interactive Diagnostics: Upon completion, the runner automatically launches an interactive Tkinter/Matplotlib viewer, natively merging `.enz` environmental parameters, chamber injection timelines, and chemical species outputs onto overlapping twin-axis graphs.

2. Core Technical Components
A. The Physical-Chemical "Bridge"

The agent must maintain the transparency of how gas-phase mixing ratios are converted to concentrations:

    Air Density (M): Dynamically calculated using the Ideal Gas Law:
    M=kB​TP​

    All ppb inputs from the GUI are converted to molec/cm3 using this dynamic M.

    Dilution Logic: Applied as a first-order decay term for all species:
    dtd[C]​=MCM_Net−[C]×Vol×3600Flow+0.6​

B. The Artifacts Engine

Unlike standard MCM models, this tool explicitly models SAPHIR "Background" physics:

    Photolytic Sources: Parameterized emissions of HONO and HCHO from Teflon walls, scaled by measured JNO2​ and user-defined Scaling Factors.

    Missing Reactivity (Y): Implementation of a virtual species Y that reacts with OH to match observed chamber reactivity.

    Wall Dark HONO: A background thermal production term (H2​O→HONO) active in dark conditions.

3. Agent Instructions for "Takeover"

When continuing development, the Gemini agent must:

    Maintain GUI Integrity: Ensure that adding new variables (e.g., Cl2​, NOx​ sources) follows the "Active Checkbox + Scaling Factor" pattern for transparency.

    Solver Stability: Use the BDF or Radau methods in solve_ivp as the chemical systems are mathematically "stiff."

    Transparency Console: Always update the logic_preview box in the GUI when a user changes a parameter so the "Science Logic" is visible before the simulation starts.

    Unit Consistency: Strictly enforce K for Temperature, hPa for Pressure, and s−1 for photolysis frequencies.

mechanism.py: Defines the data structures (Species, Reaction).

easy_parser.py: Parses dual-syntax chemistry files (Standard FACSIMILE `% k : A + B = C ;` and SAPHIR EASY `A --> B`).

mechanism_solver.py: A legacy ODE calculator. Its core functionality has been migrated to scenario_runner.py for improved performance and dynamic scoping.

saphir_gui.py: The Builder (creates the JSON Matrix).

scenario_runner.py: The Executor (interpolates measurements and runs the solver).

## SAPHIR GUI (`saphir_gui.py`)
The `saphir_gui.py` module provides a graphical user interface (GUI) for building and configuring simulation scenarios. It is built using PyQt6 and allows users to define environmental conditions, initial concentrations, chamber effects, and chemical injections, and then save this configuration to a `scenario_config.json` file.

### Purpose and Workflow
The primary goal of the GUI is to streamline the process of creating complex simulation scenarios by providing an intuitive, interactive interface. The workflow is as follows:
1.  **Load Folder**: The user starts by selecting a root folder for the scenario, which should contain the necessary mechanism files and a SAPHIR `.enz` data file.
2.  **Configure Scenario**: The user navigates through the tabs to define:
    -   **Constrained Variables**: Maps model variables to `.enz` columns or fixed constants, and defines Derived Photolysis Frequencies (J-Value scaling).
    -   **Chamber Sources**: Configures chamber-specific artifacts (light/dark emissions, Y-background) and steady Continuous Sources (e.g., mixing influxes, water impurities).
    -   **Chamber Initialization**: Sets initial concentrations, schedules UTC-based Burst Injections, and defines Custom Reactions for ad-hoc VOC chemistry.
3.  **Active Scenario Summary**: A sidebar console provides a real-time, highly explicit summary of the configured logic, mapped variables, and active equations.
4.  **Generate JSON**: Once the configuration is complete, the user saves the settings to a `scenario_config.json` file, which can then be used by the main pipeline to run the simulation.

### Key Components
-   **`SaphirGUI(QMainWindow)`**: The main application window.
-   **Folder Selection**: A simple dialog to select the working directory.
-   **Sidebar**: Displays mechanism files, detects `.enz` files, calculates exact UTC Start/Stop ranges from data Epochs, and houses the Active Scenario Summary.
-   **Tabs**:
    -   **Constrained Variables**: Tables for environment mapping and derived J-values.
    -   **Chamber Sources**: Controls for chamber artifacts and continuous species influxes.
    -   **Chamber Initialization**: Controls for initial conditions, burst injections, and custom reactions.
-   **Active Scenario Summary**: A read-only text area that details the active state of the matrix.
-   **Scientific Help**: Provides tooltips with scientific explanations for various parameters.
-   **JSON Output**: The `save_json` function serializes the configuration into a structured JSON file.

## Bug Fixes Applied (2026-03-28)

The following root-cause bugs were identified by auditing the EASY chemistry files against
the Python ODE engine and have been corrected:

### Fix 1 — COMPRC (Critical: all Troe KMT rates were evaluating to 0.0)
**Root cause:** `COMPRC` in `anorg.mcm331.dark.rate` is defined using EASY's `FUNCTION[...]`
(square bracket) syntax. The `easy_parser` only strips `FUNCTION(...)` (round brackets), so
`COMPRC` was stored as the raw IDL template string. In the ODE context, `ctx['COMPRC']` resolved
to `CallableFloat(0.0)` (SafeDict fallback). Every `KMT0x = COMPRC(k0, ki, Fc)` expression
then evaluated to `0.0`, **overwriting** the mcm_kinetics KMT lambda functions that had been
loaded just before. Result: NO+O3→NO2, OH+NO2→HNO3, HO2+NO→OH+NO2 and most other
pressure-dependent reactions were completely disabled.

**Fix:** Added `COMPRC(k0_full, ki, Fc)` as a proper Python function to `mcm_kinetics.py`.
Since its name does not start with `KMT`, it is loaded into `ctx` directly (no T/M wrapping),
matching the call signature used by all EASY rate expressions.

**JSON change:** Removed `"COMPRC"` from `active_parameters` so the broken `FUNCTION[...]`
template string never re-evaluates and overwrites the Python callable. Also removed `KFPAN`
and `KBPAN` FUNCTION[] entries are now handled correctly via the Python COMPRC.

### Fix 2 — PHTFKT photolysis callable (all J-values were 0.0)
**Root cause:** `PHTFKT`, `PHTFKT_CS`, `PHOTMMSC` are also `FUNCTION[...]` IDL templates.
They were stored as raw strings in `active_parameters` and evaluated to `0.0` in
`var_assignments`. Every J-value expression (`J1 = PHTFKT(a,b,c)` etc.) therefore produced
`0.0`, killing all photolysis rates.

**Fix:** `scenario_runner.py` now builds `PHTFKT` as a Python closure between pass 1 and
pass 2 of `var_assignments`. Pass 1 computes `COSX`, `LNO2`, `MNO2`, `NNO2`, and the
COMPRC-based KMTs. Then a Python `_phtfkt(a, b, c)` closure is created using those values
(measured-jno2 scaling: `cs_J * min(jno2_meas/jno2_cs, 1)`), and placed in `ctx['PHTFKT']`.
Pass 2 then evaluates all J-values correctly. `PHTFKT_CS` and `AMAX`/`AMIN` are also set
at this point. `PHTFKT_CS` and `PHTFKT` removed from `active_parameters` to prevent the
FUNCTION[] template from overwriting the callable in pass 2.

### Fix 3 — DIL dilution double-counting eliminated
**Root cause:** `ethanol.fac.dil.mech` and `ethanol.fac.dil.rate` define explicit `species + DIL →`
reactions for all organic species, with rate `DILUTE`. The Python ODE already applies global
first-order dilution via `dydt -= y * k_dil` to ALL species. Including both would double-count
dilution for organic species. Furthermore, inorganic species (NO, NO2, O3) have no corresponding
EASY dil reactions, so the EASY-only approach would leave them undiluted.

**Fix:** Removed `ethanol.fac.dil.mech` and `ethanol.fac.dil.rate` from `mechanism_files`.
The Python global `k_dil = (FLOW + 0.6) / (chamber_vol * 3600)` handles all species uniformly.

### Fix 4 — Spurious active_parameters removed
**Root cause:** The `easy_parser` accidentally parsed German comment text and internal Fortran
notation (e.g., `"LEBENSDAUER GEGEN WANDVERLUST UNGEFÄHR"`, `"MIT        F(K0,KI,FC)"`,
`"UND        N(FC)"`, `"IM FORMAT: K(K0,KI,F)"`) as parameter definitions and stored them.
Intermediate Troe variables (`K`, `K0`–`K4`, `F`, `FC`, `KI`) also appeared. In `var_assignments`,
these fail evaluation and set spurious `0.0` entries in `ctx`, polluting the chemical state.
Aliases (`ALIAS[OH]` etc.) and FUNCTION[] templates (`DENSITY`, `ARRH`, `AMAX`, `AMIN`) were
also incorrectly included.

**Fix:** Removed all spurious entries from `active_parameters` in `scenario_config.json`.
Clean list now contains only valid Arrhenius/COMPRC rate constants, J-values, and constants.

### Fix 5 — Photolysis auto-linker corrected and moved
**Root cause:** The auto-linker ran **before** `var_assignments`, so its seeded values were
immediately overwritten by `var_assignments`. Additionally, the J-number mappings were wrong:
`J41 ← JO1D` (J41 is organic peroxide photolysis, not JO1D) and `J1 ← JNO2` (J1 is JO1D in
MCM numbering, not JNO2).

**Fix:** Auto-linker now runs **after** both passes of `var_assignments` as a fallback for any
J-values still 0. Corrected mappings:
- `J1 ← JO1D`   (O3 → O(1D) photolysis)
- `J4 ← JNO2`   (NO2 → NO + O photolysis)
- `J7 ← JHONO`  (HONO → OH + NO photolysis)
- `J3 ← JH2O2`  (H2O2 → 2OH photolysis)
- `J11 ← JHCHO_R`, `J12 ← JHCHO_M`

## ODE Evaluation Order (after fixes)

```
system_ode(t, y):
  1. Load env_map data (TEMP, PRESS, JNO2, JO1D, JHONO, FLOW, HV, SZ, etc.)
  2. Water auto-balancer (RH → H2O)
  3. ctx.update: T, M, O2, N2, math functions
  4. Load mcm_kinetics: COMPRC (direct), KMT lambdas (overwritten below)
  5. Load species concentrations
  --- PASS 1 var_assignments ---
  6. Evaluates COSX, LNO2/MNO2/NNO2, JSCALE, KMT01-KMT14 (via COMPRC), KFPAN/KBPAN, etc.
  --- Build PHTFKT Python callable ---
  7. Creates ctx['PHTFKT'], ctx['PHTFKT_CS'], ctx['AMAX'], ctx['AMIN']
  --- PASS 2 var_assignments ---
  8. Evaluates J1-J51 (via PHTFKT), RO2 sum, etc.
  --- Auto-linker (fallback) ---
  9. Seeds J1←JO1D, J4←JNO2, J7←JHONO, etc. if still 0
```

## Next Steps & Future Enhancements

1.  **Simulation Scenario Testing** (`simulation_scenarios.py`):
    *   Bottom-up validation: Phase 0 (O3 dilution only) → Phase 1 (dark NOx) → Phase 2 (photolysis)
        → Phase 3 (SAPHIR artifacts) → Phase 4 (ethanol). See `easy_functions_definitions.py`
        for the full testing plan.
2.  **Batch Processing Pipeline**:
    *   Create a command-line interface or a "Batch Mode" in the GUI to run multiple
        `scenario_config.json` files sequentially.
3.  **Data Export & Reporting**:
    *   Add PDF/PNG export buttons directly to the interactive Matplotlib viewer.
4.  **Deprecate Obsolete Wrappers**:
    *   `mechanism_solver.py`'s calculate function has been fully migrated natively into
        `scenario_runner.py`'s ODE closure. Legacy stubs should be safely removed.
