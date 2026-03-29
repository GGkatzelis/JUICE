# SAPHIR-MCM Tool: Architecture Overview

The tool is a full-stack simulation framework replacing the legacy MCM/EASY Fortran-based workflow with a modern Python/ODE-solver ecosystem. It bridges raw chamber measurement data (`.enz` files) and the Master Chemical Mechanism (MCM) kinetics, with GUI-driven scenario building, stiff ODE solving, and interactive visualization.

## 1. Data-Centric Design Pattern

Decoupled three-stage workflow:

- **Front-End** (`saphir_gui.py`): PyQt6 "Scenario Builder." Sole purpose is to produce `scenario_config.json`.
- **The Matrix** (`scenario_config.json`): Single source of truth. Defines Measured (interpolated) vs. Fixed variables, scaling factors, custom/excluded reactions, and active parameters.
- **Back-End Engine** (`scenario_runner.py`): Headless script that consumes the JSON, solves the stiff ODE system with `scipy.integrate.solve_ivp` (Radau method), exports CSV, and launches an interactive viewer.

## 2. Module Summary

| File | Role |
|---|---|
| `mechanism.py` | `Reaction` and `Mechanism` dataclasses |
| `easy_parser.py` | Dual-syntax parser: FACSIMILE `% k : A + B = C ;` and SAPHIR EASY `A --> B` |
| `mcm_kinetics.py` | Troe fall-off KMT functions and `COMPRC` Python callable |
| `load_saphir_environment.py` | Loads `.enz` files into scipy interpolators |
| `saphir_gui.py` | GUI scenario builder |
| `scenario_runner.py` | ODE executor |
| `mechanism_solver.py` | Legacy stub — core logic fully migrated to `scenario_runner.py` |

## 3. Core Physics

**Air density M** (Ideal Gas Law): `M = (P×100) / (kB×T) × 1e-6` [molec/cm³]

**Dilution** (first-order decay for all species): `d[C]/dt -= [C] × (FLOW + 0.6) / (V × 3600)`

**SAPHIR artifacts** (modelled explicitly):
- Photolytic HONO wall source: `4.7e13 × (1 + (RH/11.6)²) × exp(−3950/T) × JNO2 × HV`
- Photolytic HCHO wall source: `2.8e6 × HV`
- Missing OH reactivity Y: virtual species, `Y + OH → HO2`
- Dark thermal HONO: `H2O → HONO` at constant rate

**Units:** K for temperature, hPa for pressure, s⁻¹ for J-values, molec/cm³ for concentrations.

## 4. Solver Configuration (current)

- **Method:** `Radau` (switched from BDF — see Fix 6 below)
- **Tolerances:** `rtol=1e-3`, `atol=1e4`, `max_step=60 s`
- **Zero-clipping:** removed from inside `system_ode`; applied only at CSV export via `np.maximum(self.sol.y, 0.0)`

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
    -   **Chamber Initialization**: Controls for initial conditions, burst injections, custom reactions, and excluded reactions.
    -   **Discovered Parameters**: Searchable tree of all parameters found in the loaded chemistry files, with UNSAFE parameters automatically flagged.
-   **Active Scenario Summary**: A read-only text area that details the active state of the matrix, including any active reaction exclusions.
-   **Scientific Help**: Provides tooltips with scientific explanations for various parameters.
-   **JSON Output**: The `save_json` function serializes the configuration into a structured JSON file.

### UNSAFE Parameter Classification (2026-03-28)
The **Discovered Parameters** tab auto-classifies parameters from the loaded EASY files into safe and unsafe categories. Unsafe parameters are automatically unchecked and coloured red with a warning tooltip. They are never written into `active_parameters` in the JSON.

**UNSAFE (always blocked):** `KMT01–04`, `KMT07–10`, `KMT12–14`, `KFPAN`, `KBPAN` — COMPRC-based Troe rates handled by `mcm_kinetics.py`. `COMPRC`, `PHTFKT`, `PHTFKT_CS`, `PHOTMMSC`, `DENSITY`, `ARRH`, `AMAX`, `AMIN` — EASY `FUNCTION[...]` IDL templates that cannot be evaluated as Python callables.

**Why:** If these were activated they would overwrite the correct Python implementations in the ODE context with `0.0` (SafeDict fallback), silently disabling large parts of the radical chain chemistry.

### Custom Reactions
The "Custom Reactions" table in the Chamber Initialization tab allows users to add ad-hoc reactions with unknown VOC species. Each row specifies Reactants, Products, a rate expression `k` (valid Python using `T`, `M`, etc.), and a comment. Species not found in the loaded mechanism files are highlighted orange. Active rows are written to `custom_reactions` in the JSON and added to `self.mechanism.reactions` in the runner before solving.

### Excluded Reactions (2026-03-28)
The "Excluded Reactions" table in the Chamber Initialization tab allows users to suppress specific reactions that exist in the loaded chemistry files. Workflow:
1. Type the Reactants and Products of the reaction to remove (e.g. `C2H5OH + OH` / `HOCH2CH2O2`).
2. The GUI validates the entry in real time against all reactions parsed from the loaded files and shows `✓ Found (N)` or `✗ Not found`.
3. Matching is stoichiometry-coefficient-agnostic and order-insensitive (uses frozensets of species names).
4. Active exclusions are written to `excluded_reactions` in the JSON.
5. `scenario_runner.py` filters out matching reactions in `_build_mechanism` **before** artifacts are added, printing how many reactions were removed.

## 5. JSON Schema (`scenario_config.json`)

Key fields and their consumers:

| Field | Type | Consumer |
|---|---|---|
| `mechanism_files` | list of filenames | `_build_mechanism` — sorted `.rate` first, then `.mech`/`.incl` |
| `env_map` | dict of var→{type,value} | `_setup_environment` — builds interpolators or constants |
| `j_scalings` | dict of J-name→expression | appended to `var_assignments` after file parse |
| `artifacts` | dict of flags+equations | `_add_artifacts` — adds HONO/HCHO/Y/dark-HONO reactions |
| `initials` | dict of species→ppb | sets `y0` before solve |
| `continuous_sources` | list of {target,rate} | `dydt += eval(rate, ctx)` each ODE step |
| `burst_injections` | list of {target,ppb,start_utc,duration_s} | Gaussian pulse injected into `dydt` |
| `custom_reactions` | list of {reactants,products,rate} | appended to `self.mechanism.reactions` |
| `excluded_reactions` | list of {reactants,products} | filtered out of `self.mechanism.reactions` before artifacts |
| `active_parameters` | list of parameter names | **DEPRECATED — no longer consumed by the runner.** Field may still exist in old JSON files but is ignored. All uncommented parameters from EASY files are now loaded unconditionally. |

## 6. Bug Fixes Applied (2026-03-28)

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

### Fix 6 — CONST() leading-whitespace crash (`TypeError: unsupported operand type(s) for *: 'function' and 'float'`)
**Root cause (three-step failure chain):**
1. EASY rate files write `CONST( expr )` with a space after the opening `(`, e.g. `KMT05 = CONST( 1.44D-13*(1+(M/4.2D+19)) )`. The old CONST-stripping regex `\((.*)\)` captured ` expr ` with a leading space.
2. `compile(' expr', '<string>', 'eval')` raises `IndentationError: unexpected indent` — Python's `eval` mode does not allow leading whitespace. The parameter silently failed to compile and was **not** added to `var_assignments`. Affected: `KMT05`, `KMT06_PARAM`, `KMT11_PARAM`, `RO2`, `RO`, `PANS`.
3. Because `var_assignments` never ran for KMT05, the `lambda *args, f=func: f(*args, T_now, M_now)` loaded from `mcm_kinetics` remained in `ctx['KMT05']`. Any reaction using `KMT05` as its rate evaluated to the lambda → `v = k = <function>` → `v *= y[idx]` → `TypeError` crash.

**Why diagnostics didn't catch it:** The pre-flight test context uses `SafeDict`, which returns `0.0` for unknown names. In that context `ctx['KMT05']` = 0.0, so `eval('KMT05') = 0.0` — reaction appears dead but does not crash. In the real ODE, mcm_kinetics has already populated `ctx['KMT05']` with a live lambda.

**Three fixes applied:**
1. **`easy_parser.py` line 69** — regex changed to `\(\s*(.*?)\s*\)\s*$` so leading/trailing whitespace inside `CONST(...)` is stripped at parse time.
2. **`scenario_runner.py` `_build_mechanism`** — added `p_clean = p_clean.strip()` before `compile()` as defense-in-depth.
3. **`scenario_runner.py` ODE reaction loop** — added `if callable(k): k = 0.0` after `eval(rxn.safe_expr)`. A live lambda in `ctx` now falls back to zero rate instead of crashing the solver.

### Phase 4 Stiffness & Radical Spikes (RESOLVED 2026-03-28)
**Problem:** When running the full ethanol chemistry (Phase 4), highly reactive radical species (like OH, HO2, and RO2) exhibited artificial "spikes" or jagged chatter in their timeseries, despite the overall macroscopic trends (e.g., O3, NO2, VOC decay) looking correct.

**Root cause diagnosed:** The `np.maximum(y, 0.0)` zero-clipping inside `system_ode` was **poisoning the BDF Jacobian**. BDF computes its Jacobian via finite differences: `J_ij ≈ (f(y + ε·eⱼ) - f(y)) / ε`. When a radical went slightly negative (normal internal predictor arithmetic), both the base and perturbed evaluations were clipped to 0, making `J_OH_* ≈ 0`. The solver went blind to radical dependencies and produced chaotic spikes.

**Fix applied (two changes to `scenario_runner.py`):**
1. **Removed `y = np.maximum(y, 0.0)` from inside `system_ode`** — the solver now works with its own internal (possibly transiently negative) values, preserving Jacobian integrity.
2. **Switched solver from `BDF` to `Radau`** — Radau is a fully implicit Runge-Kutta method, more robust than BDF for extreme stiffness ratios (radicals at µs timescales vs. VOCs at hour timescales). Slower per step but far fewer rejected steps.
3. **Applied clipping at export only** — `np.maximum(self.sol.y, 0.0)` now runs in `_export_results` so saved CSV values are never negative, but the solver's internal arithmetic is unaffected.

**Result:** Clean radical timeseries with correct OH, HO2, RO2 dynamics. PSSA implementation is no longer required.

## 7. ODE Evaluation Order

```
_build_mechanism():
  1. Sort mechanism_files: all .rate files first (in JSON order), then .mech/.incl (in JSON order)
  2. Parse each file → accumulate species, reactions, parameters into self.mechanism
  3. ALL uncommented parameters → var_assignments (no whitelist; parser already blocks RESERVED/invalid names)
  4. Append j_scalings to var_assignments
  5. Add custom_reactions from JSON
  6. Filter excluded_reactions (frozenset matching, coefficient-agnostic)
  7. _add_artifacts() → HONO/HCHO/Y/dark-HONO reactions appended

system_ode(t, y):
  1. Load env_map data (TEMP, PRESS, JNO2, JO1D, JHONO, FLOW, HV, SZ, etc.)
  2. Water auto-balancer (RH → H2O)
  3. ctx.update: T, M, O2, N2, math functions
  4. Load mcm_kinetics: COMPRC (direct callable), KMT01-KMT16 as lambdas (overwritten in pass 1)
  5. Load species concentrations into ctx
  --- PASS 1 var_assignments ---
  6. Evaluates COSX, LNO2/MNO2/NNO2, JSCALE, KMT01-KMT14 (via COMPRC), KFPAN/KBPAN, etc.
     (PHTFKT-based J-values fail here — expected, PHTFKT not yet callable)
  --- Build PHTFKT Python callable ---
  7. Creates ctx['PHTFKT'], ctx['PHTFKT_CS'], ctx['AMAX'], ctx['AMIN']
  --- PASS 2 var_assignments ---
  8. Evaluates J1-J51 (via PHTFKT), RO2 sum, etc.
  --- Auto-linker (fallback) ---
  9. Seeds J1←JO1D, J4←JNO2, J7←JHONO, J3←JH2O2, J11←JHCHO_R, J12←JHCHO_M if still 0
  --- Reaction loop ---
  10. For each reaction: eval rate → if callable(k): k=0.0 → accumulate dydt
  11. Apply global dilution: dydt -= y * k_dil
  12. Apply continuous sources and burst injections
```

## 8. Known Gotchas for Future Development

- **`active_parameters` is gone from the runner.** Parameters from EASY files are now loaded unconditionally (all uncommented entries whose names pass the parser's `^[A-Z][A-Z0-9_]*$` filter). The GUI's UNSAFE classification still prevents COMPRC-based KMTs and FUNCTION[] templates from being written into the JSON, but even if they were, the parser would not emit them as `parameters` entries because they appear in `RESERVED`. Do not reinstate a runner-level whitelist — it was the root cause of Fix 7's regression.
- **"Last definition wins" for parameters.** `var_assignments` is a list, not a dict. If the same parameter is defined in multiple files, it gets evaluated multiple times per ODE step; the last evaluation overwrites `ctx`. File order in `mechanism_files` (after the `.rate`-first sort) determines precedence.
- **Reactions never deduplicated.** If the same reaction appears in two loaded files it runs twice. Avoid loading overlapping chemistry files.
- **CONST() spaces.** EASY files frequently write `CONST( expr )` with interior spaces. The parser now strips these (`\(\s*(.*?)\s*\)\s*$`). If adding new EASY-file parsing logic, always `.strip()` extracted expressions before `compile()`.
- **Diagnostics gap.** The pre-flight test uses `SafeDict` (returns 0.0 for unknowns), so a live lambda in `ctx` appears as 0.0 and passes the check. The actual ODE has mcm_kinetics loaded and would crash. The `callable(k)` guard in the reaction loop is the runtime safety net.

### Fix 7 — `active_parameters` whitelist silently excluded valid rate constants (2026-03-29)

**Root cause:** `scenario_runner._build_mechanism` applied an `active_parameters` whitelist: only parameters whose names appeared in `cfg['active_parameters']` were appended to `var_assignments`. Any parameter not in the list — including legitimately important rate constants parsed from EASY files — was silently dropped. This caused C2H5OH to rise on injection then immediately collapse to zero because key rate constants governing its oxidation pathway were missing from the ODE context.

**Fix:** Removed the `active_parameters` filter loop entirely from `_build_mechanism`. All uncommented parameters whose names pass the parser's `^[A-Z][A-Z0-9_]*$` guard now flow into `var_assignments` unconditionally. The parser already excludes RESERVED words (`COMPRC`, `PHTFKT`, `DENSITY`, etc.) and malformed tokens, so the runner-level filter was redundant and harmful.

**Result (2026-03-29):** Simulation confirmed working. C2H5OH rises on injection and decays over hours as expected. Diagnostics report zero failed reactions. KMT01–KMT14 evaluate without warnings.

## 9. Current State (as of 2026-03-29)

All seven bug fixes are applied and the simulation is working correctly:
- Phase 4 (full ethanol oxidation in SAPHIR): confirmed working
- C2H5OH injection, dilution, and OH-driven oxidation all physically sensible
- Radical timeseries (OH, HO2, RO2) smooth — no spikes
- All KMT01–KMT14 Troe rates active and evaluating correctly
- All photolysis J-values (J1–J51) evaluating via PHTFKT Python closure
- Zero failed reactions in `[DIAGNOSTICS]` output

## 10. Next Steps & Future Enhancements

1.  **GUI update:** Cosmetic/structural improvements to `saphir_gui.py` planned now that chemistry is stable.
2.  **Batch Processing Pipeline:** CLI or GUI "Batch Mode" to run multiple `scenario_config.json` files sequentially.
3.  **Data Export & Reporting:** PDF/PNG export buttons in the interactive Matplotlib viewer.
4.  **Deprecate `mechanism_solver.py`:** Its `calculate` function has been fully migrated into `scenario_runner.py`. Legacy stubs can be safely removed.
