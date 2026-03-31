# SAPHIR-MCM Tool: Architecture Overview

The tool is a full-stack simulation framework replacing the legacy MCM/EASY Fortran-based workflow with a modern Python/ODE-solver ecosystem. It bridges raw chamber measurement data (`.enz` files) and the Master Chemical Mechanism (MCM) kinetics, with GUI-driven scenario building, stiff ODE solving, and interactive visualization.

## 1. Data-Centric Design Pattern

Decoupled three-stage workflow:

- **Front-End** (`saphir_gui.py`): PyQt6 "Scenario Builder." Sole purpose is to produce `scenario_config.json`.
- **The Matrix** (`scenario_config.json`): Single source of truth. Defines Measured (interpolated) vs. Fixed variables, scaling factors, custom/excluded/blocked reactions and parameters.
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
| `results_viewer.py` | Interactive results viewer + process analysis engine (launched by runner) |
| `Python_vs_FACSIMILE.py` | Standalone two-tab comparison GUI: timeseries overlay + slope matrix vs FACSIMILE `.cs` output |

## 3. Core Physics

**Air density M** (Ideal Gas Law): `M = (P×100) / (kB×T) × 1e-6` [molec/cm³]

**Dilution** (first-order decay for all species): `d[C]/dt -= [C] × (FLOW + 0.6) / (V × 3600)`

**SAPHIR artifacts** (modelled explicitly):
- Photolytic HONO wall source: `4.7e13 × (1 + (RH/11.6)²) × exp(−3950/T) × JNO2 × HV`
- Photolytic HCHO wall source: `2.8e6 × HV`
- Missing OH reactivity Y: virtual species, `Y + OH → HO2`
- Dark thermal HONO: zeroth-order `→ HONO` at rate `1e-14 × H2O` (molec cm⁻³ s⁻¹). The coefficient is `1e-14` not `1e-12` because FACSIMILE's formula buries a `1e-2` factor inside its H2O definition (`H2O_FACSIMILE = 1e-2 × M × VMR`), while Python's ODE context has `H2O = VMR × M` — 100× larger. Using `1e-12` would produce a dark source ~100× too strong.

**Units:** K for temperature, hPa for pressure, s⁻¹ for J-values, molec/cm³ for concentrations.

## 4. Solver Configuration (current)

- **Method:** `Radau` (switched from BDF — see Fix 6 below)
- **Tolerances:** `rtol=1e-3`, `atol=1e4`, `max_step=60 s`
- **Zero-clipping:** removed from inside `system_ode`; applied only at CSV export via `np.maximum(self.sol.y, 0.0)`

## SAPHIR GUI (`saphir_gui.py`)

PyQt6 scenario builder. Produces a named `<filename>.json` and can directly launch the runner.

### Workflow
1. **Load Folder** — select scenario root containing `chemistry/` subfolder and a `.enz` file.
2. **Configure** — four tabs (see below).
3. **Active Scenario Summary** — sidebar console live-updates with all active settings and warnings.
4. **Generate JSON + Run** — bottom bar: editable filename field, 💾 save button, ▶ Run button.

### Tabs

**Constrained Variables**
- Maps model variables (TEMP, PRESS, JNO2, …) to `.enz` columns or fixed constants.
- **J-Scalings table** (4 columns: Act / Target / Expression / Source): defines additional J-values derived from measured ones. Source column shows `"Manual"` or the chemistry filename the expression was auto-detected from. **"🔄 Sync from Files"** button scans `master_parameters` for J-named entries not already listed (skips J1–J51, env_map entries, JSCALE) and adds them with source attribution.

**Chamber Sources**
- Artifact checkboxes: HONO light source (Rohrer et al.), HCHO light source, NO2 light source, Y background, wall dark HONO.
- HONO Scaling Factor default: `"1"`.
- Continuous Sources table for background influx rates.

**Chamber Initialization**
- Initial conditions, burst injections (auto-detected from `q*` columns in `.enz`), custom reactions.
- **Custom Reactions**: species not in loaded mechanism highlighted orange. Written to `custom_reactions` in JSON; appended to `self.mechanism.reactions` in runner before solve.
- **Excluded Reactions**: inline autocomplete via `ReactionCompleterDelegate` (`QStyledItemDelegate`). As you type a species name (contains-matching, case-insensitive), the reactants column suggests all matching reactant-set strings from loaded chemistry files. Selecting reactants auto-fills products if only one product set exists. Written to `excluded_reactions`; filtered from `self.mechanism.reactions` in `_build_mechanism` before artifacts, using frozenset (coefficient-agnostic, order-insensitive) matching.

**Chemistry Overview** (previously "Discovered Parameters" — redesigned 2026-03-29)
- Read-only audit tree of all parameters parsed from loaded chemistry files.
- 4 columns: Parameter / Expression / Source File / **Status**.
- Status values: `🟢 Active (auto-loaded)`, `⚠ Backend-managed` (UNSAFE — KMT01-14, COMPRC, PHTFKT etc.), `⚫ Commented — skipped`. No checkboxes — all loading is unconditional (Fix 7).
- **Blocked Parameters** section: user can add specific parameter names to suppress. Written to `blocked_parameters` in JSON; runner strips those names from `var_assignments` after file parse, before j_scalings.

### UNSAFE parameters (informational only)
`KMT01–04`, `KMT07–10`, `KMT12–14`, `KFPAN`, `KBPAN`, `COMPRC`, `PHTFKT`, `PHTFKT_CS`, `PHOTMMSC`, `DENSITY`, `ARRH`, `AMAX`, `AMIN` — shown in orange in the Chemistry Overview tree. These are either handled by the Python backend directly or use EASY `FUNCTION[...]` syntax that fails to compile. Since Fix 7 loads all parameters unconditionally, they are listed as informational warnings, not blocking controls.

### Run Simulation button
"▶ Run Simulation" (bottom right, green) launches `scenario_runner.py <config.json>` as a subprocess via `subprocess.Popen`. If a JSON was just saved, offers "Run last saved" or "Browse for JSON". The GUI remains responsive while the runner window is open.

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
| `active_parameters` | list of parameter names | **DEPRECATED — no longer written or consumed.** Old JSON files may still contain this field; it is silently ignored by the runner. |
| `blocked_parameters` | list of parameter names | `_build_mechanism` — strips named parameters from `var_assignments` after file parse, before j_scalings. Safe blocklist (inverse of the old whitelist). |

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
  4. Apply blocked_parameters: strip named params from var_assignments (user-controlled blocklist)
  5. Append j_scalings to var_assignments (these override file-defined values — last-write-wins)
  6. Add custom_reactions from JSON
  7. Filter excluded_reactions (frozenset matching, coefficient-agnostic)
  8. _add_artifacts() → HONO/HCHO/Y/dark-HONO reactions appended

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

- **`active_parameters` is gone.** Parameters load unconditionally (parser's `^[A-Z][A-Z0-9_]*$` filter is the only gate). Do not reinstate a runner-level whitelist — it was the root cause of Fix 7's regression. Use `blocked_parameters` for targeted suppression instead.
- **`blocked_parameters` is a blocklist (safe).** Specific names can be suppressed via the Chemistry Overview tab. Applied after file parse, before j_scalings. Empty list = no effect on existing simulations.
- **"Last definition wins" for parameters.** `var_assignments` is a list, not a dict. If the same parameter is defined in multiple files, it gets evaluated multiple times per ODE step; the last evaluation overwrites `ctx`. File order in `mechanism_files` (after the `.rate`-first sort) determines precedence.
- **Reactions never deduplicated.** If the same reaction appears in two loaded files it runs twice. Avoid loading overlapping chemistry files.
- **CONST() spaces.** EASY files frequently write `CONST( expr )` with interior spaces. The parser now strips these (`\(\s*(.*?)\s*\)\s*$`). If adding new EASY-file parsing logic, always `.strip()` extracted expressions before `compile()`.
- **Diagnostics gap.** The pre-flight test uses `SafeDict` (returns 0.0 for unknowns), so a live lambda in `ctx` appears as 0.0 and passes the check. The actual ODE has mcm_kinetics loaded and would crash. The `callable(k)` guard in the reaction loop is the runtime safety net.

### Fix 7 — `active_parameters` whitelist silently excluded valid rate constants (2026-03-29)

**Root cause:** `scenario_runner._build_mechanism` applied an `active_parameters` whitelist: only parameters whose names appeared in `cfg['active_parameters']` were appended to `var_assignments`. Any parameter not in the list — including legitimately important rate constants parsed from EASY files — was silently dropped. This caused C2H5OH to rise on injection then immediately collapse to zero because key rate constants governing its oxidation pathway were missing from the ODE context.

**Fix:** Removed the `active_parameters` filter loop entirely from `_build_mechanism`. All uncommented parameters whose names pass the parser's `^[A-Z][A-Z0-9_]*$` guard now flow into `var_assignments` unconditionally. The parser already excludes RESERVED words (`COMPRC`, `PHTFKT`, `DENSITY`, etc.) and malformed tokens, so the runner-level filter was redundant and harmful.

**Result (2026-03-29):** Simulation confirmed working. C2H5OH rises on injection and decays over hours as expected. Diagnostics report zero failed reactions. KMT01–KMT14 evaluate without warnings.

## 9. Current State (as of 2026-03-31)

All seven bug fixes are applied and the simulation is working correctly:
- Phase 4 (full ethanol oxidation in SAPHIR): confirmed working
- C2H5OH injection, dilution, and OH-driven oxidation all physically sensible
- Radical timeseries (OH, HO2, RO2) smooth — no spikes
- All KMT01–KMT14 Troe rates active and evaluating correctly
- All photolysis J-values (J1–J51) evaluating via PHTFKT Python closure
- Zero failed reactions in `[DIAGNOSTICS]` output

**New tools added (2026-03-31):**
- `Python_vs_FACSIMILE.py`: standalone comparison GUI (timeseries + slope matrix)
- Franz's FACSIMILE reference output in `MCM_tested_scenarios/ethanol_scenario/MCM output/case1,case2/`
- `results_viewer.py` (renamed from `check_results.py`) massively expanded:
  - `BudgetEngine`: vectorized process analysis; mirrors runner's var_assignments/artifact logic; H2O unit auto-detection; O1D PSS reconstruction; two-pass vectorized eval; minor pathway grouping
  - `ResultsViewer`: three-folder QTreeWidget left panel (MCM Output / SAPHIR Environment / SAPHIR Instruments); single-click species → budget plot; auto-loads scenario_config.json from CSV path; stacked production/destruction chart

**Reactivity Analysis (added 2026-03-31):**
- `BudgetEngine._prepare_ctx(df)`: extracted shared helper — builds vec ctx, runs two-pass var_assignments, photolysis auto-linker, and O1D PSS. Called by both `compute` and `compute_reactivity`.
- `BudgetEngine.compute_reactivity(radical, df)` → `{label: np.array}` (s⁻¹): for each reaction consuming `radical`, contribution = k × Π[co-reactants] (one instance of radical excluded from product). Covers OH, O3, NO3 (or any species).
- `ResultsViewer` Tab 2 "Reactivity": dropdown (OH/O3/NO3) + "Compute Reactivity" button → stacked area of contributions + total k (black dashed) + HV shading + time-average bar chart. Same layout as Process Analysis. Fully post-facto — no ODE changes needed.

**SAPHIR NC Instrument Loading (added 2026-04-01):**
- `results_viewer.py` replaces `check_results.py` (import updated in `scenario_runner.py`).
- Left panel is now a 3-folder `QTreeWidget`: **MCM Output** (modelled species), **SAPHIR Environment** (`SAPHIR_*` columns from the ENZ/CSV), **SAPHIR Instruments** (loaded from NC files).
- "Load SAPHIR NC Data" button: infers experiment date from the CSV timestamp (`t0_unix` → UTC → `Z:/IEK8-SAPHIR/data/YYYY/MM/DD`), shows an editable confirmation dialog, then lists all `.nc` files in a checkable selector (`_NcFileSelector`).
- NC time axis: `"Seconds since 2000-01-01 00:00:00 UTC"` → converted via `_NC_EPOCH_OFFSET = 946684800` to Unix, then to hours relative to simulation start. No interpolation — instrument data is plotted on its own native time grid.
- `_MatchDialog`: 5-column table (NC File | NC Variable | Units | MCM Match | Load) with auto-suggested MCM match (exact → SAPHIR_-prefixed fallback). User can edit the match column and choose which variables to load.
- Instrument data stored in `self._inst_data[nc_basename/nc_var]` — NOT injected into `self.df` to keep the budget engine clean.
- Single-click on instrument node → "Instrument data" info panel in Process Analysis tab (no budget available). Double-click → plots on chosen axis via `_plot_instrument`.
- `_auto_measured_koh_from_df` extended to also scan `self._inst_data` for kOH variables.
- Reactivity tab: "Load Measured k_OH (.csv/.txt)" replaced by **"Load k_OH (NC file)"** which opens `Z:/IEK8-SAPHIR/data` (or last used NC dir) and reads `kOH` variable from the selected file.

**Fix 8 — dark_hono H2O unit mismatch (2026-03-31):**
FACSIMILE's dark HONO wall source is `H2O → HONO` with `k = 1e-12`, where FACSIMILE defines `H2O = 1e-2 × M × VMR`. Python's ODE context uses `H2O = VMR × M` (100× larger). The expression `"1e-12 * H2O"` therefore ran 100× too fast, causing ~10× offset in simulated HONO concentration (photolytic source moderated the ratio). Fixed by changing the default expression to `"1e-14 * H2O"` in `saphir_gui.py` (line 717) and `scenario_config.json`. The ODE physics (H2O conversion, kinetics, dilution) are correct — only the coefficient in this user-configurable artifact expression was wrong.

### GUI Enhancements (as of 2026-03-29)

- **Logo Integration:** The FZ-Jülich logo has been added to the top application bar for branding and a more professional appearance.

- **Chamber Initialization Refinements:** The "Chamber Initialization" tab and related logic have been significantly overhauled for better automation and accuracy:
    - **Dynamic Burst Injections:** The "Burst Injections" table is no longer pre-populated with default O3/NO injections. It is now dynamically filled by detecting columns prefixed with `Injection_` in the loaded `.enz` file.
    - **Concentration-from-data:** The ppb value for these dynamic injections is now taken directly from the value in the `Injection_*` column at the point of injection, rather than using a hardcoded default.
    - **H2O Injection Guard:** The GUI now prevents an H2O injection from being added to the "Burst Injections" table if H2O is already configured as a constrained environmental variable, avoiding conflicting definitions.
    - **Automatic Initial Conditions:** The "Initial Conditions" table is now automatically populated from the first row of the loaded `.enz` file when a project folder is loaded, streamlining setup.

- **Derived Photolysis Order:** To improve readability, the default manual photolysis scaling rules (e.g., for JCL2) are now programmatically added to the *end* of the "Derived Photolysis Frequencies" list, after any J-values that have been automatically synchronized from the loaded chemistry files.

**Fix 9 — FACSIMILE comparison files have ENZ column name mismatch (2026-03-31):**
`case1` and `case2` FACSIMILE EASY files read `input(qVCP2)` for VOC injection and `input(QNO)` for NO injection, but the actual `kammer.enz` has columns `Injection_Ethanol` and `Injection_NO` — no match. FACSIMILE's `INITIAL[*]=CONST(1.)` sets all uninitialized species to 1, so `qVOC` and `qNO` stay at 1 throughout the run, creating continuous ethanol (591 ppb/min) and NO (~7 ppb/min) injection for 18 hours. This drives a chemical blow-up at ~hour 5.5 (NO2 → 1.6e11, HONO → 1.6e10). The 11× HONO and 2× OH discrepancy seen in comparisons is entirely from this FACSIMILE bug, not a Python model deficiency. In the pre-blow-up period (hours 0–5.5), HONO agrees to within 3% and the Python model is correct.

**Fix 10 — Reactivity tab: measured k_OH overlay (2026-03-31):**
`ResultsViewer` Reactivity tab now has a "Load Measured k_OH" button that loads a two-column CSV/TXT file (time + k_OH in s⁻¹) and overlays it as a red line on the reactivity timeseries. Also auto-detects columns matching `kOH`/`K_OH`/`SAPHIR_KOH` in `self.df`. Active only when radical=OH.

## 10. Next Steps & Future Enhancements

1.  **Overlay measured k_OH** — Load actual SAPHIR CRM k_OH measurement file and overlay on Reactivity tab. Likely in the SAPHIR archive NC file (variable TBD). Check if archive has `kOH` or `OH_reactivity` variable for 09.07.2024.
2.  **Fix FACSIMILE comparison references** — case1/case2 EASY files need correct column names (`qVCP2`→`Injection_Ethanol`, `QNO`→`Injection_NO`) or a corrected ENZ file with those names. Until fixed, FACSIMILE comparison is not valid beyond t < 5.5h. Confirm with Franz whether a corrected FACSIMILE run can be produced.
3.  **Simple CONST scenario** — [scenario_config - Simple Scenario 1.json](MCM_tested_scenarios/ethanol_scenario/scenario_config%20-%20Simple%20Scenario%201.json) — awaiting colleague confirmation on whether species should be held constant or evolved freely.
4.  **Test excluded/custom reactions end-to-end** — user to verify add/remove reactions produce correct ODE changes.
2.  **Test excluded/custom reactions end-to-end** — user to verify add/remove reactions produce correct ODE changes; report results.
3.  **Batch Processing Pipeline:** CLI or GUI "Batch Mode" to run multiple JSON files sequentially.
4.  **Data Export & Reporting:** PDF/PNG export buttons in the interactive Matplotlib viewer.
5.  **Deprecate `mechanism_solver.py`:** Its `calculate` function has been fully migrated into `scenario_runner.py`. Legacy stubs can be safely removed.
