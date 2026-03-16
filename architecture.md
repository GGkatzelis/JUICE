# MCM_SAPHIR_Tool Architecture & Status

## 1. Project Overview
The **MCM_SAPHIR_Tool** is a Python-based pipeline designed to translate chemical mechanism definitions (from `.easy` format) into FACSIMILE compatible scripts (`.fac`) and header files (`.hd`). It then automates the execution of the legacy FACSIMILE Fortran executable (`facsim_b.exe`) to run simulations and generate output data (`.cs`, `.kv`, `.pd`).

## 2. System Architecture & Data Flow
1. **Input:** `.easy` file (e.g., `NOx_O3.easy`) containing initial conditions, reactions, rate constants, and file definitions.
2. **Parsing:** `easy_parser.py` reads the `.easy` file and populates a central data structure (`mechanism.py`).
3. **Generation:** 
   - `fac_generator.py` reads the `Mechanism` data and translates it into FACSIMILE syntax (`.fac` file).
   - `hd_generator.py` reads the `Mechanism` data and builds the output column headers (`.hd` file).
4. **Execution:** `main.py` uses `subprocess` to feed the generated `.fac` file into `facsim_b.exe`.
5. **Output:** The FACSIMILE executable reads auxiliary data (e.g., `kammer.enz`) and generates time-series outputs (`.cs` for concentrations, `.kv` for reaction rates, `.pd` for production/destruction).

## 3. File Directory & Roles

- **`mechanism.py`**: The core data model. Contains dataclasses (`Mechanism`, `Reaction`, `Species`) that hold the parsed data in a structured, accessible format.
- **`easy_parser.py`**: The parser. Uses Regular Expressions to extract parameters, initial conditions, input variables, and reactions from the `.easy` file. It populates the `Mechanism` object.
- **`fac_generator.py`**: The FACSIMILE script generator. Handles topological sorting of variables (to prevent dependency errors) and maps the `Mechanism` object into the specific blocks (`COMPILE initial`, `COMPILE equations`, etc.) required by FACSIMILE.
- **`hd_generator.py`**: The Header generator. Creates the `.hd` file which defines the structure and columns for the `.cs`, `.kv`, and `.pd` output files.
- **`main.py`**: The entry point/CLI. Orchestrates the pipeline: Argument parsing -> easy_parser -> fac/hd generators -> file writing -> subprocess execution of `facsim_b.exe`.

## 4. Current Status
**Overall Status:** The translation pipeline is built, functioning, and successfully executes the simulation from end-to-end! IDE default run support is added. File headers and layout structures successfully map 1:1 with the legacy outputs.

**Current Issue (Logic Alignment - Mixing Ratios):**
The generated Python `.fac` file is missing the mixing ratio conversion logic in the `COMPILE parset` block compared to the legacy `.fac` files. 
Specifically, the original tool converts concentrations to mixing ratios, performs density/T/P assignments, and converts back to concentrations before assigning the reaction rate arrays.

**Immediate Next Steps:**
1. Update `fac_generator.py` to inject the mixing ratio conversion `DO` loops (`* convert to mixing ratios;` and `* convert back to conc;`) inside the `_generate_parset_compile_block()` function.
2. Verify the generated output data exactly matches the legacy values after the equations are added.