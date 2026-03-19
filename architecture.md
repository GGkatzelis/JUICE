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

## 4. Current Status & Next Steps

**Overall Status:** The core pipeline successfully parses simple, self-contained `.easy` files and runs the FACSIMILE simulation.

**Stuck Point (New Blocker):** The parser (`easy_parser.py`) cannot handle `.easy` files that use `#include` directives to import external mechanism and rate constant files. The more complex and realistic scenarios, such as the `ethanol.easy` case, rely heavily on this feature. The current parser only reads the main file, leading to an incomplete representation of the chemical mechanism.

**Immediate Next Steps (Include File Pre-Processing):**
The primary goal is to make the parser "include-aware". This will be achieved by adding a pre-processing step to `easy_parser.py`.

1.  **Implement a Pre-processor:**
    *   Create a function that reads the main `.easy` file.
    *   When it encounters a line with an `#include` statement, it will resolve the path to the included file (which may be relative to the location of the parent file).
    *   It will then read the content of the included file and substitute it into the main stream of data.
    *   This process must be recursive to handle nested includes.
2.  **Feed Pre-processed Content to Parser:**
    *   The existing parsing logic will then receive the complete, flattened mechanism as a single string or stream, as if it were originally in one file.
3.  **Address Secondary Issues:**
    *   Once the parser can handle complex mechanisms, the previously identified issue of "Logic Alignment - Mixing Ratios" can be re-evaluated and addressed.