# MCM_SAPHIR_Tool Architecture & Status

## 1. Project Overview
The **MCM_SAPHIR_Tool** is a Python-based pipeline designed to translate chemical mechanism definitions (from `.easy` format) into FACSIMILE compatible scripts (`.fac`) and header files (`.hd`). It then automates the execution of the legacy FACSIMILE Fortran executable (`facsim_b.exe`) to run simulations and generate output data (`.cs`, `.kv`, `.pd`).

## 2. System Architecture & Data Flow
1. **Input:** `.easy` file (e.g., `NOx_O3.easy`, `ethanol.easy`) containing initial conditions, reactions, rate constants, and file definitions. These files can contain `#include` directives to import other files.
2. **Parsing:** `easy_parser.py` reads the `.easy` file and recursively processes any `#include` directives. It parses the combined content and populates a central data structure (`mechanism.py`).
3. **Generation:**
   - `fac_generator.py` reads the `Mechanism` data and translates it into FACSIMILE syntax (`.fac` file). This includes expanding user-defined functions and variables.
   - `hd_generator.py` reads the `Mechanism` data and builds the output column headers (`.hd` file).
4. **Execution:** `main.py` uses `subprocess` to feed the generated `.fac` file into `facsim_b.exe`.
5. **Output:** The FACSIMILE executable reads auxiliary data (e.g., `kammer.enz`) and generates time-series outputs (`.cs` for concentrations, `.kv` for reaction rates, `.pd` for production/destruction).

## 3. File Directory & Roles

- **`mechanism.py`**: The core data model. Contains dataclasses (`Mechanism`, `Reaction`, `Species`) that hold the parsed data in a structured, accessible format.
- **`easy_parser.py`**: The parser. Handles `#include` directives and uses Regular Expressions to extract parameters, initial conditions, input variables, and reactions from the `.easy` file and included files. It populates the `Mechanism` object.
- **`fac_generator.py`**: The FACSIMILE script generator. Handles topological sorting of variables and maps the `Mechanism` object into the specific blocks (`COMPILE initial`, `COMPILE equations`, etc.) required by FACSIMILE. It contains the `_expand_expr` function for expanding user-defined functions.
- **`hd_generator.py`**: The Header generator. Creates the `.hd` file which defines the structure and columns for the `.cs`, `.kv`, and `.pd` output files.
- **`main.py`**: The entry point/CLI. Orchestrates the pipeline: Argument parsing -> easy_parser -> fac/hd generators -> file writing -> subprocess execution of `facsim_b.exe`.

**4. Mimicking Original File Output:**
- **Goal:** To produce a `.fac` file that is more similar to the original `ethanol.fac` file, which uses pre-set values and simpler assignments.
- **Problem:** The `_expand_expr` function in `fac_generator.py` was expanding all functions and variables, leading to a much larger and more complex output file than the original.
- **Solution:** The `_expand_expr` function was modified to return the expression as-is, without any expansion. This will prevent the expansion of functions and variables, and the generated file will be more similar to the original `ethanol.fac`.


## 5. Include File Reference
This section provides a summary of the include files used in the `ethanol.easy` scenario.

-   **`anorg.mcm331.dark.rate`**: Defines rate constants for inorganic "dark" reactions (reactions that occur without light). It uses a `COMPRC` function for complex rate coefficients.
-   **`anorg.mcm331.dark.mech`**: Lists the inorganic "dark" chemical reactions.
-   **`anorg.mcm331.phot.rate`**: Defines photolysis frequencies (light-dependent reaction rates) for inorganic reactions, using a `PHTFKT` (photolysis function).
-   **`anorg.mcm331.phot.mech`**: Lists the inorganic photolysis reactions.
-   **`ethanol.fac.saphir.incl`**: Contains constants and functions specific to the MCMv3.3.1 ethanol mechanism. This includes definitions for `T` (temperature), `P` (pressure), aliases for chemical species, and photolysis-related functions. It also defines a version of the `PHTFKT` function.
-   **`ethanol.fac.ro2.incl`**: Defines the `RO2` variable as the sum of various peroxy radical species.
-   **`ethanol.fac.mech`**: Contains the core reaction mechanism for the degradation of ethanol.
-   **`ethanol.fac.rate`**: Provides the corresponding rate constants for the reactions listed in `ethanol.fac.mech`.
-   **`j_clear_sky.incl`**: Defines functions for clear-sky photolysis calculations. It includes definitions for `AMAX`, `AMIN`, and `COSX`.
-   **`Cl2_related.incl`**: A large file containing reactions and rate constants for chlorine chemistry, based on a 2015 publication by Xue et al. It defines many species and reactions and appears to be a self-contained chemical mechanism for chlorine-related atmospheric chemistry.
