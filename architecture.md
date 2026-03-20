# MCM_SAPHIR_Tool Architecture & Status

## 1. Project Overview
The **MCM_SAPHIR_Tool** is a Python-based pipeline designed to translate chemical mechanism definitions (from `.easy` format) into FACSIMILE compatible scripts (`.fac`) and header files (`.hd`). It then automates the execution of the legacy FACSIMILE Fortran executable (`facsim_b.exe`) to run simulations and generate output data (`.cs`, `.kv`, `.pd`).

## 2. System Architecture & Data Flow
1. **Input:** `.easy` file (e.g., `NOx_O3.easy`, `ethanol.easy`) containing initial conditions, reactions, rate constants, and file definitions. These files can contain `#include` directives to import other files.
2. **Parsing:** `easy_parser.py` reads the `.easy` file and recursively processes any `#include` directives. It parses the combined content and populates a central data structure (`mechanism.py`).
3. **Generation:**
   - `fac_generator.py` reads the `Mechanism` data and translates it into FACSIMILE syntax (`.fac` file). This includes expanding user-defined functions and variables.
   - `hd_generator.py` (to be created) will read the `Mechanism` data and build the output column headers (`.hd` file).
4. **Execution:** `main.py` (to be created) will use `subprocess` to feed the generated `.fac` file into `facsim_b.exe`.
5. **Output:** The FACSIMILE executable reads auxiliary data (e.g., `kammer.enz`) and generates time-series outputs (`.cs` for concentrations, `.kv` for reaction rates, `.pd` for production/destruction).

## 3. File Directory & Roles

The project is currently undergoing a reverse engineering process. The original Python files have been deleted, and new files are being created from scratch to ensure a clean and well-structured implementation.

- **`mechanism.py`**: The core data model, re-created from scratch. Contains dataclasses (`Mechanism`, `Reaction`, `Species`) that hold the parsed data in a structured, accessible format.
- **`easy_parser.py`**: The parser, re-created from scratch. Handles `#include` directives and uses Regular Expressions to extract parameters, initial conditions, input variables, and reactions from the `.easy` file and included files. It populates the `Mechanism` object.
- **`fac_generator.py`**: The FACSIMILE script generator, re-created from scratch. It is responsible for translating the `Mechanism` object into a `.fac` file. The current implementation is focused on generating the core logic and handling expression expansion.

## 4. Mimicking Original File Output
- **Goal:** To produce a `.fac` file that is a functionally equivalent translation of the definitions in the `.easy` file and its includes.
- **Corrected Goal:** The new goal is to correctly expand all user-defined EASY functions and variables, while leaving the built-in FACSIMILE functions (like `AMAX`, `EXP`, etc.) as they are.

## 5. Status Update
The project is being reverse-engineered from the ground up. The original Python files were deleted to facilitate a clean start. The following files have been re-created:
- `mechanism.py`: Defines the core data structures.
- `easy_parser.py`: Parses the `.easy` files.
- `fac_generator.py`: Generates the `.fac` files.

The current focus is on implementing the expression expansion logic in `fac_generator.py`. This has proven to be a challenging task. Several approaches have been attempted, including recursive and iterative methods, but they have so far failed to correctly handle the nested and complex nature of the expressions. The main difficulty lies in the correct order of expansion for nested functions and parameters.

The current state is that the `fac_generator.py` script can generate a basic `.fac` file, but the expression expansion is not yet functional. This is the main blocker for the project.

## 6. Include File Reference
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
