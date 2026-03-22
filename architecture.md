🏗️ SAPHIR-MCM Tool: Architecture Overview

The tool is designed as a Modular Simulation Framework that replaces the legacy MCM/EASY Fortran-based workflow with a modern Python/ODE-solver ecosystem. It bridges the gap between raw chamber measurement data (.enz files) and the Master Chemical Mechanism (MCM) kinetics.
1. Data-Centric Design Pattern

The architecture follows a Decoupled Workflow:

    Front-End (saphir_gui.py): A PyQt6-based "Scenario Builder." Its sole purpose is to generate a scenario_config.json. This JSON acts as the "Chemical Matrix," containing user-selected mechanism files, coordinate mappings for environmental data, and chamber-specific artifact constants.

    The Matrix (scenario_config.json): The single source of truth. It defines what is "Measured" (dynamic interpolation) vs. what is "Fixed" (constants) and stores the "Scaling Factors" for chamber emissions.

    Back-End Engine (scenario_runner.py): A headless simulation script that consumes the JSON. It uses scipy.integrate.solve_ivp (BDF method) to solve the stiff chemical ODEs.

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

easy_parser.py: Reads MCM chemistry files.

mechanism_solver.py: The ODE calculator (handles T, M, RH, and Artifact Equations).

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

## To-Dos

Based on the review of the codebase, here is a list of tasks that need to be completed:

1.  **Implement missing artifacts in `scenario_runner.py`**:
    *   Implement the "Y Background" artifact.
    *   Implement the "Wall Dark HONO" artifact.
2.  **Implement injection logic in `scenario_runner.py`**:
    *   Implement continuous injections.
    *   Implement burst injections.
3.  **Complete the `mechanism_solver.py`**:
    *   The `calculate_rhs` function is very basic. It should be extended to handle all the features described in the `architecture.md` file.
4.  **Improve the `check_results.py` and `test_pipeline.py` scripts**:
    *   These scripts are very basic. They should be improved to provide a more comprehensive testing and analysis framework.
5.  **Update the GUI to support all features**:
    *   The GUI should be updated to allow the user to configure all the features of the simulation, including all artifacts and injections.
6.  **Create a proper command-line interface**:
    *   A proper command-line interface should be created to allow the user to run simulations without using the GUI.
7.  **Add documentation and examples**:
    *   The codebase should be documented, and examples should be provided to help users get started.