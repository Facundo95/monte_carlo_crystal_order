# Monte Carlo Simulation of Alloy Thermodynamics and Magnetism

## 1. Overview

This project implements a Monte Carlo simulation using the **Metropolis algorithm** to study the magnetic and structural properties (specifically Long-Range Order, LRO) of an alloy system, likely a Cu-Al-Mn alloy, on a superlattice structure under varying magnetic fields ($H$) and temperatures ($T$).

The simulation tracks atomic species/structural states (`red`) and magnetic spins/moments (`magn`) on a large 3D lattice, calculating LRO parameters and total magnetization as the system evolves toward equilibrium.

## 2. Project Structure

The code is organized into modular C++ files for clarity, debugging, and maintainability.

| **File** | **Description**|
| :--- | :--- |
| `main.cpp` | Program entry point. Handles command-line arguments, reads simulation parameters from input.txt, and initiates the SimulationLoop. |
| `simulation.h` / `simulation.cpp` | Contains the core simulation logic, the SimulationParameters struct, the SimulationLoop, and the MonteCarloStep (Metropolis algorithm).|
|`lattice.h` / `lattice.cpp` | , Contais the Lattice class (for managing the grid and neighbor lookups) |
| `rng.h` / `rng.cpp` | Implements the custom Linear Congruential Generator (LCG) used for generating pseudo-random numbers in the Metropolis acceptance test. |
| `file_handler.h` / `file_handler.cpp` | Provides helper functions for managing file streams, including opening the output file and checking for errors.|
| `input_parser.h` / `input_parser.cpp` | Provides helper functions for parsing the input file before passing to the SimulationParameters struct. |
| `Makefile` | Script for automated compilation, linking, cleaning, and documentation generation.|
| `input.txt` | Required input file for defining runtime parameters and file names (see Section 4). |

## 3. Requirements and Compilation

**Requirements**

* A `C++` compiler (e.g., g++) supporting C++17 standards or newer.

* The `make` utility.

* *Doxygen* (optional, for generating documentation).

**Compilation**

The project uses a `Makefile` to automate the compilation process.
  1. Open your terminal in the project directory.
  2. Run the default target:

```
make
```

This compiles all source files, links them, and creates the executable file named mc_simulation.

**Cleaning**

To remove all generated object files (`.o`), the executable, and the documentation, run:

```
make clean
```

## 4. Running the Simulation

The simulation requires two types of input files: a configuration file for parameters and an initial configuration file for the lattice state.

**A. Parameter Input (`input.txt`)**

You must create an `input.txt` file in the project directory using the format `KEY VALUE`. This file dictates the simulation run settings. It's not necessary to name the file as the example, you can use any name and extension you want.

```
# Example content for input.txt

# --- File Information ---
FILENAME_PREFIX cu-al-mn_0.67-0.25-0.08
INITIAL_STATE_FILE cu-al-mn_0.67-0.25-0.08_365.txt

# --- Simulation Parameters ---
STEP_NUMBER 100000 
J_M3 150.0
J_M6 100.0

# --- Temperature (T) Sweep ---
T_UPPER 365.0
T_LOWER 365.0
STEP_T 10.0

# --- Magnetic Field (H) Sweep ---
H_UPPER 200.0
H_LOWER -200.0
STEP_H 1.0

# --- Output Flag (0 = disabled, 1 = enabled) ---
SAVE_FINAL_CONFIG 0
```

**B. Initial Lattice State**

Ensure the file specified by `INITIAL_STATE_FILE` (e.g., `cu-al-mn_0.67-0.25-0.08_365.txt`) is present in the directory.

**C. Execution**

Ensure the `input.txt` and initial state file are ready.

Execute the compiled program:

```
./mc_simulation -in input.txt
```

**D. Output**

The results will be written to a file named using the `FILENAME_PREFIX` and the suffix `_out.txt` (e.g., `cu-al-mn_0.67-0.25-0.08_out.txt`). This file contains the calculated LRO parameters, magnetization, and accumulated energy for the measured steps.

## 5. Generating Documentation

If you have *Doxygen* installed and have created a `Doxyfile` in the project root, you can generate the HTML documentation:

```
make doc
```

The resulting documentation files will be placed in the `./html` directory.
