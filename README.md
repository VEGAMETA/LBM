
# Lattice Boltzmann Method (LBM) Simulator

This project implements a Lattice Boltzmann Method (LBM) simulator for fluid dynamics simulations in 1D, 2D, and 3D.

## Features

- Supports 1D, 2D, and 3D LBM simulations
- Implements D3Q19 lattice model for 3D simulations
- Includes collision and streaming steps
- Applies boundary conditions for various geometries
- Utilizes OpenMP for parallel processing in 3D simulations
- Provides functions for initialization, equilibrium distribution calculation, and result output

## Requirements

- C compiler (gcc recommended)
- OpenMP support (for parallel processing in 3D simulations)
- Math library (libm)

## Compilation

To compile the project, use the following command:

```bash
gcc -o lbm_simulator main.c lbm.c -lm -fopenmp
```

Replace `main.c` with your main file that includes the LBM simulation loop.

You can use the following command to build, compile and run

```bash
make
```

## Usage

After compilation, run the simulator with:

```bash
./lbm_simulator
```

The simulation parameters can be adjusted in the `lbm.h` header file or passed as arguments to the main function.

## File Structure

- `lbm.h`: Header file containing LBM parameters and function declarations
- `lbm.c`: Implementation of LBM functions
- `main.c`: Main file to run the simulation (not provided in the given code snippet)

## Key Functions

- `create_Xd_array`: Functions to create 1D, 2D, 3D, and 4D arrays
- `free_Xd_array`: Functions to free allocated memory for arrays
- `initialize_Xd`: Functions to initialize the simulation for 1D, 2D, and 3D
- `equilibrium_Xd`: Functions to calculate equilibrium distribution
- `apply_boundary_conditions_Xd`: Functions to apply boundary conditions
- `collide_and_stream_Xd`: Functions to perform collision and streaming steps
- `swap_distributions_Xd`: Functions to swap distribution arrays
- `write_to_file_Xd`: Functions to output simulation results

## Output

The simulation results are written to a file named `density.dat` in the current directory.
