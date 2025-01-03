
# Lattice Boltzmann Method (LBM) Simulator

This project implements a Lattice Boltzmann Method (LBM) simulator for fluid dynamics simulations in 1D, 2D, and 3D.

## Features

- Supports 1D, 2D, and 3D LBM simulations
- Implements D2Q9 lattice model for 3D simulations
- Includes collision and streaming steps
- Applies boundary conditions for various geometries
- Utilizes OpenMP for parallel processing in 3D simulations
- Provides functions for initialization, equilibrium distribution calculation, and result output

## Requirements

- C compiler (gcc recommended)
- OpenMP support (for parallel processing in 3D simulations)
- Math library (libm)

## Compilation

You can use the following command to build, compile and run

```bash
make
```

## Usage

After compilation, run the simulator with:

```bash
./lbm
```

The simulation parameters can be adjusted in the `lbm.h` header file or passed as arguments to the main function.

## File Structure

- `./include`: Header files containing LBM parameters and function declarations and Raylib visualization
- `./src: Implementation of LBM functions and Raylib visualization
- `main.c`: Main file to run the simulation (not provided in the given code snippet)
