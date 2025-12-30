# HydroPlas: AI-Driven High-Performance Hydrodynamic Plasma Simulation Framework

HydroPlas is a 1D/2D hydrodynamic plasma simulation code built on top of PETSc. It is designed to model low-temperature collisional plasmas using the Drift-Diffusion-Reaction approximation.

## Features

- **Geometry**: 1D and 2D Cartesian grids.
- **Physics**: 
  - Species continuity equations (Drift-Diffusion).
  - Electron Energy transport.
  - Poisson equation for self-consistent electrostatic field.
- **Numerics**:
  - Scharfetter-Gummel exponential flux scheme for stability.
  - Fully Implicit Time Integration (PETSc TS: BDF/Backward Euler).
- **Chemistry**:
  - Flexible lookup table support.
  - Interface design for BOLSIG+ integration.
- **Configuration**: JSON-based runtime configuration.

## Dependencies

- **PETSc** (3.19+ recommended)
- **MPI**
- **nlohmann/json**
- **CMake** (3.14+)
- **C++17** Compiler

## Building

```bash
mkdir build
cd build
cmake ..
make
```

## Running

Ensure PETSc environment variables are set (if not using system install).

```bash
mpirun -n 1 ./HydroPlas -snes_fd
```

Note: The current solver implementation uses a Matrix-Free Finite Difference Jacobian (`-snes_fd`) which may struggle with convergence for stiff plasma sheaths. For production runs, an analytic Jacobian implementation is recommended.

## Configuration

Configuration is loaded from `config/default_config.json`. You can modify domain size, voltage, time steps, and chemistry settings there.

## Project Structure

- `src/main.cpp`: Entry point.
- `src/config`: Configuration parser.
- `src/mesh`: DMDA mesh generation.
- `src/solver`: PETSc TS solver and residual physics (FluxSG, Poisson).
- `src/boundary`: Boundary condition logic (DC, RF, SEE).
- `src/chemistry`: Transport data management and BOLSIG+ interface.
