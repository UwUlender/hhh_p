#!/bin/bash
# export PETSC_DIR=/usr/lib/petscdir/petsc3.19/x86_64-linux-gnu-real
cd /workspace/HydroPlas/build
# Use mpirun.openmpi provided by the system
mpirun -n 1 ./HydroPlas -snes_fd
