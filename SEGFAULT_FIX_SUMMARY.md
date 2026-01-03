# Segmentation Fault Fix Summary

## Problem
The HydroPlas simulation was experiencing a segmentation fault (SEGV signal 11) when running with the `main_test.yaml` configuration file. The error occurred immediately after starting the simulation:

```
[0]PETSC ERROR: Caught signal number 11 SEGV: Segmentation Violation, probably memory access out of range
```

## Root Cause
The segmentation fault was caused by **out-of-bounds array access** in the `FormFunction` routine in `PlasmaSolver.cpp`. Specifically:

### Location 1: Poisson Equation Calculation (lines 420-434)
The code was accessing ghost cell indices without proper boundary checks:
```cpp
// Original problematic code:
for (int i=xs; i<xs+xm; ++i) {
    double phi_R = x[j][i+1][idx_phi]; // Out of bounds when i == nx-1
    double phi_L = x[j][i-1][idx_phi]; // Out of bounds when i == 0
    ...
}
```

When `i == 0`, accessing `x[j][i-1][idx_phi]` attempts to read `x[j][-1][idx_phi]`, which is invalid.
When `i == nx-1`, accessing `x[j][i+1][idx_phi]` attempts to read `x[j][nx][idx_phi]`, which is also invalid (nx is the last valid index, nx+1 is beyond the array bounds).

## Solution

### Main Fix: Added Boundary Checks
Modified the Poisson equation calculation to **skip boundary cells** where Dirichlet boundary conditions are applied:

```cpp
// Fixed code:
for (int j=ys; j<ys+ym; ++j) {
    for (int i=xs; i<xs+xm; ++i) {
         // Skip boundary cells - they have Dirichlet BCs
         if (i == 0 || i == nx-1) continue;
         
         // POISSON: -Div(eps Grad phi) = rho
         double phi_C = x[j][i][idx_phi];
         double phi_R = x[j][i+1][idx_phi];  // Now safe: i < nx-1
         double phi_L = x[j][i-1][idx_phi];  // Now safe: i > 0
         ...
    }
}
```

### Secondary Fix: Matrix Assembly
Also fixed an issue with PETSc's matrix-free Jacobian by:
1. Ensuring proper matrix assembly calls (`MatAssemblyBegin/End`)
2. Simplifying the Jacobian setup to use the same matrix for both J and P
3. Adding proper handling for the matrix-free case

## Files Modified
- `/workspace/HydroPlas/src/solver/PlasmaSolver.cpp`:
  - Line ~392-435: Fixed Poisson equation loop to skip boundary cells
  - Line ~65-111: Simplified Jacobian setup
  - Line ~612-656: Added proper matrix assembly in FormJacobian

## Verification
After the fix, the simulation runs successfully:
```bash
$ ./build/HydroPlas --config /workspace/HydroPlas/config/main_test.yaml
Reading configuration from /workspace/HydroPlas/config/main_test.yaml
Initializing Grid...
Initializing Chemistry...
Initializing Boundary...
Initializing Solver...
Initializing Output...
Solver Configuration:
  Type: JFNK
  Time step: 1.84366e-11 s
  End time: 7.37463e-05 s
  Tolerance: 1e-08
  Max iterations: 200
Starting Simulation...
Step 100, Time 1.84366e-09
Step 200, Time 3.68732e-09
Step 300, Time 5.53097e-09
...
```

## Technical Details

### Why Boundary Cells Are Skipped
In the numerical scheme:
1. **Boundary cells (i=0 and i=nx-1)** have **Dirichlet boundary conditions** for the potential (φ)
2. These are handled separately in the boundary condition section of the code (lines 442-601)
3. For these cells, we directly set `f[j][i][idx_phi] = x[j][i][idx_phi] - V` (where V is the electrode voltage)
4. Therefore, the Poisson equation discretization is **only needed for interior cells** (0 < i < nx-1)

### Array Layout
- The DMDA (Distributed Array) structure in PETSc provides:
  - **Owned cells**: `i ∈ [xs, xs+xm)` 
  - **Ghost cells**: One layer of neighboring cells for stencil operations
- However, at physical boundaries (i=0 or i=nx-1), there are **no ghost cells beyond** the domain
- Attempting to access x[j][-1] or x[j][nx] results in undefined behavior and segmentation faults

## Recommendations
1. **Always add bounds checking** when accessing array elements near boundaries
2. **Separate handling** of boundary conditions from interior discretization
3. Use PETSc's `DMDAGetGhostCorners` if ghost cell access is needed
4. Consider using PETSc's built-in boundary condition utilities (`DMAddBoundary`)

## Status
✅ **FIXED** - The segmentation fault has been resolved and the simulation runs successfully.
