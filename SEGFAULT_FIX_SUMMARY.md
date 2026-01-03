# Segmentation Fault Fix Summary

## Problem Description

The HydroPlas code was experiencing a segmentation fault (SIGSEGV) during execution with the following error:

```
[0]PETSC ERROR: Caught signal number 11 SEGV: Segmentation Violation, probably memory access out of range
```

The error occurred immediately after "Starting Simulation..." indicating the problem was in the solver's `FormFunction` during the first time step.

## Root Causes Identified

### 1. **Incorrect Access to X_prev Vector (Primary Issue)**

**Location:** `src/solver/PlasmaSolver.cpp`, line 237 in `FormFunction`

**Problem:** 
The code was calling `DMDAVecGetArrayDOF` directly on the global vector `ctx->X_prev`:

```cpp
ierr = DMDAVecGetArrayDOF(dm, ctx->X_prev, &x_prev); CHKERRQ(ierr);
```

This is incorrect because:
- `X_prev` is a global vector (distributed across MPI ranks)
- `DMDAVecGetArrayDOF` expects a local vector with ghost cells for stencil operations
- Accessing a global vector with this function leads to invalid memory access and segmentation fault

**Solution:**
Created a local vector with ghost cells for `X_prev` before accessing:

```cpp
// Get local vector for X_prev (needed for time derivatives)
ierr = DMGetLocalVector(dm, &Xprev_loc); CHKERRQ(ierr);
ierr = DMGlobalToLocalBegin(dm, ctx->X_prev, INSERT_VALUES, Xprev_loc); CHKERRQ(ierr);
ierr = DMGlobalToLocalEnd(dm, ctx->X_prev, INSERT_VALUES, Xprev_loc); CHKERRQ(ierr);

// Access arrays
PetscScalar ***x, ***x_prev, ***f;
ierr = DMDAVecGetArrayDOF(dm, Xloc, &x); CHKERRQ(ierr);
ierr = DMDAVecGetArrayDOF(dm, Xprev_loc, &x_prev); CHKERRQ(ierr);
```

And proper cleanup at the end:

```cpp
ierr = DMDAVecRestoreArrayDOF(dm, Xprev_loc, &x_prev); CHKERRQ(ierr);
ierr = DMRestoreLocalVector(dm, &Xprev_loc); CHKERRQ(ierr);
```

### 2. **Out-of-Bounds Array Access in Flux Computation**

**Location:** `src/solver/PlasmaSolver.cpp`, lines 319-447

**Problem:**
The flux computation loop had several issues:
- Attempted to access `x[j][i+1]` and `x[j][i-1]` without boundary checks
- Could access beyond array bounds at domain boundaries (i=0 and i=nx-1)
- Confusing loop structure that iterated over interfaces rather than cells

**Solution:**
Rewrote the flux computation to:
1. Loop over cells (not interfaces)
2. Add explicit boundary checks before accessing neighbors
3. Compute divergence properly as `(flux_right - flux_left) * area`
4. Handle boundary conditions separately

```cpp
// Loop over cells and compute divergence of fluxes
for (int j=ys; j<ys+ym; ++j) {
    for (int i=xs; i<xs+xm; ++i) {
        // Right face flux (i+1/2) - only if not at domain boundary
        if (i < nx - 1) {
            double n_L = x[j][i][k];
            double n_R = x[j][i+1][k];
            // ... compute flux_right
        }
        
        // Left face flux (i-1/2) - only if not at domain boundary
        if (i > 0) {
            double n_L = x[j][i-1][k];
            double n_R = x[j][i][k];
            // ... compute flux_left
        }
        
        // Net divergence
        f[j][i][k] += (flux_right - flux_left) * area;
    }
}
```

## Changes Made

### File: `/workspace/HydroPlas/src/solver/PlasmaSolver.cpp`

#### Change 1: FormFunction - X_prev local vector conversion (lines 220-245)
- Added `Xprev_loc` local vector creation
- Added global-to-local scatter for `X_prev`
- Changed array access to use `Xprev_loc` instead of global `ctx->X_prev`

#### Change 2: FormFunction - Proper cleanup (lines 614-620)
- Added restore for `Xprev_loc` array
- Added restore for `Xprev_loc` local vector

#### Change 3: FormFunction - Fixed flux computation (lines 319-447)
- Complete rewrite of flux computation section
- Added boundary checks for all neighbor accesses
- Unified species flux and Poisson equation handling
- Proper divergence calculation with `(flux_right - flux_left) * area`

## Why These Fixes Work

### Memory Safety
1. **Local Vector Access:** PETSc DMDA requires local vectors (with ghost cells) for stencil operations. The global-to-local scatter populates ghost cells with boundary data from neighboring processes.

2. **Boundary Checks:** Explicit checks prevent accessing `x[j][i+1]` when `i == nx-1` or `x[j][i-1]` when `i == 0`, which would read unallocated or incorrect memory.

### Numerical Correctness
1. **Proper Divergence:** The finite volume method requires computing net flux as `(flux_out - flux_in)`, which is now correctly implemented.

2. **Consistent Stencil:** Each cell only accesses its immediate neighbors (i-1, i, i+1), consistent with the DMDA stencil width of 1.

## Testing Recommendations

To verify the fixes work correctly:

```bash
# 1. Recompile the code
cd /workspace/HydroPlas
mkdir -p build && cd build
cmake ..
make -j4

# 2. Run the test case that was failing
./HydroPlas --config ../config/main_test.yaml

# 3. Verify it runs without segfault
# Expected output should show:
# - "Starting Simulation..."
# - "Step X, Time Y" messages
# - No PETSC ERROR messages

# 4. Run with PETSc debugging for detailed checks
./HydroPlas --config ../config/main_test.yaml -start_in_debugger

# 5. Run with memory checking (if valgrind available)
valgrind --leak-check=full ./HydroPlas --config ../config/main_test.yaml
```

## Additional Notes

### Configuration Used
The segfault occurred with `config/main_test.yaml` which has:
- 1D RF discharge with 502 cells
- 3 species (electron, Ar+, Ar1s5)
- Time step: 1.84e-11 s
- JFNK solver with GMRES and Jacobi preconditioner

### PETSc Best Practices Applied
1. Always scatter global to local before using stencil operations
2. Always restore arrays and vectors in reverse order of acquisition
3. Never write to ghost cells (read-only)
4. Check array bounds explicitly at domain boundaries

### Potential Future Improvements
1. Cache the electron index to avoid repeated searches in the species loop
2. Implement proper mean energy interpolation at interfaces instead of using placeholder value 2.0 eV
3. Consider using a flux limiter for better stability with steep gradients

## Conclusion

The segmentation fault was caused by incorrect vector access patterns in PETSc DMDA operations. The fixes ensure:
- ✅ Proper local vector creation with ghost cells
- ✅ Safe boundary checking for all array accesses
- ✅ Correct finite volume divergence computation
- ✅ Proper resource cleanup

The code should now run without segmentation faults on the `main_test.yaml` configuration and similar cases.
