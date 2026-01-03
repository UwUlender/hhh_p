# Quick Fix Reference for Segmentation Fault

## TL;DR - What Was Fixed

The segmentation fault in HydroPlas was caused by **improper vector access** in PETSc DMDA operations. Two critical fixes were applied to `src/solver/PlasmaSolver.cpp`:

1. **X_prev vector access** - Now properly converts global to local before array access
2. **Boundary checks** - Added explicit checks to prevent out-of-bounds array access

## Files Modified

- `/workspace/HydroPlas/src/solver/PlasmaSolver.cpp` - FormFunction (lines 220-583)

## To Rebuild and Test

```bash
# Navigate to HydroPlas directory
cd /workspace/HydroPlas

# Clean previous build
rm -rf build

# Build with cmake
mkdir build && cd build
cmake ..
make -j4

# Test with the configuration that was failing
./HydroPlas --config ../config/main_test.yaml

# Should now run without segfault
```

## Expected Output (Success)

```
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
...
```

## What Changed (Technical)

### Before (BUGGY):
```cpp
// WRONG: Accessing global vector with array interface
PetscScalar ***x_prev;
ierr = DMDAVecGetArrayDOF(dm, ctx->X_prev, &x_prev); // SEGFAULT HERE
```

### After (FIXED):
```cpp
// CORRECT: Convert to local vector first
Vec Xprev_loc;
ierr = DMGetLocalVector(dm, &Xprev_loc);
ierr = DMGlobalToLocalBegin(dm, ctx->X_prev, INSERT_VALUES, Xprev_loc);
ierr = DMGlobalToLocalEnd(dm, ctx->X_prev, INSERT_VALUES, Xprev_loc);

PetscScalar ***x_prev;
ierr = DMDAVecGetArrayDOF(dm, Xprev_loc, &x_prev); // NOW SAFE

// ... use x_prev ...

// Cleanup
ierr = DMDAVecRestoreArrayDOF(dm, Xprev_loc, &x_prev);
ierr = DMRestoreLocalVector(dm, &Xprev_loc);
```

### Flux Computation Fix:
```cpp
// Added boundary checks for all neighbor access
if (i < nx - 1) {
    // Safe to access x[j][i+1]
    double n_R = x[j][i+1][k];
    // ... compute flux_right
}

if (i > 0) {
    // Safe to access x[j][i-1]
    double n_L = x[j][i-1][k];
    // ... compute flux_left
}
```

## Verification Checklist

- [ ] Code compiles without errors
- [ ] Simulation starts without immediate crash
- [ ] No PETSC ERROR messages appear
- [ ] Step output is printed periodically
- [ ] HDF5 output file is created
- [ ] Results are physically reasonable (densities > 0)

## If Still Having Issues

1. **Enable PETSc debugging:**
   ```bash
   ./HydroPlas --config config/main_test.yaml -start_in_debugger
   ```

2. **Check memory with valgrind:**
   ```bash
   valgrind --leak-check=full ./HydroPlas --config config/main_test.yaml
   ```

3. **Reduce problem size for testing:**
   Edit `config/main_test.yaml` and reduce `x_nodes` array to fewer points

4. **Check PETSc installation:**
   ```bash
   mpirun --version
   pkg-config --modversion PETSc
   ```

## Summary

✅ Fixed improper global vector access → now uses local vectors  
✅ Fixed out-of-bounds array access → added boundary checks  
✅ Maintained PETSc best practices → proper resource management  
✅ No API changes → existing configurations still work  

The code should now run successfully on all test cases including `main_test.yaml`.

---
**Date:** January 3, 2026  
**Status:** Fixed and Ready for Testing
