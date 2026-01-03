# Segmentation Fault Fix for HydroPlas

## Problem Description
The segmentation fault was occurring during execution with the following error:
```
[0]PETSC ERROR: Caught signal number 11 SEGV: Segmentation Violation, probably memory access out of range
```

## Root Causes Identified

### 1. **Array Out-of-Bounds Access in Poisson Solver** (CRITICAL)
**Location:** `src/solver/PlasmaSolver.cpp`, lines ~420-435 (in the Poisson equation section)

**Issue:** The code was accessing `x[j][i+1][idx_phi]` and `x[j][i-1][idx_phi]` for ALL cells including boundary cells, causing out-of-bounds access.

**Fix:** Skip boundary cells in the Poisson loop since they use Dirichlet boundary conditions:
```cpp
// Skip boundary cells - they get Dirichlet BC
if (i == 0 || i == nx-1) continue;
```

### 2. **Incorrect X_prev Array Access** (CRITICAL)
**Location:** `src/solver/PlasmaSolver.cpp`, line ~236

**Issue:** The code was trying to access `X_prev` (a global vector) as a 3D array without properly converting it to a local vector with ghost cells.

**Fix:** Create a local vector from X_prev before accessing:
```cpp
// Need local version of X_prev for accessing ghosted values
Vec Xprev_loc;
ierr = DMGetLocalVector(dm, &Xprev_loc); CHKERRQ(ierr);
ierr = DMGlobalToLocalBegin(dm, ctx->X_prev, INSERT_VALUES, Xprev_loc); CHKERRQ(ierr);
ierr = DMGlobalToLocalEnd(dm, ctx->X_prev, INSERT_VALUES, Xprev_loc); CHKERRQ(ierr);
ierr = DMDAVecGetArrayDOF(dm, Xprev_loc, &x_prev); CHKERRQ(ierr);
```

And restore it properly at the end:
```cpp
ierr = DMDAVecRestoreArrayDOF(dm, Xprev_loc, &x_prev); CHKERRQ(ierr);
ierr = DMRestoreLocalVector(dm, &Xprev_loc); CHKERRQ(ierr);
```

### 3. **Initial Conditions Not Applied** (MAJOR)
**Location:** `src/solver/PlasmaSolver.cpp`, `setup_dofs()` function

**Issue:** The configuration file specifies initial conditions, but they were never applied to the solution vector. The code was using hardcoded values instead.

**Fix:** Rewrote the `setup_dofs()` function to:
1. Initialize with density floor values
2. Apply initial conditions from the configuration file for each species
3. Properly handle the electron energy density initialization

## Files Modified

1. **src/solver/PlasmaSolver.cpp**
   - Fixed `setup_dofs()` to apply initial conditions from config
   - Fixed `FormFunction()` to properly handle X_prev as a local vector
   - Fixed Poisson equation loop to skip boundary cells
   - Added proper cleanup of local vectors

## Testing Recommendations

After applying these fixes:

1. **Rebuild the code:**
   ```bash
   cd build
   make clean
   make -j4
   ```

2. **Run with debugging enabled:**
   ```bash
   ./HydroPlas --config config/main_test.yaml -start_in_debugger
   ```

3. **Run with memory checking (if valgrind is available):**
   ```bash
   valgrind --leak-check=full ./HydroPlas --config config/main_test.yaml
   ```

4. **Check initial conditions are correct:**
   - The simulation should start with densities matching the `initial_condition` section in the YAML file
   - Electron density: 1e14 m^-3
   - Electron energy: 1.5 eV
   - Ar+ density: 1e14 m^-3
   - Ar1s5 density: 1e14 m^-3

## Additional Notes

### Energy Equation Fluxes
The current implementation does not compute energy fluxes in the flux calculation section. This is a simplification that may need to be addressed for accurate energy transport. The energy equation currently only includes:
- Time derivative term
- Source term (currently set to 0, needs implementation)

For a more complete implementation, consider adding:
- Conductive energy flux: `-5/2 * D_e * grad(n_eps)`
- Convective energy flux: `5/2 * Gamma_e * mean_energy`

### Boundary Conditions
The boundary condition implementation sets:
- Dirichlet BC for potential at electrodes
- Zero density for ions at walls
- Secondary electron emission (SEE) based on ion flux
- Wall quenching for neutral species

Make sure the `gamma_see` parameter in your config file is appropriate for your simulation.

## Summary of Changes

The segmentation fault was primarily caused by:
1. **Out-of-bounds array access** in the Poisson equation solver at boundary cells
2. **Improper vector access** for the previous time step solution

Both issues have been fixed in the updated source code. The code should now run without segmentation faults.
