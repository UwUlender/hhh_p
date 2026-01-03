# Summary: PETSc Configuration Errors - RESOLVED

## Issue Description

The HydroPlas program was crashing with the following PETSc errors when executed:

```
[0]PETSC ERROR: Unknown type. Check for miss-spelling or missing package
[0]PETSC ERROR: Unable to find requested KSP type GMRES
[0]PETSC ERROR: Unable to find requested PC type PBP
```

Followed by a segmentation fault (SIGSEGV).

## Root Cause Analysis

The errors occurred because:

1. **Invalid KSP Type:** The configuration used `"GMRES"` (uppercase), but PETSc type names are **case-sensitive** and must be lowercase: `"gmres"`

2. **Invalid PC Type:** The configuration used `"PBP"` (Physics-Based Preconditioner), which is **not a valid PETSc preconditioner type**. Valid types include `pbjacobi`, `bjacobi`, `asm`, `ilu`, etc.

3. **Segmentation Fault:** After failing to set the KSP and PC types, PETSc was in an invalid state, leading to the segmentation fault when trying to solve.

## Files Modified

### Source Code Changes

#### 1. `src/config/ConfigParser.hpp`
- **Line 85-86:** Changed default values
  - `preconditioner = "PBP"` → `preconditioner = "pbjacobi"`
  - `ksp_type = "GMRES"` → `ksp_type = "gmres"`

#### 2. `src/config/ConfigParser.cpp`
- **Lines 217-221:** Updated default values in parser
  - `preconditioner = "PBP"` → `preconditioner = "pbjacobi"`
  - `ksp_type = "GMRES"` → `ksp_type = "gmres"`
  
- **Lines 222-230:** Updated fallback defaults
  - `preconditioner = "PBP"` → `preconditioner = "pbjacobi"`
  - `ksp_type = "GMRES"` → `ksp_type = "gmres"`

### Configuration File Updates

#### 3. `config/default_config.yaml`
- Added explicit `solver` section with correct types:
  ```yaml
  solver:
    type: JFNK
    tolerance: 1.0e-6
    max_iterations: 50
    time_step: 1.0e-12
    end_time: 1.0e-9
    ksp_type: gmres
    preconditioner: pbjacobi
  ```

#### 4. `config/complete_feature_demo.yaml`
- **Line 147-148:** Updated solver configuration
  - `preconditioner: PBP` → `preconditioner: pbjacobi`
  - `ksp_type: GMRES` → `ksp_type: gmres`

### Documentation Created

#### 5. `PETSC_SOLVER_FIX.md`
- Detailed explanation of the problem and solution
- List of valid PETSc solver types
- Rebuild instructions

#### 6. `config/test_petsc_fix.yaml`
- New test configuration file
- Demonstrates correct solver configuration
- Includes helpful comments about PETSc requirements

#### 7. `docs/PETSC_SOLVER_GUIDE.md`
- Comprehensive guide to PETSc solver configuration
- Examples for different use cases
- Troubleshooting section
- Performance optimization tips

## Solution Summary

### Required Changes

**All PETSc type names must be lowercase:**
- ✅ Use `gmres` instead of `GMRES`
- ✅ Use `fgmres` instead of `FGMRES`
- ✅ Use `pbjacobi` instead of `PBP` or `PBJACOBI`

**Use valid PETSc preconditioner types:**
- ✅ `pbjacobi` - Point-block Jacobi (parallel-safe, **default**)
- ✅ `bjacobi` - Block Jacobi (parallel-safe)
- ✅ `asm` - Additive Schwarz Method (parallel-safe)
- ✅ `ilu` - Incomplete LU (serial only)
- ✅ `jacobi` - Simple Jacobi (parallel-safe)
- ✅ `none` - No preconditioning

## Rebuilding Instructions

After these changes, rebuild the project:

```bash
cd /path/to/HydroPlas
rm -rf build
mkdir build
cd build
cmake ..
make -j$(nproc)
```

## Testing the Fix

Run the program with any configuration file:

```bash
# Test with default config
./HydroPlas --config config/default_config.yaml

# Test with the new test config
./HydroPlas --config config/test_petsc_fix.yaml

# Test with complete demo
./HydroPlas --config config/complete_feature_demo.yaml
```

The program should now:
1. ✅ Initialize PETSc solvers without errors
2. ✅ Run without segmentation faults
3. ✅ Produce expected output files

## Valid Solver Configurations

### Recommended for Parallel Runs
```yaml
solver:
  ksp_type: gmres
  preconditioner: pbjacobi  # or 'bjacobi', 'asm'
```

### Recommended for Serial Runs
```yaml
solver:
  ksp_type: gmres
  preconditioner: ilu  # Best serial performance
```

### For Better Convergence
```yaml
solver:
  ksp_type: fgmres
  preconditioner: asm
  tolerance: 1.0e-8
  max_iterations: 100
```

## Verification Checklist

- ✅ All source code defaults changed to valid PETSc types
- ✅ All configuration files updated
- ✅ Documentation created
- ✅ Test configuration provided
- ✅ Comprehensive solver guide written

## Additional Notes

1. **Case Sensitivity:** PETSc is **very strict** about case sensitivity. Always use lowercase for type names.

2. **Parallel vs Serial:** Some preconditioners like `ilu` and `sor` only work in serial mode (1 MPI process). For parallel runs, use `pbjacobi`, `bjacobi`, or `asm`.

3. **Configuration Override:** You can override solver settings at runtime:
   ```bash
   ./HydroPlas --config config.yaml -ksp_type fgmres -pc_type asm
   ```

4. **Monitoring:** Enable solver monitoring to debug convergence:
   ```bash
   ./HydroPlas --config config.yaml -ksp_monitor -snes_monitor
   ```

## References

- PETSc Documentation: https://petsc.org/release/docs/
- KSP Manual: https://petsc.org/release/docs/manual/ksp/
- PC Manual: https://petsc.org/release/docs/manual/pc/

## Status

🟢 **RESOLVED** - All issues fixed, code updated, documentation complete.
