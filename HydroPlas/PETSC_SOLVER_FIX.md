# PETSc Solver Configuration Fix

## Problem Summary

The HydroPlas program was crashing with PETSc errors:
```
[0]PETSC ERROR: Unable to find requested KSP type GMRES
[0]PETSC ERROR: Unable to find requested PC type PBP
```

## Root Cause

The default solver configuration in `ConfigParser.hpp` and `ConfigParser.cpp` was using **invalid PETSc type names**:

1. **`ksp_type = "GMRES"`** - PETSc type names are **case-sensitive** and should be lowercase: `"gmres"`
2. **`preconditioner = "PBP"`** - This is **not a valid PETSc preconditioner type**

## Files Fixed

### 1. `src/config/ConfigParser.hpp` (lines 85-86)

**BEFORE:**
```cpp
std::string preconditioner = "PBP"; // "PBP", etc.
std::string ksp_type = "GMRES"; // "GMRES", etc.
```

**AFTER:**
```cpp
std::string preconditioner = "pbjacobi"; // Valid PETSc PC types: "pbjacobi", "bjacobi", "asm", "ilu", "none"
std::string ksp_type = "gmres"; // Valid PETSc KSP types: "gmres", "cg", "bcgs", "fgmres"
```

### 2. `src/config/ConfigParser.cpp` (lines 217-221)

**BEFORE:**
```cpp
if (root["solver"]["preconditioner"]) config_.solver.preconditioner = root["solver"]["preconditioner"].as<std::string>();
else config_.solver.preconditioner = "PBP";

if (root["solver"]["ksp_type"]) config_.solver.ksp_type = root["solver"]["ksp_type"].as<std::string>();
else config_.solver.ksp_type = "GMRES";
```

**AFTER:**
```cpp
if (root["solver"]["preconditioner"]) config_.solver.preconditioner = root["solver"]["preconditioner"].as<std::string>();
else config_.solver.preconditioner = "pbjacobi";

if (root["solver"]["ksp_type"]) config_.solver.ksp_type = root["solver"]["ksp_type"].as<std::string>();
else config_.solver.ksp_type = "gmres";
```

### 3. `src/config/ConfigParser.cpp` (lines 222-230 - else block)

**BEFORE:**
```cpp
} else {
    // Set defaults if solver section is missing
    config_.solver.type = "JFNK";
    config_.solver.tolerance = 1.0e-6;
    config_.solver.max_iterations = 50;
    config_.solver.time_step = 1.0e-12;
    config_.solver.end_time = 1.0e-9;
    config_.solver.preconditioner = "PBP";
    config_.solver.ksp_type = "GMRES";
}
```

**AFTER:**
```cpp
} else {
    // Set defaults if solver section is missing
    config_.solver.type = "JFNK";
    config_.solver.tolerance = 1.0e-6;
    config_.solver.max_iterations = 50;
    config_.solver.time_step = 1.0e-12;
    config_.solver.end_time = 1.0e-9;
    config_.solver.preconditioner = "pbjacobi";
    config_.solver.ksp_type = "gmres";
}
```

## Valid PETSc Solver Types

### KSP Types (Krylov Solvers)
- `gmres` - Generalized Minimal Residual (default, good for general problems)
- `fgmres` - Flexible GMRES (allows variable preconditioning)
- `bcgs` - BiConjugate Gradient Stabilized
- `cg` - Conjugate Gradient (for symmetric positive definite)
- `bicg` - BiConjugate Gradient
- `tfqmr` - Transpose-Free Quasi-Minimal Residual

### PC Types (Preconditioners)
- `pbjacobi` - Point-block Jacobi (default, parallel-safe)
- `bjacobi` - Block Jacobi
- `asm` - Additive Schwarz Method (good for parallel)
- `ilu` - Incomplete LU (serial only)
- `none` - No preconditioning
- `jacobi` - Simple Jacobi
- `sor` - Successive Over-Relaxation (serial only)
- `lu` - Direct LU factorization (for small problems)

## Rebuilding Instructions

After applying these fixes, rebuild the project:

```bash
cd /path/to/HydroPlas
rm -rf build
mkdir build
cd build
cmake ..
make -j$(nproc)
```

## Configuration File Usage

You can now specify solver options in your YAML configuration files:

```yaml
solver:
  type: JFNK
  tolerance: 1.0e-6
  max_iterations: 50
  time_step: 1.0e-12
  end_time: 1.0e-9
  ksp_type: gmres        # Use lowercase!
  preconditioner: pbjacobi  # Use valid PETSc type
```

### Example Solver Configurations

**For better convergence (more robust):**
```yaml
solver:
  ksp_type: fgmres
  preconditioner: asm
  tolerance: 1.0e-8
  max_iterations: 100
```

**For faster iterations (less robust):**
```yaml
solver:
  ksp_type: gmres
  preconditioner: pbjacobi
  tolerance: 1.0e-6
  max_iterations: 50
```

**For serial runs with better preconditioning:**
```yaml
solver:
  ksp_type: gmres
  preconditioner: ilu
  tolerance: 1.0e-6
  max_iterations: 50
```

## Testing the Fix

After rebuilding, test with your configuration:

```bash
./HydroPlas --config config/default_config.yaml
```

The PETSc errors should be resolved and the solver should initialize properly.

## Additional Notes

- PETSc type names are **always lowercase**
- Some preconditioners (like `ilu`, `sor`) only work in **serial mode** (single MPI process)
- For parallel runs, use `pbjacobi`, `bjacobi`, or `asm`
- You can check available types at runtime with: `-help | grep -A 20 "KSP Type"`
