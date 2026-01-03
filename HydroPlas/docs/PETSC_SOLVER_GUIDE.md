# PETSc Solver Configuration Guide

## Overview

HydroPlas uses PETSc (Portable, Extensible Toolkit for Scientific Computation) for solving the coupled plasma fluid equations. This guide explains how to properly configure the solver to avoid common errors.

## Common Errors and Solutions

### Error: "Unable to find requested KSP type GMRES"

**Cause:** PETSc type names are case-sensitive and must be lowercase.

**Solution:** Use `gmres` instead of `GMRES` in your configuration.

### Error: "Unable to find requested PC type PBP"

**Cause:** "PBP" is not a valid PETSc preconditioner type.

**Solution:** Use a valid preconditioner type like `pbjacobi`, `bjacobi`, `asm`, etc.

## Solver Configuration Parameters

### Basic Parameters

```yaml
solver:
  type: JFNK                    # Solver algorithm (currently only JFNK supported)
  tolerance: 1.0e-6             # Nonlinear convergence tolerance
  max_iterations: 50            # Maximum nonlinear iterations per timestep
  time_step: 1.0e-12           # Time step size in seconds
  end_time: 1.0e-9             # Simulation end time in seconds
  ksp_type: gmres              # Linear solver type (must be lowercase!)
  preconditioner: pbjacobi     # Preconditioner type (must be lowercase!)
```

## Valid Solver Types

### KSP Types (Linear Solvers)

PETSc provides various Krylov subspace methods. **All type names must be lowercase.**

| Type | Description | When to Use |
|------|-------------|------------|
| `gmres` | Generalized Minimal Residual | **Default**, good for general problems |
| `fgmres` | Flexible GMRES | Better for variable preconditioning |
| `bcgs` | BiConjugate Gradient Stabilized | Alternative to GMRES |
| `cg` | Conjugate Gradient | Only for symmetric positive definite systems |
| `bicg` | BiConjugate Gradient | Simple iterative method |
| `tfqmr` | Transpose-Free QMR | Memory efficient |
| `richardson` | Richardson iteration | Simple fixed-point iteration |

**Recommended:** Start with `gmres` or `fgmres`.

### Preconditioner Types

Preconditioners accelerate convergence. **All type names must be lowercase.**

| Type | Description | Parallel? | Recommended For |
|------|-------------|-----------|----------------|
| `pbjacobi` | Point-block Jacobi | ✅ Yes | **Default**, parallel runs |
| `bjacobi` | Block Jacobi | ✅ Yes | Parallel, better than pbjacobi |
| `asm` | Additive Schwarz Method | ✅ Yes | Parallel, domain decomposition |
| `jacobi` | Simple Jacobi | ✅ Yes | Simplest option |
| `ilu` | Incomplete LU | ❌ Serial only | Single-process runs |
| `sor` | Successive Over-Relaxation | ❌ Serial only | Single-process runs |
| `lu` | Direct LU factorization | ❌ Serial only | Small problems only |
| `none` | No preconditioning | ✅ Yes | Testing/debugging |

**Recommended:** 
- Parallel runs: `pbjacobi` (default), `bjacobi`, or `asm`
- Serial runs: `ilu` (best convergence)

## Configuration Examples

### Example 1: Default Configuration (Robust)

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

### Example 2: Better Convergence (Parallel)

For difficult problems that don't converge with default settings:

```yaml
solver:
  type: JFNK
  tolerance: 1.0e-8              # Tighter tolerance
  max_iterations: 100            # More iterations
  time_step: 5.0e-13            # Smaller timestep
  end_time: 1.0e-9
  ksp_type: fgmres              # More flexible
  preconditioner: asm           # Better preconditioning
```

### Example 3: Fast Iterations (Less Robust)

For quick testing or well-conditioned problems:

```yaml
solver:
  type: JFNK
  tolerance: 1.0e-5             # Looser tolerance
  max_iterations: 30            # Fewer iterations
  time_step: 2.0e-12           # Larger timestep
  end_time: 1.0e-9
  ksp_type: bcgs                # Faster convergence
  preconditioner: pbjacobi
```

### Example 4: Serial Run (Best Preconditioning)

For single-process runs on a workstation:

```yaml
solver:
  type: JFNK
  tolerance: 1.0e-6
  max_iterations: 50
  time_step: 1.0e-12
  end_time: 1.0e-9
  ksp_type: gmres
  preconditioner: ilu           # Best serial preconditioner
```

## Advanced: Runtime PETSc Options

You can override solver settings at runtime using PETSc command-line options:

```bash
# Use different KSP type
./HydroPlas --config config.yaml -ksp_type fgmres

# Use different preconditioner
./HydroPlas --config config.yaml -pc_type asm

# View solver details
./HydroPlas --config config.yaml -ksp_monitor -ksp_view

# Increase Krylov subspace size
./HydroPlas --config config.yaml -ksp_gmres_restart 50

# Set KSP tolerance
./HydroPlas --config config.yaml -ksp_rtol 1e-8
```

## Troubleshooting

### Problem: Solver doesn't converge

**Solutions:**
1. Reduce time step size
2. Increase `max_iterations`
3. Try different KSP type (`fgmres`, `bcgs`)
4. Try better preconditioner (`asm`, `ilu`)
5. Tighten `tolerance`

### Problem: Solver is too slow

**Solutions:**
1. Use better preconditioner (`asm` instead of `pbjacobi`)
2. Increase time step (if stable)
3. Use `fgmres` instead of `gmres`
4. Run in parallel with multiple MPI processes

### Problem: "Unknown type" error

**Solution:** Make sure all PETSc type names are **lowercase**. PETSc type names are case-sensitive!

### Problem: Segmentation fault after solver initialization

**Possible causes:**
1. Using serial-only preconditioner (`ilu`, `sor`) in parallel run
2. Invalid solver configuration
3. Memory issues

**Solutions:**
1. Use parallel-safe preconditioner (`pbjacobi`, `bjacobi`, `asm`)
2. Check all configuration parameters
3. Run with smaller problem size for testing

## Checking Available Solver Types

To see all available solver types in your PETSc installation:

```bash
# List all KSP types
./HydroPlas --config config.yaml -help | grep -A 30 "KSP Type"

# List all PC types
./HydroPlas --config config.yaml -help | grep -A 50 "PC Type"

# View PETSc configuration
./HydroPlas --config config.yaml -help | grep "Configure options"
```

## Performance Tips

1. **For parallel runs:** Use `pbjacobi` or `asm` preconditioner
2. **For serial runs:** Use `ilu` preconditioner for best performance
3. **For memory-limited systems:** Use `bcgs` or `tfqmr` instead of `gmres`
4. **For difficult problems:** Start with smaller timesteps and tighter tolerances
5. **For production runs:** Use `fgmres` with `asm` for best balance

## References

- [PETSc Documentation](https://petsc.org/release/docs/)
- [PETSc KSP Manual](https://petsc.org/release/docs/manual/ksp/)
- [PETSc PC Manual](https://petsc.org/release/docs/manual/pc/)
- [PETSc SNES Manual](https://petsc.org/release/docs/manual/snes/)

## Summary Checklist

- ✅ All PETSc type names must be **lowercase**
- ✅ Use valid KSP type: `gmres`, `fgmres`, `bcgs`, etc.
- ✅ Use valid PC type: `pbjacobi`, `bjacobi`, `asm`, `ilu`, etc.
- ✅ For parallel runs, use parallel-safe preconditioners
- ✅ For serial runs, `ilu` gives best preconditioning
- ✅ Start with default settings, then tune for your problem
