# Visual Guide to the Segmentation Fault Fix

## The Problem Visualized

### Before Fix: Out-of-Bounds Access

```
Array cells:  [0] [1] [2] ... [nx-2] [nx-1]
                ↓                       ↓
            boundary               boundary

Poisson loop was accessing:
- At i=0:    x[j][-1][phi]  ← OUT OF BOUNDS! 💥
- At i=nx-1: x[j][nx][phi]  ← OUT OF BOUNDS! 💥
```

### After Fix: Boundary Skip

```
Array cells:  [0] [1] [2] ... [nx-2] [nx-1]
                ↓    ↓             ↓      ↓
             skip  process      process  skip

Poisson loop now:
- At i=0:    SKIP (Dirichlet BC) ✓
- At i=1..nx-2: Process (Laplacian) ✓
- At i=nx-1: SKIP (Dirichlet BC) ✓
```

## The Vector Access Fix

### Before: Direct Access (WRONG)

```
X_prev (Global Vector)
         ↓ [Direct access - NO GHOSTS!]
    x_prev[j][i][k]  ← UNDEFINED BEHAVIOR! 💥
```

### After: Proper Local Vector (CORRECT)

```
X_prev (Global Vector)
         ↓ [Global to Local scatter]
Xprev_loc (Local Vector with ghosts)
         ↓ [Get array representation]
    x_prev[j][i][k]  ← SAFE ACCESS! ✓
```

## Initial Conditions Flow

### Before: Hardcoded Values

```
config.yaml              Code
-----------              ----
e: 1e14        ---X--->  x[i] = 1e14 (hardcoded)
Ar+: 1e14      ---X--->  x[i] = 1e14 (hardcoded)
e_energy: 1.5  ---X--->  [IGNORED!] 💥
```

### After: Config Applied

```
config.yaml              Code
-----------              ----
e: 1e14        ------>  x[j][i][e_idx] = 1e14 ✓
Ar+: 1e14      ------>  x[j][i][Ar_idx] = 1e14 ✓
e_energy: 1.5  ------>  x[j][i][eps] = 1.5 * n_e ✓
```

## Memory Layout Comparison

### Problem Area: Poisson at Boundary

```
BEFORE (causing segfault):
┌───────────────────────────────────┐
│ Valid Memory                      │
├─────┬─────┬─────┬─────┬─────┬────┤
│ [0] │ [1] │ ... │[n-2]│[n-1]│     │
└─────┴─────┴─────┴─────┴─────┴─────┘
  ↓                            ↓
  Access i-1=-1 💥         Access i+1=n 💥
  
AFTER (fixed):
┌───────────────────────────────────┐
│ Valid Memory                      │
├─────┬─────┬─────┬─────┬─────┬────┤
│ [0] │ [1] │ ... │[n-2]│[n-1]│     │
└─────┴─────┴─────┴─────┴─────┴─────┘
  ↓     ↓             ↓      ↓
 Skip Process      Process  Skip ✓
```

## Code Structure Comparison

### setup_dofs() - Before vs After

```
BEFORE:                          AFTER:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Create vectors                   Create vectors
  ↓                                ↓
Use 1D array access              Use 3D DMDA array access
  ↓                                ↓
Set hardcoded 1e14               Set density_floor defaults
  ↓                                ↓
[END] ← No config read! 💥       Loop through config.initial_conditions
                                   ↓
                                 Apply each species density ✓
                                   ↓
                                 Apply electron energy ✓
                                   ↓
                                 [END] ← Config applied! ✓
```

## The Three Fixes in Context

```
Main Simulation Loop
     │
     ├─→ initialize()
     │        │
     │        ├─→ setup_dofs() ←─────────── FIX #3: Apply initial conditions
     │        └─→ setup_solver()
     │
     └─→ solve_step(dt, t)
              │
              └─→ SNESSolve()
                       │
                       └─→ FormFunction() ←─ FIX #1: Boundary check
                                │             FIX #2: X_prev local vector
                                │
                                ├─→ Compute Sources
                                ├─→ Compute Fluxes
                                ├─→ Poisson Equation ← FIX #1 applied here
                                └─→ Boundary Conditions
```

## Summary: What Each Fix Does

```
┌─────────────────────────────────────────────────────────────┐
│ FIX #1: Boundary Check (if i==0 || i==nx-1) continue       │
│ ───────────────────────────────────────────────────────────│
│ Prevents: Accessing array elements beyond allocated memory │
│ Impact:   PRIMARY cause of segfault eliminated            │
│ Location: Poisson equation loop in FormFunction()         │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│ FIX #2: Local Vector for X_prev                            │
│ ───────────────────────────────────────────────────────────│
│ Prevents: Undefined behavior from accessing global vector  │
│ Impact:   Ensures correct time derivative calculation     │
│ Location: Beginning and end of FormFunction()             │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│ FIX #3: Apply Initial Conditions from Config               │
│ ───────────────────────────────────────────────────────────│
│ Prevents: Starting with incorrect/uninitialized values    │
│ Impact:   Simulation uses user-specified initial state    │
│ Location: Complete rewrite of setup_dofs()                │
└─────────────────────────────────────────────────────────────┘
```

## Expected Simulation Behavior

### Before Fixes
```
./HydroPlas --config main_test.yaml
  ↓
Reading configuration... ✓
Initializing Grid... ✓
Initializing Chemistry... ✓
Initializing Solver... ✓
Starting Simulation...
  ↓
[0]PETSC ERROR: Segmentation Violation 💥
CRASH!
```

### After Fixes
```
./HydroPlas --config main_test.yaml
  ↓
Reading configuration... ✓
Initializing Grid... ✓
Initializing Chemistry... ✓
Initializing Solver... ✓
Starting Simulation... ✓
  ↓
Step 100, Time 1.84e-09 ✓
Step 200, Time 3.69e-09 ✓
Step 300, Time 5.53e-09 ✓
...
[Simulation continues successfully]
```

---

**Key Takeaway:** All three fixes work together to ensure memory safety and correct initialization, eliminating the segmentation fault.
