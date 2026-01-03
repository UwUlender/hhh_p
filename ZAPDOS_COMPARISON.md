# Comparison with Zapdos (MOOSE Framework)

## Overview
Zapdos is a plasma simulation application built on the **MOOSE Framework** (Multiphysics Object-Oriented Simulation Environment). Unlike **HydroPlas**, which is a standalone solver using PETSc directly with a custom Finite Difference/Volume implementation, Zapdos leverages the Finite Element Method (FEM) and the extensive modularity of MOOSE.

## Key Differences & What Zapdos Does "Better"

### 1. Boundary Conditions (Physics Correctness)
*   **Zapdos:** Implements a wide variety of kinetic and sheath boundary conditions (e.g., `HagelaarElectronBC`, `SakiyamaIonAdvectionBC`). These models correctly account for:
    *   Thermal motion of electrons near walls.
    *   Secondary Electron Emission (SEE) with energy dependence.
    *   Surface charge accumulation (dielectrics).
    *   Reflection coefficients.
*   **HydroPlas:** Currently uses simplified hardcoded boundary conditions in `PlasmaSolver.cpp`. It mainly applies Dirichlet conditions for potential and a basic drift-flux approximation for species, which may not be accurate for low-pressure discharges where kinetic effects at the sheath are important.

### 2. Modularity & Architecture
*   **Zapdos:** Highly modular. New physics (e.g., a new reaction type or boundary condition) can be added as a separate class (`Kernel`, `BC`, `Material`) without modifying the core solver.
*   **HydroPlas:** Monolithic `PlasmaSolver` class. Adding new physics requires editing the main residual loop (`FormFunction`), which increases complexity and risk of bugs.

### 3. Material Properties
*   **Zapdos:** Uses a `Material` system to compute properties (mobility, diffusion, permittivity) that can depend on any variable (T_gas, E-field, etc.).
*   **HydroPlas:** Properties are largely handled by `Species` and `LookupTable` classes, which are somewhat rigid.

## Implementation Plan for HydroPlas

To bridge the gap in correctness and flexibility, we can adopt the **Boundary Condition Strategy** from Zapdos. We cannot easily switch to MOOSE/FEM, but we can refactor our solver to use polymorphic Boundary Condition objects.

### Proposed Architecture

1.  **`BoundaryCondition` (Base Class):** An abstract interface for applying BCs.
2.  **`HagelaarBC` (New Class):** Port the logic from Zapdos' `HagelaarElectronBC` to our Finite Volume context.
3.  **`BoundaryManager` Update:** instead of just storing config, it should instantiate and manage a list of `BoundaryCondition` objects for each boundary/species pair.
4.  **`PlasmaSolver` Refactor:** Replace the hardcoded BC loop in `FormFunction` with calls to `boundary_manager->apply_bcs(...)`.

### Example Implementation: Hagelaar BC
The Hagelaar boundary condition (Hagelaar et al., 2000) is standard for fluid plasma models. It expresses the electron flux $\Gamma_e$ to the wall as:

$$ \Gamma_e = \frac{1-r}{1+r} \left( \pm \mu_e E_n n_e + \frac{1}{4} v_{th} n_e \right) $$

where $r$ is reflection, $\mu_e E_n$ is the drift component (only if directed towards wall), and $v_{th}$ is thermal velocity.

HydroPlas currently lacks the sophisticated thermal velocity handling and reflection coefficients found in Zapdos.
