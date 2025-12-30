# Updates for Excited Species Transport

## Overview
This update implements the "explicit transport protocols for excited species" as requested. The simulation now supports an arbitrary number of excited neutral species, solving the Advection-Diffusion-Reaction (ADR) equation for each.

## Key Changes

### 1. Configuration (`src/config/`)
- **`ConfigParser`**: Updated to parse `excited_species`, `reactions`, `gas_velocity`, and `gas_temperature` from the JSON configuration.
- **`excited_species`**: A list of objects defining properties like diffusion coefficient, mass, and wall interaction probabilities.

### 2. Solver Architecture (`src/solver/`)
- **Dynamic DOFs**: The solver now dynamically allocates Degrees of Freedom based on the number of excited species defined in the config.
- **FieldSplit Preconditioning**: The `PCFIELDSPLIT` setup was updated to include excited species in the transport block (Split 0), keeping the Poisson equation separate (Split 1).
- **`FormIFunction`**: 
    - Implemented the **Advection-Diffusion-Reaction** equation for neutral species.
    - Used **Scharfetter-Gummel** discretization for neutral fluxes, using the background `gas_velocity` to define the local Péclet number (drift).
    - Implemented **Robin Boundary Conditions** for neutrals, accounting for wall quenching ($\gamma_k$) and advective outflow.
    - Coupled neutral fluxes to the electron equation via **Secondary Electron Emission (SEE)**.

### 3. Mesh (`src/mesh/`)
- **`MeshGenerator`**: Updated to accept the total number of DOFs dynamically.

## Usage
To use the new features, add an `excited_species` section to your JSON config:

```json
"chemistry": {
    "gas_velocity": 50.0,
    "excited_species": [
        {
            "name": "Ar_m",
            "diffusion_coeff": 1.5e-4,
            "mass": 6.63e-26,
            "wall_quenching_prob": 1.0,
            "wall_see_prob": 0.05
        }
    ]
}
```

## Theory Verification
- **Advection**: The flux term uses $u_{gas}$ as the convective velocity.
- **Diffusion**: $D_k$ is used for diffusive transport.
- **SG Scheme**: The numerical flux is calculated using the Bernoulli function (SG scheme) to ensure stability in advection-dominated regimes ($Pe \gg 1$).
- **Coupling**: The excited species are fully coupled in the implicit system (Newton-Krylov), allowing for consistent time integration.
