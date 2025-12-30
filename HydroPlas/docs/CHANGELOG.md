# Changelog

All notable changes to the HydroPlas project.

## [1.0.0] - 2025-12-30

### Major Features Added

#### Explicit Excited Species Transport
- **NEW:** Complete implementation of Advection-Diffusion-Reaction (ADR) equations for excited neutral species
- **NEW:** Scharfetter-Gummel discretization extended to neutral species with gas velocity
- **NEW:** Péclet number calculation for automatic regime detection (diffusion vs. advection)

#### Comprehensive Reaction Chemistry
- **NEW:** `ReactionHandler` class with 9 reaction mechanisms:
  - Direct ionization (e + A → 2e + A⁺)
  - Excitation (e + A → e + A*)
  - Stepwise ionization (e + A* → 2e + A⁺)
  - Penning ionization (A* + A* → A⁺ + A + e)
  - Metastable pooling (memory effects)
  - Superelastic collisions (electron heating)
  - Radiative decay (species-dependent)
  - Collisional quenching
  - Three-body recombination
- **NEW:** Automatic energy accounting for electron energy equation
- **NEW:** Fallback Arrhenius forms when lookup table incomplete

#### Advanced I/O System
- **NEW:** `OutputManager` class with dual-mode output:
  - Text output (ASCII, always available)
  - HDF5/OpenPMD output (optional, visualization-ready)
- **NEW:** OpenPMD-compliant hierarchical structure:
  - `/data/{iteration}/meshes/{field_name}`
  - SI units metadata
  - Time and geometry attributes
- **NEW:** Automatic field naming for arbitrary number of species
- **NEW:** TSMonitor integration for periodic output

#### Boundary Conditions
- **NEW:** Robin-type boundary conditions for excited neutrals:
  - Surface quenching with sticking coefficient γ_k
  - Thermal flux formulation: (γ_k v_th / 4) n*
  - Advection contribution for flowing gas
- **NEW:** Secondary electron emission from excited species:
  - Coupling of neutral flux to electron equation
  - Memory effect implementation for DBDs

#### Configuration System
- **NEW:** Extended JSON schema for excited species:
  - `diffusion_coeff`: Binary diffusion coefficient [m²/s]
  - `mass`: Species mass [kg]
  - `energy_level`: Excitation energy [eV]
  - `wall_quenching_prob`: Surface sticking coefficient
  - `wall_see_prob`: Secondary emission yield
- **NEW:** Gas velocity parameter for advection
- **NEW:** Gas temperature for thermal velocity calculations

#### Example Configurations
- **NEW:** `config/argon_complete.json` - Full Argon chemistry (Ar_m, Ar_r, Ar2*)
- **NEW:** `config/dbd_argon.json` - Dielectric barrier discharge with memory effect
- **NEW:** `config/plasma_jet_argon.json` - APPJ with high Péclet number
- **NEW:** `config/penning_mixture.json` - Ne/Ar mixture demonstrating non-local ionization
- **UPDATED:** `config/excited_test.json` - Single metastable test case

#### Documentation
- **NEW:** `docs/THEORY.md` (15 pages) - Comprehensive theoretical framework:
  - ADR system derivation
  - Reaction mechanisms with rate expressions
  - Scharfetter-Gummel scheme derivation
  - Boundary condition formulations
  - Solver architecture details
  - References to key literature
- **NEW:** `docs/USER_GUIDE.md` (12 pages) - Practical user manual:
  - Quick start guide
  - Configuration file reference
  - Example use cases
  - Troubleshooting
  - Performance optimization
  - Advanced topics
- **NEW:** `IMPLEMENTATION_SUMMARY.md` - Complete implementation documentation
- **MAJOR UPDATE:** `README.md` - Professional presentation with physics overview

### Enhancements

#### Solver Module
- **ENHANCED:** `Solver.cpp` now uses modular `ReactionHandler` instead of hardcoded sources
- **ENHANCED:** AppCtx extended with `reactions` and `output` pointers
- **ENHANCED:** Automatic FieldSplit configuration includes excited species in transport block
- **ENHANCED:** MonitorOutput callback for periodic data writing

#### Build System
- **ENHANCED:** CMakeLists.txt with optional HDF5 detection
- **ENHANCED:** Conditional compilation flag `-DUSE_HDF5`
- **ENHANCED:** New source files integrated into build

### Code Quality
- **IMPROVED:** Modular design separates concerns (physics, numerics, I/O)
- **IMPROVED:** Extensive inline documentation
- **IMPROVED:** Comprehensive error handling in ReactionHandler
- **IMPROVED:** Minimum density enforcement to prevent numerical issues

### Physics Accuracy
- **IMPROVED:** Energy loss/gain properly tracked for all reactions
- **IMPROVED:** Second-order reactions (Penning, pooling) correctly implemented
- **IMPROVED:** Wall boundary conditions physically realistic (thermal + advection flux)
- **IMPROVED:** Secondary emission coupling enables memory effects

### Performance
- **OPTIMIZED:** Vectorized excited species handling
- **OPTIMIZED:** FieldSplit preconditioner leverages block structure
- **OPTIMIZED:** Efficient Bernoulli function with Taylor series fallback

---

## [0.1.0] - Previous Version

### Initial Features
- Basic drift-diffusion equations for electrons and ions
- Scharfetter-Gummel for charged species (electric drift)
- Simple ionization source term
- PETSc TS implicit time integration
- BOLSIG+ interface (structure only)
- Basic boundary conditions
- Minimal documentation

---

## Future Releases (Planned)

### [1.1.0] - Analytic Jacobian
- Hand-coded Jacobian for faster Newton convergence
- Anticipated 10× speed improvement

### [1.2.0] - Multi-Gas Chemistry
- Ne/Ar, He/N₂, air chemistry
- Cross-species reactions
- Impurity tracking

### [1.3.0] - Photoionization
- τ-integral for far UV transport
- Helmholtz solver option
- Radiation trapping

### [2.0.0] - 3D Support
- Full 3D geometry
- Adaptive mesh refinement
- Parallel I/O with MPI-HDF5

---

## Version Numbering

- **Major version (X.0.0):** Breaking API changes, major physics additions
- **Minor version (0.X.0):** New features, backward-compatible
- **Patch version (0.0.X):** Bug fixes, documentation updates

---

## Contributors

- HydroPlas Development Team
- AI Assistant (Documentation & Implementation)

---

**Last Updated:** December 30, 2025
