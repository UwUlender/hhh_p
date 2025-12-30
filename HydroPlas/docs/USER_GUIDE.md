# HydroPlas User Guide: Simulating Plasmas with Excited Species Transport

## Quick Start

### 1. Building the Code

```bash
cd HydroPlas
mkdir build && cd build
cmake ..
make -j4
```

**Dependencies:**
- PETSc 3.19+ (with MPI)
- nlohmann/json
- HDF5 (optional, for OpenPMD output)

### 2. Running Your First Simulation

```bash
# Basic run with default config
./HydroPlas config/excited_test.json

# With PETSc options
./HydroPlas config/argon_complete.json -ts_monitor -snes_monitor
```

### 3. Visualizing Results

**Text output:**
```bash
cd output
gnuplot
> plot "ne_000100.txt" using 1:3 with lines
```

**HDF5 output (if compiled with HDF5):**
```bash
# ParaView with OpenPMD plugin
paraview output/hydroplas_*.h5

# Python analysis
import h5py
f = h5py.File('output/hydroplas_000100.h5', 'r')
ne = f['/data/100/meshes/ne'][:]
```

---

## Configuration File Structure

### Domain

```json
"domain": {
    "Lx": 0.01,      // Domain length [m]
    "Nx": 200        // Number of grid points
}
```

### Time Integration

```json
"time": {
    "dt": 1e-11,           // Initial time step [s]
    "t_end": 1e-8,         // End time [s]
    "output_interval": 50  // Write output every N steps
}
```

### Boundary Conditions

```json
"boundary": {
    "voltage_type": "DC",           // "DC", "RF", or "Pulsed"
    "voltage_amplitude": 500.0,     // [V]
    "frequency": 13.56e6,           // [Hz] (for RF)
    "gamma_see": 0.1,               // Secondary emission coeff (ions)
    "dielectric_permittivity": 4.0, // Relative εr
    "dielectric_thickness": 1e-3    // [m]
}
```

### Chemistry

#### Pre-calculated Mode (Recommended)

```json
"chemistry": {
    "mode": "Pre-calculated",
    "transport_table_file": "data/transport.dat",
    "gas_velocity": 50.0,        // [m/s], 0 for static
    "gas_temperature": 300.0,    // [K]
    "excited_species": [ ... ]
}
```

#### Inline BOLSIG+ Mode

```json
"chemistry": {
    "mode": "Inline BOLSIG+",
    "cross_section_file": "data/argon.txt",
    // ... rest same as above
}
```

---

## Excited Species Configuration

### Minimal Example (Single Metastable)

```json
"excited_species": [
    {
        "name": "Ar_m",
        "diffusion_coeff": 1.5e-4,      // [m²/s]
        "mass": 6.63e-26,               // [kg]
        "energy_level": 11.55,          // [eV]
        "wall_quenching_prob": 0.001,   // γ_quench
        "wall_see_prob": 0.01           // γ_see (electrons per quenched atom)
    }
]
```

### Complete Argon System

```json
"excited_species": [
    {
        "name": "Ar_m",
        "diffusion_coeff": 1.5e-4,
        "mass": 6.63e-26,
        "energy_level": 11.55,
        "wall_quenching_prob": 0.001,
        "wall_see_prob": 0.01,
        "comment": "Metastable 1s5 - long lifetime"
    },
    {
        "name": "Ar_r",
        "diffusion_coeff": 1.5e-4,
        "mass": 6.63e-26,
        "energy_level": 11.72,
        "wall_quenching_prob": 0.01,
        "wall_see_prob": 0.001,
        "comment": "Resonant 1s4 - fast radiative decay"
    },
    {
        "name": "Ar2*",
        "diffusion_coeff": 7.5e-5,
        "mass": 1.326e-25,
        "energy_level": 9.8,
        "wall_quenching_prob": 0.1,
        "wall_see_prob": 0.0,
        "comment": "Excimer - formed from Ar* + Ar + M"
    }
]
```

---

## Example Use Cases

### 1. Dielectric Barrier Discharge (DBD)

**Config:** `config/dbd_argon.json`

**Physics highlights:**
- Memory effect: Ar* accumulates at dielectric
- Lower breakdown voltage in subsequent pulses
- Secondary emission from metastables important

**Key parameters to vary:**
- `dielectric_permittivity` (2.0 - 10.0 for glass/ceramic)
- `wall_see_prob` (0.001 - 0.1)
- `voltage_type = "RF"` with `frequency = 13.56e6`

### 2. Atmospheric Pressure Plasma Jet (APPJ)

**Config:** `config/plasma_jet_argon.json`

**Physics highlights:**
- High Péclet number (Pe ~ 500)
- Advection-dominated transport
- Plasma bullet propagation via excited species

**Key parameters to vary:**
- `gas_velocity` (10 - 200 m/s)
- Domain length `Lx` (0.01 - 0.1 m)
- Number of grid points for resolution

### 3. Penning Mixture (Ne/Ar)

**Config:** `config/penning_mixture.json`

**Physics highlights:**
- Non-local ionization
- Ne* + Ar → Ne + Ar+ + e
- Demonstrates importance of excited transport

**Key species:**
- Ne_m (E = 16.62 eV > E_iz(Ar) = 15.76 eV)

---

## Understanding Output

### Field Names

| Variable | Units | Description |
|----------|-------|-------------|
| `ne` | m⁻³ | Electron density |
| `ni` | m⁻³ | Ion density |
| `neps` | eV·m⁻³ | Electron energy density |
| `phi` | V | Electrostatic potential |
| `sigma` | C/m² | Surface charge (dielectrics only) |
| `Ar_m` | m⁻³ | Metastable density |
| `Ar_r` | m⁻³ | Resonant state density |
| ... | | (other excited species) |

### Derived Quantities

**Mean electron energy:**
```
mean_energy [eV] = neps / ne
```

**Electric field:**
```
E [V/m] = -dφ/dx
```

**Reduced field:**
```
E/N [Td] = E / (N_gas * 1e-21)
where N_gas ~ 3.22e22 m⁻³ at 1 atm, 273 K
```

---

## Troubleshooting

### Convergence Issues

**Symptom:** SNES fails to converge

**Solutions:**
1. Reduce time step: `dt = 1e-13` (start small)
2. Improve initial condition (closer to steady state)
3. Add damping: `-snes_linesearch_type bt`
4. Increase Newton iterations: `-snes_max_it 50`

### Negative Densities

**Symptom:** `ne < 0` or `n* < 0`

**Causes:**
1. Time step too large
2. Reaction rates too stiff

**Solutions:**
1. Reduce `dt`
2. Check that Scharfetter-Gummel is being used (not central differences)
3. Enforce minimum densities in source terms (already in ReactionHandler)

### Unphysical Results

**Symptom:** Plasma doesn't ignite, or densities explode

**Causes:**
1. Incorrect chemistry data
2. Wrong boundary conditions
3. Unrealistic excited species parameters

**Solutions:**
1. Validate transport table with BOLSIG+
2. Check `gamma_see` (typical: 0.01 - 0.2)
3. Verify diffusion coefficients (D ~ 10⁻⁴ m²/s for Ar at atm)

---

## Performance Optimization

### Parallel Execution

```bash
mpirun -n 4 ./HydroPlas config/argon_complete.json
```

### PETSc Options for Speed

```bash
./HydroPlas config.json \
  -ts_type bdf \
  -ts_bdf_order 2 \
  -pc_type fieldsplit \
  -pc_fieldsplit_type multiplicative \
  -fieldsplit_0_pc_type ilu \
  -fieldsplit_1_pc_type lu
```

### Reducing Output Frequency

In config: `"output_interval": 1000` (instead of 10)

---

## Advanced Topics

### Custom Reaction Rates

Edit `src/chemistry/ReactionHandler.cpp` to add specific reactions:

```cpp
void compute_custom_reaction(...) {
    double k_custom = 1e-15; // m³/s
    double R = k_custom * n_species1 * n_species2;
    S_product += R;
    S_reactant1 -= R;
    S_reactant2 -= R;
}
```

### Multi-Dimensional Simulations

Currently 1D+time. 2D extension:
1. Set `domain.Ny > 1` in config
2. Solver automatically handles 2D DMDA
3. May need more grid points (Nx=200, Ny=100)

### Coupling to Flow Solvers

For self-consistent EHD:
1. Compute body force: F = ρE (already in code structure)
2. Pass to external Navier-Stokes solver
3. Update `u_gas` field at each time step

---

## Getting Help

**Issue tracker:** (Add your repository URL)

**Documentation:** See `docs/THEORY.md` for detailed physics

**Examples:** All config files in `config/` are documented

---

**Last Updated:** December 30, 2025
