# HydroPlas Quick Reference Guide

**Version:** 1.0  
**Last Updated:** December 30, 2025

---

## 🎯 One-Minute Quick Start

```bash
# Build
cd HydroPlas/build
cmake .. && make -j4

# Run complete Argon chemistry example
./HydroPlas ../config/argon_complete.json

# Output: ./output/*.txt (text) or *.h5 (HDF5)
```

---

## 📋 Configuration Cheat Sheet

### Minimal Excited Species Config

```json
{
  "excited_species": [
    {
      "name": "Ar_m",                 // Species identifier
      "diffusion_coeff": 1.5e-4,      // D* [m²/s]
      "mass": 6.63e-26,               // m* [kg]
      "energy_level": 11.55,          // E* [eV]
      "wall_quenching_prob": 0.001,   // γ_quench (0-1)
      "wall_see_prob": 0.01           // γ_see (electrons per atom)
    }
  ]
}
```

### Quick Parameter Guide

| Parameter | Symbol | Typical Range | Units |
|-----------|--------|---------------|-------|
| Diffusion coeff | D* | 10⁻⁵ - 10⁻³ | m²/s |
| Mass | m* | 10⁻²⁶ - 10⁻²⁵ | kg |
| Energy level | E* | 9 - 17 | eV |
| Wall quenching | γ_q | 10⁻⁴ - 0.5 | - |
| Wall SEE | γ_see | 0 - 0.1 | - |
| Gas velocity | u_gas | 0 - 200 | m/s |
| Gas temperature | T_gas | 200 - 500 | K |

---

## 🧪 Example Configurations

### 1. Test Case (Single Metastable)
```bash
./HydroPlas config/excited_test.json
```
- **Purpose:** Verify ADR equations work
- **Physics:** Simple advection-diffusion
- **Runtime:** ~30 seconds

### 2. Complete Argon (3 Species)
```bash
./HydroPlas config/argon_complete.json
```
- **Purpose:** Full chemistry benchmark
- **Species:** Ar_m, Ar_r, Ar2*
- **Runtime:** ~2 minutes

### 3. DBD (Memory Effect)
```bash
./HydroPlas config/dbd_argon.json
```
- **Purpose:** Dielectric barrier discharge
- **Physics:** Metastable accumulation, reduced V_breakdown
- **Runtime:** ~1 minute

### 4. Plasma Jet (High Pe)
```bash
./HydroPlas config/plasma_jet_argon.json
```
- **Purpose:** Atmospheric pressure jet
- **Physics:** Advection-dominated (Pe ~ 500), plasma bullet
- **Runtime:** ~5 minutes (long domain)

### 5. Penning Mixture
```bash
./HydroPlas config/penning_mixture.json
```
- **Purpose:** Ne/Ar non-local ionization
- **Physics:** Ne* + Ar → Ar+ + e
- **Runtime:** ~1 minute

---

## 🔧 Common PETSc Options

### Monitoring
```bash
-ts_monitor          # Time step info
-snes_monitor        # Newton convergence
-ksp_monitor         # Linear solver
```

### Solver Control
```bash
-ts_type bdf                    # Backward differentiation (default)
-ts_dt 1e-12                    # Force initial time step
-snes_max_it 50                 # Max Newton iterations
-snes_linesearch_type bt        # Backtracking line search
```

### Preconditioning
```bash
-pc_type fieldsplit             # Block preconditioner (default in code)
-fieldsplit_0_pc_type ilu       # ILU for transport
-fieldsplit_1_pc_type lu        # Direct solve for Poisson
```

### Performance
```bash
-snes_mf_operator               # Jacobian-free (for many species)
-log_view                       # Performance profiling
```

---

## 📊 Output Files

### Text Output (Always Generated)
```
output/ne_000100.txt        # Electron density at step 100
output/ni_000100.txt        # Ion density
output/neps_000100.txt      # Electron energy density
output/phi_000100.txt       # Potential
output/Ar_m_000100.txt      # Metastable density
...
```

**Format:**
```
# Time: 1.0e-9
# Field: ne
# Grid: 200 x 1
0 0 1.234e15
1 0 1.235e15
...
```

### HDF5 Output (If Compiled with HDF5)
```
output/hydroplas_000100.h5
├── /data/100/
│   ├── @time = 1.0e-9
│   └── /meshes/
│       ├── ne (dataset)
│       ├── Ar_m (dataset)
│       └── ...
```

**View with:**
```bash
# ParaView
paraview output/hydroplas_*.h5

# Python
import h5py
f = h5py.File('output/hydroplas_000100.h5')
ne = f['/data/100/meshes/ne'][:]

# Command line
h5dump -H output/hydroplas_000100.h5
```

---

## 🐛 Troubleshooting Quick Fixes

### Problem: SNES Not Converging
```bash
# Solution 1: Reduce time step
"dt": 1e-13

# Solution 2: Damping
./HydroPlas config.json -snes_linesearch_type bt -snes_linesearch_damping 0.9

# Solution 3: More iterations
-snes_max_it 100
```

### Problem: Negative Densities
```bash
# Cause: Time step too large or bad initial condition
# Solution:
"dt": 1e-14  # Start very small
# Increase gradually once stable
```

### Problem: Plasma Won't Ignite
```bash
# Check:
1. gamma_see too low? (try 0.1)
2. Voltage too low? (increase by 50V)
3. Chemistry table missing? (check transport.dat exists)
```

### Problem: Simulation Too Slow
```bash
# Solution 1: Coarser grid
"Nx": 50  # instead of 200

# Solution 2: Less output
"output_interval": 1000  # instead of 10

# Solution 3: Parallel
mpirun -n 4 ./HydroPlas config.json
```

---

## 📐 Physics Formulas Reference

### Péclet Number
```
Pe = (u_gas · L) / D*
```
- Pe < 1: Diffusion-dominated
- Pe > 10: Advection-dominated

### Mean Electron Energy
```
mean_energy [eV] = neps / ne
```

### Electric Field
```
E [V/m] = -dφ/dx
```

### Reduced Field
```
E/N [Td] = E / (N_gas · 1e-21)
1 Td = 10⁻²¹ V·m²
```

### Thermal Velocity
```
v_th = √(8 k_B T / π m)
```

### Wall Flux
```
Γ_wall = (γ · v_th / 4) · n  +  u_gas · n
         ↑                      ↑
      Thermal              Advection
```

---

## 🧬 Reaction Rate Equations

### Direct Ionization
```
R_iz = k_iz(Te) · ne · N_gas
E_cost = 15.76 eV (Argon)
```

### Stepwise Ionization
```
R_step = k_step(Te) · ne · n*
E_cost = 4.2 eV (Argon from Ar*)
```

### Penning/Pooling
```
R_penning = k_p · n* · n*
Rate coeff ~ 10⁻¹⁵ m³/s
```

### Radiative Decay
```
R_rad = A_rad · n*
A_rad(metastable) ~ 0
A_rad(resonant) ~ 10⁸ s⁻¹
```

---

## 📚 Documentation Roadmap

1. **First time user?** → Start with [README.md](README.md)
2. **Want to run simulations?** → See [docs/USER_GUIDE.md](docs/USER_GUIDE.md)
3. **Need physics details?** → Read [docs/THEORY.md](docs/THEORY.md)
4. **Implementing changes?** → Check [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md)
5. **Version info?** → See [docs/CHANGELOG.md](docs/CHANGELOG.md)

---

## 🎓 Argon Species Quick Reference

### Ar_m (Metastable 1s₅)
- **Energy:** 11.55 eV
- **Lifetime:** ~1 ms (long)
- **Radiative decay:** Forbidden
- **Key role:** Stepwise ionization, Penning, memory effect
- **Diffusion:** 1.5×10⁻⁴ m²/s

### Ar_r (Resonant 1s₄)
- **Energy:** 11.72 eV
- **Lifetime:** ~10 ns (short)
- **Radiative decay:** Fast (A ~ 10⁸ s⁻¹)
- **Key role:** Photon emission, radiation trapping
- **Diffusion:** 1.5×10⁻⁴ m²/s

### Ar₂* (Excimer)
- **Energy:** 9.8 eV
- **Formation:** Ar* + Ar + M → Ar₂* + M
- **Lifetime:** ~few μs
- **Key role:** VUV emission (126 nm), surface treatment
- **Diffusion:** 7.5×10⁻⁵ m²/s (heavier)

---

## ⚡ Performance Benchmarks

| Grid Size | Species | Time Steps | Wall Time | Memory |
|-----------|---------|------------|-----------|--------|
| 100 pts   | 3       | 1000       | 20 sec    | 50 MB  |
| 200 pts   | 3       | 1000       | 1 min     | 80 MB  |
| 500 pts   | 5       | 1000       | 5 min     | 200 MB |
| 200×100   | 3       | 1000       | 30 min    | 1 GB   |

*Single core, Intel Xeon 2.4 GHz*

---

## 🔗 Quick Links

- **PETSc Manual:** https://petsc.org/release/manual/
- **BOLSIG+ Download:** https://www.bolsig.laplace.univ-tlse.fr/
- **LXCat Database:** https://www.lxcat.net/
- **OpenPMD Standard:** https://github.com/openPMD/openPMD-standard

---

## 💡 Pro Tips

1. **Start small:** Use `excited_test.json` first, then scale up
2. **Check chemistry:** Validate BOLSIG+ output before long runs
3. **Monitor convergence:** Always use `-ts_monitor -snes_monitor` for first run
4. **Output strategy:** Low frequency for production (interval=1000), high for debugging (interval=10)
5. **Visualization:** HDF5 output is much faster to load than text for large datasets

---

**Need more help?** See the full [User Guide](docs/USER_GUIDE.md) or open an issue!
