# HydroPlas Validation & Benchmark Suite

**Version:** 1.0  
**Last Updated:** December 30, 2025

---

## Overview

This document describes the validation strategy and benchmark test cases for HydroPlas, following the architectural specification requirements. The suite consists of three canonical discharge types that test different aspects of the simulation framework:

1. **Benchmark 1: 1D DC Glow Discharge** - Tests secondary electron emission
2. **Benchmark 2: 1D RF Capacitively Coupled Plasma** - Tests time-dependent boundaries and RF physics
3. **Benchmark 3: 2D Dielectric Barrier Discharge** - Tests surface charging and multi-dimensional transport

---

## Validation Strategy

### Physical Accuracy
- ✅ Scharfetter-Gummel scheme ensures positivity and stability
- ✅ Implicit time integration handles stiff timescales
- ✅ Poisson equation couples transport self-consistently
- ✅ Boundary conditions capture essential physics (SEE, charging, RF)

### Numerical Verification
- ✅ Grid convergence tests (refine Nx, verify solution unchanged)
- ✅ Time step convergence (reduce dt, verify convergence)
- ✅ Mass conservation checks (total charge constant in steady state)
- ✅ Energy conservation (power balance: input = losses)

### Code Verification
- ✅ Unit tests for Scharfetter-Gummel scheme (known analytical solutions)
- ✅ Chemistry integration tests (verify reaction rates)
- ✅ Boundary condition tests (verify flux calculations)

---

## Benchmark 1: 1D DC Glow Discharge in Argon

### Configuration
**File:** `config/benchmark_1_dc.json`

### Physical Setup
- **Discharge type:** DC glow discharge (cathode-anode)
- **Gas:** Argon at 1 Torr (~133 Pa)
- **Geometry:** 1D, gap = 5 cm
- **Voltage:** 300 V DC
- **Critical physics:** Secondary electron emission (γ-mode sustenance)

### Purpose
Tests the **Secondary Electron Emission (SEE)** boundary condition, which is mandatory for sustaining DC discharges. Without SEE, the discharge would extinguish as ions strike the cathode without releasing electrons to continue the avalanche.

### Key Parameters
```json
"domain": {
    "Lx": 0.05,        // 5 cm gap
    "Nx": 200          // Adequate resolution for sheaths
}
"boundary": {
    "voltage_type": "DC",
    "voltage_amplitude": 300.0,   // V
    "gamma_see": 0.01              // SEE coefficient (Ar+ on metal)
}
```

### Expected Results

#### Spatial Structure
1. **Cathode Fall Region** (0 - 5 mm):
   - High electric field (E ~ 10⁴ V/m)
   - Low electron density (ne ~ 10¹⁴ m⁻³)
   - Negative glow emission

2. **Positive Column** (5 mm - 45 mm):
   - Low electric field (E ~ 10² V/m)
   - Quasi-neutral plasma (ne ≈ ni ~ 10¹⁶ m⁻³)
   - Uniform glow

3. **Anode Sheath** (45 mm - 50 mm):
   - Moderate field (E ~ 10³ V/m)
   - Thin compared to cathode fall

#### Validation Metrics
| Quantity | Expected Range | Physical Significance |
|----------|----------------|----------------------|
| Cathode fall voltage | 100-200 V | Energy for ion acceleration |
| Cathode fall thickness | 0.5-2 mm | Collisional mean free path |
| Positive column E-field | 50-500 V/m | Losses compensated by ohmic heating |
| Discharge current | 0.1-10 mA | Depends on γ_see and gap |

#### Success Criteria
✅ Cathode fall exhibits characteristic negative glow  
✅ Positive column shows quasi-neutrality (|ne - ni|/ne < 0.01)  
✅ Steady-state reached (∂n/∂t → 0 after ~10 μs)  
✅ Current-voltage characteristic follows Child-Langmuir law scaling  

---

## Benchmark 2: 1D RF Capacitively Coupled Plasma (CCP)

### Configuration
**File:** `config/benchmark_2_rf.json`

### Physical Setup
- **Discharge type:** RF capacitive (13.56 MHz standard)
- **Gas:** Argon at 100 mTorr (~13 Pa)
- **Geometry:** 1D, gap = 2 cm
- **Voltage:** 200 V RF + 0 V DC bias
- **Critical physics:** Time-dependent boundary, sheath dynamics, stochastic heating

### Purpose
Tests the **Time-Dependent Voltage** boundary condition and RF physics, essential for semiconductor processing plasmas (etching, deposition).

### Key Parameters
```json
"domain": {
    "Lx": 0.02,        // 2 cm gap
    "Nx": 256          // High resolution for thin sheaths
}
"time": {
    "dt": 1e-12,       // Initially, auto-adjusted to T_RF/100
    "t_end": 5e-6      // Multiple RF cycles (~70 periods)
}
"boundary": {
    "voltage_type": "RF",
    "voltage_amplitude": 200.0,   // V
    "frequency": 13.56e6,         // Hz (standard RF freq)
    "bias": 0.0                   // V (symmetric case)
}
```

### Automatic Time Step Control
HydroPlas automatically detects RF mode and adjusts time step:
```
T_RF = 1 / 13.56 MHz = 73.7 ns
dt_auto = T_RF / 100 = 0.737 ns
```

This ensures ~100 points per RF cycle, resolving:
- Sheath oscillations
- Electron heating phases
- Ion transit time effects

### Expected Results

#### Time-Averaged Structure
1. **Sheaths** (oscillating, 0-2 mm each side):
   - Time-varying thickness: d_sheath(t) ~ λ_D √(V_RF sin(ωt))
   - Electron heating occurs during sheath collapse

2. **Bulk Plasma** (2 mm - 18 mm):
   - Quasi-uniform density (ne ~ 10¹⁶ m⁻³)
   - Low time-averaged field

#### RF-Specific Phenomena
1. **Stochastic Heating:**
   - Electrons gain energy as sheaths collapse
   - Mean energy oscillates at 2×RF frequency

2. **DC Self-Bias:**
   - Even with symmetric electrodes, slight asymmetry develops
   - Causes time-averaged potential offset

3. **Sheath Voltage Division:**
   - Most voltage drops across sheaths, not bulk

#### Validation Metrics
| Quantity | Expected Range | Physical Significance |
|----------|----------------|----------------------|
| Time-avg plasma density | 10¹⁵-10¹⁷ m⁻³ | Power coupling efficiency |
| Electron temperature | 2-5 eV | Collisionless heating |
| Sheath thickness | 0.1-1 mm | Debye length scaling |
| Ionization rate peak location | Sheath edge | Spatial pattern of etching |

#### Success Criteria
✅ Density oscillates at 2×f_RF (ion inertia filters RF)  
✅ Electron temperature sustained despite E passing through zero  
✅ Time-averaged spatial profiles match literature  
✅ Power balance: P_RF = P_losses + P_output  

---

## Benchmark 3: 2D Dielectric Barrier Discharge (DBD)

### Configuration
**File:** `config/benchmark_3_dbd.json`

### Physical Setup
- **Discharge type:** Atmospheric pressure DBD
- **Gas:** Argon at 1 atm
- **Geometry:** 2D rectangular, 1 cm × 1 mm
- **Voltage:** 5 kV AC, 1 kHz
- **Critical physics:** Dielectric surface charging, self-pulsing, memory effect

### Purpose
Tests the **Surface Charge Accumulation** boundary condition, mandatory for DBDs. The discharge self-extinguishes due to charge buildup on the dielectric, then reignites when polarity reverses—a key feature for ozone generation and surface treatment.

### Key Parameters
```json
"domain": {
    "Lx": 0.01,        // 1 cm length
    "Ly": 0.001,       // 1 mm gap
    "Nx": 100,
    "Ny": 50
}
"time": {
    "dt": 1e-11,       // High pressure requires small dt
    "t_end": 5e-3      // Multiple AC cycles
}
"boundary": {
    "voltage_type": "RF",          // Using RF for AC (low freq)
    "voltage_amplitude": 5000.0,   // V
    "frequency": 1000.0,           // Hz (1 kHz)
    "dielectric_permittivity": 5.0,  // Glass/ceramic
    "dielectric_thickness": 0.001    // 1 mm
}
```

### Expected Results

#### Temporal Behavior
1. **First Voltage Cycle:**
   - Discharge ignites at V_breakdown ~ 3-4 kV
   - Microdischarges form randomly across electrode
   - Charge accumulates on dielectric (σ → +)
   - Discharge extinguishes as field drops

2. **Subsequent Cycles:**
   - V_breakdown decreases (memory effect!)
   - V_breakdown,2nd ~ 2 kV (40% lower due to residual charge)
   - Discharge becomes more uniform spatially

#### Spatial Structure (2D)
1. **Microdischarge Filaments:**
   - Diameter ~ 0.1-0.5 mm
   - Peak density ~ 10²⁰ m⁻³ (very high due to atm pressure)
   - Propagate toward dielectric at ~10⁵ m/s

2. **Surface Charge Pattern:**
   - σ(x,y) non-uniform, follows filament pattern
   - Persists between pulses (ms timescale)
   - Shields electric field, causing extinction

#### Validation Metrics
| Quantity | Expected Range | Physical Significance |
|----------|----------------|----------------------|
| Breakdown voltage (1st) | 3-5 kV | Paschen curve (pd product) |
| Breakdown voltage (2nd) | 1.5-3 kV | Memory effect magnitude |
| Pulse duration | 10-100 ns | Limited by charge accumulation |
| Surface charge density | 10⁻⁶-10⁻⁴ C/m² | Dielectric charging efficiency |

#### Success Criteria
✅ Discharge self-extinguishes due to surface charging  
✅ Breakdown voltage decreases in subsequent pulses (memory)  
✅ 2D structure shows filamentary microdischarges  
✅ Power dissipation matches Lissajous figure measurements  

---

## Validation Workflow

### Step 1: Grid Convergence Study
For each benchmark, run with progressively finer grids:

```bash
# Coarse
./HydroPlas config/benchmark_1_dc.json

# Medium (modify config: Nx = 400)
./HydroPlas config/benchmark_1_dc_fine.json

# Fine (Nx = 800)
./HydroPlas config/benchmark_1_dc_finest.json
```

**Acceptance:** Solutions differ by < 5% between medium and fine grids

---

### Step 2: Time Step Convergence
Run with dt, dt/2, dt/4. Verify transient behavior converges.

---

### Step 3: Physical Checks

#### Mass Conservation
```python
# Python post-processing
import h5py
import numpy as np

f = h5py.File('output/hydroplas_final.h5', 'r')
ne = f['/data/final/meshes/ne'][:]
ni = f['/data/final/meshes/ni'][:]

total_charge = np.sum(ni - ne) * dx * dy
print(f"Total charge: {total_charge} (should be ~ 0)")
```

#### Energy Balance
```python
# For RF: Average power in = Average power out
P_RF = V_RF * I_RF * cos(phase)
P_losses = Integral(n_e * nu_coll * (3/2) * k * T_e) * volume
assert abs(P_RF - P_losses) / P_RF < 0.1
```

#### Boundary Condition Verification
```python
# SEE check: Verify γ electrons per ion
flux_i_cathode = # Extract from simulation
flux_e_cathode = # Extract from simulation
gamma_measured = -flux_e_cathode / flux_i_cathode
assert abs(gamma_measured - gamma_config) < 0.05
```

---

## Comparison with Literature

### DC Glow Discharge
**Reference:** Phelps & Petrović (1999), "Cold-cathode discharges and breakdown in argon"

| Parameter | HydroPlas | Literature | Match |
|-----------|-----------|------------|-------|
| Cathode fall voltage | TBD | 150 V | ✓/✗ |
| Positive column E-field | TBD | 200 V/m | ✓/✗ |
| Current density | TBD | 1 mA/cm² | ✓/✗ |

---

### RF CCP
**Reference:** Turner et al. (2013), "Simulation benchmarks for low-pressure plasmas"

| Parameter | HydroPlas | GEC Reference | Match |
|-----------|-----------|---------------|-------|
| Plasma density | TBD | 2×10¹⁶ m⁻³ | ✓/✗ |
| Electron temperature | TBD | 3 eV | ✓/✗ |
| Sheath width | TBD | 0.5 mm | ✓/✗ |

---

### DBD
**Reference:** Becker et al. (2005), "Microplasmas and applications"

| Parameter | HydroPlas | Experimental | Match |
|-----------|-----------|--------------|-------|
| Pulse duration | TBD | 10-50 ns | ✓/✗ |
| V_breakdown ratio | TBD | 0.5-0.7 | ✓/✗ |
| Filament diameter | TBD | 0.1-0.3 mm | ✓/✗ |

---

## Running the Benchmark Suite

### Quick Validation
```bash
cd HydroPlas/build

# Benchmark 1: DC Glow (~1 minute)
./HydroPlas ../config/benchmark_1_dc.json -ts_monitor

# Benchmark 2: RF CCP (~5 minutes)
./HydroPlas ../config/benchmark_2_rf.json -ts_monitor

# Benchmark 3: DBD 2D (~30 minutes)
./HydroPlas ../config/benchmark_3_dbd.json -ts_monitor
```

### With Detailed Monitoring
```bash
./HydroPlas config.json \
  -ts_monitor \
  -snes_monitor \
  -ksp_monitor \
  -log_view
```

### Parallel Execution (for 2D)
```bash
mpirun -n 4 ./HydroPlas ../config/benchmark_3_dbd.json
```

---

## Automated Testing (Future)

### Continuous Integration
```bash
# Proposed CI/CD pipeline
./run_benchmarks.sh  # Runs all 3 benchmarks
python validate_results.py  # Checks against reference data
```

### Regression Tests
- Ensure updates don't break existing functionality
- Compare key metrics (density, voltage profiles) against baseline

---

## Troubleshooting Benchmark Failures

### Benchmark 1: DC Discharge Won't Ignite
**Symptoms:** Electron density remains at seed level

**Possible causes:**
1. γ_see too low → Increase to 0.05-0.1
2. Voltage too low → Check Paschen curve for your pd
3. Grid too coarse → Cathode fall requires Nx > 100

---

### Benchmark 2: RF Plasma Unstable
**Symptoms:** SNES convergence failures, density oscillates wildly

**Possible causes:**
1. Time step too large → Use auto-adjusted dt (T_RF/100)
2. Pressure too high/low → RF works best at 10-100 mTorr
3. Preconditioning issue → Try `-pc_fieldsplit_type additive`

---

### Benchmark 3: DBD Never Extinguishes
**Symptoms:** Discharge continuous, no pulsing behavior

**Possible causes:**
1. Dielectric BC not applied → Check `dielectric_permittivity > 1`
2. Surface charge equation not solved → Verify σ field in output
3. Capacitance too low → Increase `dielectric_thickness`

---

## References

1. **DC Glow Discharges:**
   - Raizer, Y.P. (1991). *Gas Discharge Physics*. Springer.
   - Phelps & Petrović (1999). *Plasma Sources Sci. Technol.* 8, R21.

2. **RF Plasmas:**
   - Lieberman & Lichtenberg (2005). *Principles of Plasma Discharges*. Wiley.
   - Turner et al. (2013). *Phys. Plasmas* 20, 013507.

3. **Dielectric Barrier Discharges:**
   - Kogelschatz (2003). *Plasma Chem. Plasma Process.* 23, 1.
   - Becker et al. (2005). *J. Phys. D: Appl. Phys.* 39, R55.

---

## Conclusion

The HydroPlas validation suite provides:
✅ Three canonical test cases covering all key physics  
✅ Clear success criteria and acceptance metrics  
✅ Comparison framework with literature benchmarks  
✅ Debugging procedures for common issues  

**Status:** Test cases configured, awaiting execution and result comparison.

---

**Next Steps:**
1. Execute all three benchmarks
2. Extract key metrics from output
3. Compare with literature values
4. Document discrepancies and refine models
5. Establish baseline for regression testing

**Validation Target Date:** Q1 2026
