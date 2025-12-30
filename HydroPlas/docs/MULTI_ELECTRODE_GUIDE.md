# Multi-Electrode Boundary Conditions Guide

## Overview

HydroPlas now supports flexible per-electrode voltage boundary conditions, allowing you to specify independent voltage waveforms for each electrode in your plasma simulation. This enables accurate modeling of complex discharge configurations.

## Features

### Supported Voltage Types

1. **DC (Direct Current)**: Constant voltage
2. **RF (Radio Frequency)**: Sinusoidal AC voltage
3. **AC (Alternating Current)**: Same as RF, sinusoidal waveform
4. **PULSE**: Square wave with configurable duty cycle

### Per-Electrode Parameters

Each electrode can be configured independently with:
- **Voltage type** (DC, RF, AC, PULSE)
- **Amplitude** (V)
- **Frequency** (Hz, for RF/AC/PULSE)
- **DC bias** (V, added to all waveforms)
- **Phase** (radians, for RF/AC)
- **Duty cycle** (0-1, for PULSE type)
- **Secondary electron emission coefficient** (γ_see)
- **Dielectric properties** (enable/disable, permittivity, thickness)

## Configuration Format

### Multi-Electrode Configuration

```json
{
    "boundary": {
        "electrodes": [
            {
                "name": "left",
                "voltage_type": "RF",
                "voltage_amplitude": 300.0,
                "frequency": 13.56e6,
                "bias": 0.0,
                "phase": 0.0,
                "duty_cycle": 0.5,
                "gamma_see": 0.1,
                "is_dielectric": false,
                "dielectric_permittivity": 1.0,
                "dielectric_thickness": 0.0
            },
            {
                "name": "right",
                "voltage_type": "DC",
                "voltage_amplitude": 0.0,
                "frequency": 0.0,
                "bias": 0.0,
                "phase": 0.0,
                "duty_cycle": 0.5,
                "gamma_see": 0.1,
                "is_dielectric": false,
                "dielectric_permittivity": 1.0,
                "dielectric_thickness": 0.0
            }
        ]
    }
}
```

### Legacy Single-Electrode Configuration (Still Supported)

The original configuration format is still supported for backward compatibility:

```json
{
    "boundary": {
        "voltage_type": "DC",
        "voltage_amplitude": 300.0,
        "frequency": 0.0,
        "bias": 0.0,
        "gamma_see": 0.1,
        "dielectric_permittivity": 1.0,
        "dielectric_thickness": 0.0
    }
}
```

## Voltage Waveform Equations

### DC Voltage
```
V(t) = V_amplitude + V_bias
```

### RF/AC Voltage
```
V(t) = V_amplitude × sin(2π × f × t + φ) + V_bias
```
Where:
- f = frequency (Hz)
- φ = phase (radians)

### PULSE Voltage
```
V(t) = V_amplitude + V_bias  (when t_mod < T × duty_cycle)
V(t) = V_bias                 (otherwise)
```
Where:
- T = 1/f (period)
- t_mod = t mod T

## Example Configurations

### 1. DC-DC Configuration
Two electrodes with different constant voltages:
- Left: 500 V DC
- Right: 200 V DC

See: `config/multi_electrode_dc_dc.json`

### 2. RF-Ground Configuration
Typical capacitively coupled plasma (CCP) configuration:
- Left: 250 V RF at 13.56 MHz
- Right: Grounded (0 V)

See: `config/multi_electrode_rf_ground.json`

### 3. Push-Pull RF Configuration
Two RF electrodes with 180° phase shift:
- Left: 300 V RF at 13.56 MHz, phase = 0
- Right: 300 V RF at 13.56 MHz, phase = π

See: `config/multi_electrode_rf_rf_phase.json`

### 4. Dual-Frequency Discharge
Different RF frequencies on each electrode:
- Left: 200 V at 13.56 MHz (HF)
- Right: 150 V at 2 MHz (LF)

See: `config/multi_electrode_dual_freq.json`

### 5. Dielectric Barrier Discharge (DBD)
Both electrodes with dielectric barriers:
- Left: 3000 V RF, dielectric (ε_r=4.0, d=1mm)
- Right: Grounded, dielectric (ε_r=3.5, d=1.5mm)

See: `config/multi_electrode_dbd_dual_dielectric.json`

### 6. Pulsed Discharge
Pulsed voltage on one electrode:
- Left: 800 V PULSE at 100 kHz, duty cycle = 30%, bias = -100V
- Right: 50 V DC

See: `config/multi_electrode_pulse_dc.json`

## Technical Implementation

### Electrode Indexing
- **Left boundary (i=0)**: Electrode index 0
- **Right boundary (i=M-1)**: Electrode index 1

### Boundary Condition Application

For each electrode, the code applies:

1. **Voltage boundary condition** in Poisson equation:
   - Conducting electrode: φ = V(t)
   - Dielectric electrode: ε₀E_plasma = ε_d E_dielectric + σ

2. **Particle flux boundary conditions**:
   - Ion flux to wall
   - Electron flux (thermal + field-driven)
   - Secondary electron emission from ions: γ_see × Γ_ion
   - Secondary electron emission from excited species: γ_see_k × Γ_k

3. **Surface charge accumulation** (for dielectrics):
   - dσ/dt = J_wall

## API Reference

### BoundaryManager Methods

```cpp
// Get voltage for specific electrode
double get_electrode_voltage_by_index(int electrode_index, double t);
double get_electrode_voltage(const std::string& electrode_name, double t);

// Get secondary emission coefficient
double get_electrode_gamma_see_by_index(int electrode_index);
double get_electrode_gamma_see(const std::string& electrode_name);

// Check dielectric properties
bool is_electrode_dielectric_by_index(int electrode_index);
bool is_electrode_dielectric(const std::string& electrode_name);
double get_electrode_dielectric_permittivity(const std::string& electrode_name);
double get_electrode_dielectric_thickness(const std::string& electrode_name);

// Get number of configured electrodes
int get_num_electrodes();
```

## Validation and Best Practices

### Physical Considerations

1. **Voltage Amplitude**: 
   - Low-pressure discharges: typically 100-1000 V
   - Atmospheric DBD: 1-10 kV
   - Ensure breakdown voltage is exceeded

2. **Frequency Selection**:
   - Low frequency (50-500 kHz): More ion heating
   - RF (1-100 MHz): More electron heating
   - Dual frequency: Independent control of ion and electron energies

3. **Phase Relationships**:
   - In-phase (0°): Symmetric heating
   - Out-of-phase (180°): Asymmetric heating, DC self-bias
   - 90°: Circular/elliptical electric field (2D)

4. **Dielectric Barriers**:
   - Limits current, prevents arcing
   - Creates memory effect (surface charges persist)
   - Typical materials: glass (ε_r ≈ 4-7), ceramic (ε_r ≈ 9-10)

### Simulation Guidelines

1. **Time step**: Must resolve fastest timescale
   - RF period: dt < T/50 (e.g., dt < 1 ns for 13.56 MHz)
   - Plasma frequency: dt < 1/ω_pe
   - Dielectric relaxation: dt < ε/σ

2. **Spatial resolution**:
   - Resolve sheath: dx < λ_D (Debye length)
   - Typically dx ~ 10-100 μm

3. **Output frequency**:
   - Sample multiple RF cycles
   - For dual frequency, sample enough to see beating patterns

## Troubleshooting

### Common Issues

1. **Simulation instability with RF**:
   - Reduce time step
   - Check that dt << 1/(2πf)

2. **No discharge ignition**:
   - Check voltage amplitude (must exceed breakdown)
   - Verify initial conditions (seed electrons)
   - Check gas pressure and gap distance

3. **Unphysical results with dielectrics**:
   - Verify dielectric thickness > 0
   - Check that ε_r > 1
   - Ensure surface charge equation is being solved

4. **Phase effects not visible**:
   - Simulation time must be long enough (>> 1/Δf for beat frequency)
   - Output interval must be fine enough to resolve phase

## References

1. Lieberman & Lichtenberg, "Principles of Plasma Discharges and Materials Processing", 2nd Ed.
2. Raizer, "Gas Discharge Physics"
3. Adamovich et al., "The 2017 Plasma Roadmap: Low temperature plasma science and technology", J. Phys. D: Appl. Phys. 50 323001 (2017)

## Version History

- **v1.0** (2025-12-30): Initial multi-electrode implementation
  - Support for DC, RF, AC, and PULSE waveforms
  - Per-electrode dielectric properties
  - Phase control for multi-electrode RF
  - Secondary electron emission per electrode
