# Multi-Electrode Quick Start Guide

## TL;DR

HydroPlas now supports **independent voltage control for each electrode**. You can set different voltages (DC, RF, AC, or PULSE) on the left and right boundaries.

## Quick Example

### Old Way (Single Voltage)
```json
{
    "boundary": {
        "voltage_type": "DC",
        "voltage_amplitude": 300.0
    }
}
```

### New Way (Multi-Electrode)
```json
{
    "boundary": {
        "electrodes": [
            {
                "name": "left",
                "voltage_type": "RF",
                "voltage_amplitude": 250.0,
                "frequency": 13.56e6,
                "bias": 0.0,
                "phase": 0.0,
                "gamma_see": 0.1,
                "is_dielectric": false
            },
            {
                "name": "right",
                "voltage_type": "DC",
                "voltage_amplitude": 0.0,
                "gamma_see": 0.08,
                "is_dielectric": false
            }
        ]
    }
}
```

## Voltage Types

| Type | Description | Use Case |
|------|-------------|----------|
| `DC` | Constant voltage | Standard DC discharge, grounded electrode |
| `RF` | Sinusoidal (13.56 MHz typical) | Capacitive coupling, plasma etching |
| `AC` | Same as RF | Any frequency sinusoid |
| `PULSE` | Square wave | Pulsed discharges, DBD |

## Common Configurations

### 1. Capacitively Coupled Plasma (CCP)
**One RF electrode, one grounded**
```json
"electrodes": [
    {"name": "left", "voltage_type": "RF", "voltage_amplitude": 200.0, "frequency": 13.56e6},
    {"name": "right", "voltage_type": "DC", "voltage_amplitude": 0.0}
]
```
File: `config/multi_electrode_rf_ground.json`

### 2. Dual-Frequency CCP
**Different frequencies for ion vs electron control**
```json
"electrodes": [
    {"name": "left", "voltage_type": "RF", "voltage_amplitude": 200.0, "frequency": 13.56e6},
    {"name": "right", "voltage_type": "RF", "voltage_amplitude": 150.0, "frequency": 2.0e6}
]
```
File: `config/multi_electrode_dual_freq.json`

### 3. Push-Pull (180° Phase Shift)
**Maximum voltage swing across gap**
```json
"electrodes": [
    {"name": "left", "voltage_type": "RF", "voltage_amplitude": 300.0, "phase": 0.0},
    {"name": "right", "voltage_type": "RF", "voltage_amplitude": 300.0, "phase": 3.14159}
]
```
File: `config/multi_electrode_rf_rf_phase.json`

### 4. Dielectric Barrier Discharge
**High voltage with dielectric barriers**
```json
"electrodes": [
    {
        "name": "left",
        "voltage_type": "RF",
        "voltage_amplitude": 3000.0,
        "is_dielectric": true,
        "dielectric_permittivity": 4.0,
        "dielectric_thickness": 0.001
    },
    {
        "name": "right",
        "voltage_type": "DC",
        "voltage_amplitude": 0.0,
        "is_dielectric": true,
        "dielectric_permittivity": 3.5,
        "dielectric_thickness": 0.0015
    }
]
```
File: `config/multi_electrode_dbd_dual_dielectric.json`

### 5. Pulsed Discharge
**Square wave with duty cycle**
```json
"electrodes": [
    {
        "name": "left",
        "voltage_type": "PULSE",
        "voltage_amplitude": 800.0,
        "frequency": 1e5,
        "duty_cycle": 0.3,
        "bias": -100.0
    },
    {"name": "right", "voltage_type": "DC", "voltage_amplitude": 50.0}
]
```
File: `config/multi_electrode_pulse_dc.json`

## All Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `name` | string | required | Electrode identifier ("left", "right") |
| `voltage_type` | string | "DC" | "DC", "RF", "AC", or "PULSE" |
| `voltage_amplitude` | double | 0.0 | Amplitude in Volts |
| `frequency` | double | 0.0 | Frequency in Hz (RF/AC/PULSE) |
| `bias` | double | 0.0 | DC offset in Volts |
| `phase` | double | 0.0 | Phase in radians (RF/AC) |
| `duty_cycle` | double | 0.5 | Fraction ON (PULSE, 0-1) |
| `gamma_see` | double | 0.1 | Secondary emission coefficient |
| `is_dielectric` | bool | false | Enable dielectric barrier |
| `dielectric_permittivity` | double | 1.0 | Relative permittivity (ε_r) |
| `dielectric_thickness` | double | 0.0 | Thickness in meters |

## Voltage Formulas

**DC:**
```
V(t) = voltage_amplitude + bias
```

**RF/AC:**
```
V(t) = voltage_amplitude × sin(2π × frequency × t + phase) + bias
```

**PULSE:**
```
V(t) = voltage_amplitude + bias  (when ON, fraction = duty_cycle)
V(t) = bias                       (when OFF)
```

## Running

```bash
# Use any multi-electrode config
./HydroPlas config/multi_electrode_rf_ground.json

# With PETSc options
./HydroPlas config/multi_electrode_dual_freq.json -ts_monitor -snes_monitor

# Parallel
mpirun -n 4 ./HydroPlas config/multi_electrode_rf_rf_phase.json
```

## Backward Compatibility

✅ **Old configurations still work!** No need to update existing files unless you want multi-electrode features.

Old format:
```json
{"boundary": {"voltage_type": "DC", "voltage_amplitude": 300.0}}
```
This is automatically interpreted as a single electrode on the left boundary.

## Tips

1. **Time step for RF**: Use `dt < 1/(50 × frequency)`
   - For 13.56 MHz: `dt < 1.5e-9` seconds

2. **Phase for symmetric heating**: Use same phase (0) on both electrodes

3. **Phase for asymmetric**: 180° shift (`phase = 3.14159`) for maximum asymmetry

4. **Dual-frequency**: Choose frequencies with ratio > 5:1 for clear separation

5. **Dielectric**: Always set `thickness > 0` when `is_dielectric = true`

6. **Secondary emission**: Typical values 0.01 - 0.2 for most materials

## More Information

- **Full Guide**: `docs/MULTI_ELECTRODE_GUIDE.md`
- **Implementation Details**: `MULTI_ELECTRODE_IMPLEMENTATION.md`
- **Example Configs**: `config/multi_electrode_*.json`

## Support

Open an issue if you encounter problems or need help with custom configurations.

---

**Feature Status**: ✅ Production Ready  
**Added**: December 30, 2025  
**Backward Compatible**: Yes
