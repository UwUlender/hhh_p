# Configuration Parameter Management Analysis

## Summary
This document analyzes all parameters from the user's configuration file and identifies which ones are currently managed by the ConfigParser and which are missing.

## Parameters Status

### ✅ MESH Section - FULLY IMPLEMENTED
- `type` - ✅ Parsed
- `x_nodes` - ✅ Parsed
- `y_nodes` - ✅ Parsed

### ⚠️ ELECTRODES Section - PARTIALLY IMPLEMENTED
- `name` - ✅ Parsed
- `location` - ✅ Parsed
- `voltage_type` - ⚠️ **FIELD NAME MISMATCH**: Config uses `voltage_type` but parser looks for `type`
- `voltage_amplitude` - ⚠️ **FIELD NAME MISMATCH**: Config uses `voltage_amplitude` but parser looks for `amplitude`
- `frequency` - ✅ Parsed (but as optional)
- `phase` - ✅ Parsed (but as optional)
- `bias` - ✅ Parsed (but as optional)
- `duty_cycle` - ✅ Parsed (but as optional)
- `voltage_expression` - ✅ Parsed
- `is_dielectric` - ✅ Parsed
- `dielectric_permittivity` - ⚠️ **FIELD NAME MISMATCH**: Config uses `dielectric_permittivity` but parser looks for `permittivity`
- `dielectric_thickness` - ⚠️ **FIELD NAME MISMATCH**: Config uses `dielectric_thickness` but parser looks for `thickness`
- `gamma_see` - ✅ Parsed

### ✅ PLASMA Section - FULLY IMPLEMENTED
- `gas_pressure` - ✅ Parsed
- `background_temp` - ✅ Parsed
- `background_gas` - ✅ Parsed

### ⚠️ SPECIES Section - PARTIALLY IMPLEMENTED
- `name` - ✅ Parsed
- `type` - ✅ Parsed
- `charge` - ✅ Parsed
- `mass` - ✅ Parsed
- `diffusion_coeff` - ✅ Parsed
- `mobility_file` - ✅ Parsed
- `mobility_coeff` - ❌ **MISSING**: Numeric mobility coefficient (alternative to file)

### ❌ REACTIONS Section - CRITICAL ISSUES
Current struct only supports:
- `equation` - ✅ Parsed
- `type` - ✅ Parsed
- `rate_type` - ✅ Parsed
- `a`, `b`, `e_a` - ✅ Parsed (for Arrhenius)
- `table_file` - ✅ Parsed

**Missing parameters:**
- `name_for_output` - ❌ **MISSING**: Friendly name for output
- `energy_change` - ❌ **MISSING**: Energy loss/gain in reaction (eV)
- `equation_constants` - ❌ **MISSING**: Names of constants for equation-based rates
- `equation_values` - ❌ **MISSING**: Values of equation constants
- `equation_variables` - ❌ **MISSING**: Variable names for equation-based rates
- `equation` - ❌ **MISSING**: Mathematical equation for reaction rate

### ❌ INITIAL_CONDITION Section - NOT IMPLEMENTED
**Entire section missing from ConfigParser!**

Required fields:
- `name` - Species name or "e_energy"
- `value` - Initial value (density in 1/m³ or energy in eV)

Example from config:
```yaml
initial_condition:
  - name: e
    value: 1e14
  - name: e_energy
    value: 1.5
  - name: Ar+
    value: 1e14
```

### ❌ SOLVER Section - NOT IMPLEMENTED
**Entire section missing from ConfigParser!**

Required fields:
- `type` - Solver type (e.g., "JFNK")
- `tolerance` - SNES tolerance
- `max_iterations` - Maximum solver iterations
- `time_step` - Time step size in seconds
- `end_time` - Simulation end time in seconds
- `preconditioner` - Preconditioner type (e.g., "PBP")
- `ksp_type` - Krylov solver type (e.g., "GMRES")

**Current workaround**: Hardcoded in main.cpp (dt = 1e-12, t_end = 1e-9)

### ⚠️ OUTPUT Section - PARTIALLY IMPLEMENTED
Current implementation:
- `format` - ✅ Parsed
- `frequency_step` - ✅ Parsed
- `save_rates` - ✅ Parsed
- `filename` - ✅ Parsed (but hardcoded default)

**Missing parameters:**
- `minimum_time_interval` - ❌ **MISSING**: Minimum time between outputs (alternative to frequency_step)
- `save_fields` - ❌ **MISSING**: Array of field names to save
- `rates` - ❌ **MISSING**: Array of reaction rate names to save
- `checkpoint_frequency` - ❌ **MISSING**: Checkpoint save frequency
- `checkpoint_filename` - ❌ **MISSING**: Checkpoint file name
- `add_timestep` - ❌ **MISSING**: Whether to add timestep number to filename

### ❌ RESTART Section - PARTIALLY IMPLEMENTED
**Not in ConfigParser struct, only handled via command-line flags**

Should include:
- `enabled` - Whether restart is enabled
- `file` - Restart file path
- `step` - Restart step number

**Current workaround**: Command-line flags (-restart, -restart_step)

### ❌ ADVANCED Section - NOT IMPLEMENTED
**Entire section missing from ConfigParser!**

Required fields:
- `adaptive_dt` - Enable adaptive time stepping
- `dt_min` - Minimum time step
- `dt_max` - Maximum time step
- `mpi_decomposition` - MPI domain decomposition
- `density_floor` - Minimum density (1/m³)
- `energy_floor` - Minimum electron energy (eV)
- `wall_quenching_probability` - Wall quenching for metastables
- `thermal_velocity_factor` - Factor for thermal flux BCs

## Priority Recommendations

### HIGH PRIORITY (Required for correct simulation):
1. **SOLVER section** - Currently hardcoded, should be configurable
2. **INITIAL_CONDITION section** - Essential for setting up simulations
3. **REACTIONS.energy_change** - Needed for energy balance
4. **SPECIES.mobility_coeff** - Alternative to table-based mobility

### MEDIUM PRIORITY (Improves functionality):
5. **REACTIONS equation-based rates** - For complex reaction kinetics
6. **OUTPUT.save_fields and rates** - Better control over outputs
7. **ADVANCED section** - For simulation tuning

### LOW PRIORITY (Nice to have):
8. **RESTART section in config** - Already works via command-line
9. **OUTPUT.checkpoint settings** - Can use defaults
10. **Field name consistency** - Fix mismatches

## Recommended Actions

1. **Add missing struct definitions** to ConfigParser.hpp
2. **Implement parsing logic** in ConfigParser.cpp
3. **Update main.cpp** to use solver parameters from config
4. **Add validation** for required parameters
5. **Update example configs** to use consistent field names
