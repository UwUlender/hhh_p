# Configuration Parameter Management - Implementation Complete

## Overview
This document summarizes the comprehensive update to HydroPlas's configuration management system. All parameters from the user's configuration file specification are now properly parsed and managed.

## Changes Made

### 1. Updated ConfigParser.hpp

Added the following new struct definitions and fields:

#### New Structs Added:
- **InitialConditionConfig** - For setting initial species densities and electron energy
- **SolverConfig** - For solver parameters (previously hardcoded in main.cpp)
- **RestartConfig** - For restart/checkpoint configuration
- **AdvancedConfig** - For advanced simulation tuning parameters

#### Enhanced Existing Structs:

**SpeciesConfig:**
- Added `mobility_coeff` - Numeric mobility coefficient (alternative to mobility_file)

**ReactionConfig:**
- Added `name_for_output` - Friendly name for output files
- Added `energy_change` - Energy loss/gain in eV for electron energy balance
- Added `equation_constants` - For equation-based reaction rates
- Added `equation_values` - Values for equation constants
- Added `equation_variables` - Variable names for equation-based rates
- Added `rate_equation` - Mathematical equation string for rate calculation

**OutputConfig:**
- Added `minimum_time_interval` - Minimum time between outputs (alternative to frequency_step)
- Added `save_fields` - Vector of field names to save
- Added `rates` - Vector of reaction rate names to save
- Added `add_timestep` - Boolean to add timestep number to filenames
- Added `checkpoint_frequency` - Save checkpoint every N steps
- Added `checkpoint_filename` - Checkpoint file name

**SimulationConfig:**
- Added `initial_conditions` - Vector of initial condition configurations
- Added `solver` - Solver configuration struct
- Added `restart` - Restart configuration struct
- Added `advanced` - Advanced configuration struct

### 2. Updated ConfigParser.cpp

Implemented parsing for all new parameters with the following features:

#### Electrode Section Enhancements:
- Support for both `voltage_type` and `type` field names (backward compatible)
- Support for both `voltage_amplitude` and `amplitude` field names
- Support for both `dielectric_permittivity` and `permittivity` field names
- Support for both `dielectric_thickness` and `thickness` field names
- Support for both `is_dielectric` and `dielectric` field names

#### Species Section:
- Added parsing for `mobility_coeff` (numeric coefficient)

#### Reactions Section:
- Parse `name_for_output` with fallback to empty string
- Parse `energy_change` with proper handling of empty strings
- Parse equation-based rate parameters:
  - `equation_constants`
  - `equation_values`
  - `equation_variables`
  - `rate_equation` (with backward compatibility for `equation` field)

#### Initial Conditions Section:
- Full implementation of initial condition parsing
- Supports multiple species and electron energy initialization

#### Solver Section:
- Parse all solver parameters with sensible defaults:
  - `type` (default: "JFNK")
  - `tolerance` (default: 1.0e-6)
  - `max_iterations` (default: 50)
  - `time_step` (default: 1.0e-12 s)
  - `end_time` (default: 1.0e-9 s)
  - `preconditioner` (default: "PBP")
  - `ksp_type` (default: "GMRES")

#### Output Section:
- Parse `minimum_time_interval`
- Parse `save_fields` array
- Parse `rates` array
- Parse `add_timestep` (supports both boolean and "yes"/"no" strings)
- Parse checkpoint configuration

#### Restart Section:
- Parse `enabled`, `file`, and `step` parameters

#### Advanced Section:
- Parse all advanced parameters:
  - `adaptive_dt`
  - `dt_min` and `dt_max`
  - `mpi_decomposition`
  - `density_floor`
  - `energy_floor`
  - `wall_quenching_probability`
  - `thermal_velocity_factor`

### 3. Updated main.cpp

Modified main.cpp to use configuration parameters instead of hardcoded values:

- Use `config.solver.time_step` instead of hardcoded `dt = 1e-12`
- Use `config.solver.end_time` instead of hardcoded `t_end = 1e-9`
- Support for command-line overrides via `-dt` and `-tend` flags
- Enhanced restart logic that respects both config file and command-line options
- Added informational output showing solver configuration

## Configuration File Parameter Support

### ✅ Fully Supported Sections

1. **mesh** - All parameters supported
2. **electrodes** - All parameters supported (with field name flexibility)
3. **plasma** - All parameters supported
4. **species** - All parameters supported (including new `mobility_coeff`)
5. **reactions** - All parameters supported (including equation-based rates)
6. **initial_condition** - Fully implemented
7. **solver** - Fully implemented
8. **output** - All parameters supported
9. **restart** - Fully implemented
10. **advanced** - Fully implemented

### Field Name Compatibility

The parser now supports multiple field name variations for backward compatibility:

| User's Config Field | Alternative Field Name | Status |
|---------------------|------------------------|--------|
| `voltage_type` | `type` | ✅ Both supported |
| `voltage_amplitude` | `amplitude` | ✅ Both supported |
| `dielectric_permittivity` | `permittivity` | ✅ Both supported |
| `dielectric_thickness` | `thickness` | ✅ Both supported |
| `is_dielectric` | `dielectric` | ✅ Both supported |

## Example Usage

### Basic Configuration
```yaml
solver:
  type: JFNK
  tolerance: 1.0e-8
  max_iterations: 200
  time_step: 1.8436578e-11
  end_time: 7.37463e-5
  preconditioner: PBP
  ksp_type: GMRES
```

### Initial Conditions
```yaml
initial_condition:
  - name: e
    value: 1e14  # 1/m^3
  - name: e_energy
    value: 1.5  # eV
  - name: Ar+
    value: 1e14
```

### Equation-Based Reactions
```yaml
reactions:
  - equation: "Ar1s3 -> Ar1s5"
    name_for_output: radia1
    type: radiation
    energy_change: 0.0
    rate_type: equation
    equation_constants: "rho_L M_gas Tgas g1s5"
    equation_values: "0.03 39.95 300.0 5.0"
    equation_variables: "Ar1s5_density Ar1s4_density"
    rate_equation: "{18900000 * (2.0 - exp(...))}"
```

### Advanced Configuration
```yaml
advanced:
  adaptive_dt: false
  dt_min: 1.0e-13
  dt_max: 1.0e-10
  density_floor: 1.0e6
  energy_floor: 0.1
  wall_quenching_probability: 0.1
  thermal_velocity_factor: 0.25
```

## Default Values

All configuration parameters have sensible default values, so minimal configurations will still work. Default values are set both in the struct definitions (as member initializers) and in the parsing logic (as fallback values).

## Benefits

1. **Complete Parameter Coverage** - All parameters from the specification are now managed
2. **Backward Compatibility** - Multiple field name variations supported
3. **Robust Defaults** - Sensible defaults for all optional parameters
4. **Type Safety** - Proper C++ types with initialization
5. **Flexibility** - Command-line overrides still supported for key parameters
6. **Maintainability** - Clear structure and documentation

## Files Modified

1. `/workspace/HydroPlas/src/config/ConfigParser.hpp` - Updated struct definitions
2. `/workspace/HydroPlas/src/config/ConfigParser.cpp` - Implemented parsing logic
3. `/workspace/HydroPlas/src/main.cpp` - Use config parameters instead of hardcoded values

## Verification

- ✅ All struct fields have default initializers
- ✅ All parsing logic handles missing optional parameters
- ✅ Field name variations supported for backward compatibility
- ✅ Syntax verified with g++ compiler
- ✅ Default values match user's expected ranges

## Next Steps (Optional Enhancements)

1. Add validation for parameter ranges (e.g., time_step > 0)
2. Add error messages for required but missing parameters
3. Implement config file versioning for future compatibility
4. Add unit tests for ConfigParser
5. Update documentation with all supported parameters

## Conclusion

The HydroPlas configuration management system now supports all parameters specified in the user's configuration file. The implementation is robust, flexible, and maintains backward compatibility with existing configuration files.
