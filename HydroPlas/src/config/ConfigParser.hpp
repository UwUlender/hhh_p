#pragma once

#include <string>
#include <vector>
#include <map>
#include <yaml-cpp/yaml.h>

namespace HydroPlas {

struct MeshConfig {
    std::string type; // "1D_rectilinear" or "2D_rectilinear"
    std::vector<double> x_nodes; // For non-uniform
    std::vector<double> y_nodes; // For non-uniform (if 2D)
    // Supports uniform definition via size/length if nodes are empty (logic in parser)
};

struct ElectrodeConfig {
    std::string name;
    std::string location; // "x_min", "x_max", "y_min", "y_max" or specific coords
    
    // Voltage definition
    std::string voltage_type; // "DC", "RF", "PULSE", "Expression"
    double voltage_amplitude = 0.0;
    double frequency = 0.0;
    double phase = 0.0;
    double bias = 0.0;
    double duty_cycle = 0.5;
    
    std::string voltage_expression; // For "Expression" type
    
    // Dielectric
    bool is_dielectric = false;
    double dielectric_permittivity = 1.0;
    double dielectric_thickness = 0.0;
    
    // Species BC
    double gamma_see = 0.0;
};

struct PlasmaConfig {
    double gas_pressure; // Pa
    double background_temp; // K
    std::string background_gas; // e.g., "Ar"
};

struct SpeciesConfig {
    std::string name;
    std::string type; // "electron", "ion", "neutral"
    double charge = 0.0; // In elementary charges (e.g., -1, 1, 0)
    double mass = 0.0; // kg
    double diffusion_coeff = 0.0; // For neutrals (m^2/s)
    double mobility_coeff = 0.0; // Numeric mobility coefficient (alternative to file)
    std::string mobility_file; // For charged (BOLSIG+)
};

struct ReactionConfig {
    std::string equation; // "e + Ar -> 2e + Ar+"
    std::string name_for_output; // Friendly name for output
    std::string type; // "ionization", "excitation", "elastic", "radiation", etc.
    double energy_change = 0.0; // Energy loss/gain in eV (for electron energy balance)
    std::string rate_type; // "constant", "arrhenius", "table", "equation"
    double a = 0.0; // For Arrhenius or Constant
    double b = 0.0; // For Arrhenius
    double e_a = 0.0; // For Arrhenius (Activation Energy)
    std::string table_file; // If table
    
    // For equation-based rates
    std::string equation_constants; // Space-separated names
    std::string equation_values; // Space-separated values
    std::string equation_variables; // Space-separated variable names
    std::string rate_equation; // Mathematical equation string for rate calculation
};

struct InitialConditionConfig {
    std::string name; // Species name or "e_energy"
    double value = 0.0; // Initial density (1/m³) or energy (eV)
};

struct SolverConfig {
    std::string type = "JFNK"; // "JFNK", etc.
    double tolerance = 1.0e-6; // SNES tolerance
    int max_iterations = 50; // Maximum solver iterations
    double time_step = 1.0e-12; // Time step size (seconds)
    double end_time = 1.0e-9; // Simulation end time (seconds)
    std::string preconditioner = "pbjacobi"; // Valid PETSc PC types: "pbjacobi", "bjacobi", "asm", "ilu", "none"
    std::string ksp_type = "gmres"; // Valid PETSc KSP types: "gmres", "cg", "bcgs", "fgmres"
};

struct OutputConfig {
    std::string format = "hdf5"; // "hdf5"
    int frequency_step = 100; // Output every N steps
    double minimum_time_interval = 0.0; // Minimum time between outputs (alternative)
    bool save_rates = false;
    std::string filename = "hydroplas_out.h5";
    std::vector<std::string> save_fields; // List of fields to save
    std::vector<std::string> rates; // List of reaction rates to save
    bool add_timestep = false; // Add timestep number to filename
    
    // Checkpoint configuration
    int checkpoint_frequency = 1000; // Save checkpoint every N steps
    std::string checkpoint_filename = "checkpoint.h5";
};

struct RestartConfig {
    bool enabled = false;
    std::string file = "checkpoint.h5";
    int step = 0;
};

struct AdvancedConfig {
    // Adaptive time stepping
    bool adaptive_dt = false;
    double dt_min = 1.0e-15;
    double dt_max = 1.0e-9;
    
    // Parallel configuration
    std::string mpi_decomposition = "auto"; // "auto" or specific decomposition
    
    // Convergence criteria
    double density_floor = 1.0e6; // Minimum density (1/m³)
    double energy_floor = 0.1; // Minimum electron energy (eV)
    
    // Boundary condition details
    double wall_quenching_probability = 0.0; // For metastables
    double thermal_velocity_factor = 0.25; // For thermal flux BCs
};

struct SimulationConfig {
    MeshConfig mesh;
    std::vector<ElectrodeConfig> electrodes;
    PlasmaConfig plasma;
    std::vector<SpeciesConfig> species;
    std::vector<ReactionConfig> reactions;
    std::vector<InitialConditionConfig> initial_conditions;
    SolverConfig solver;
    OutputConfig output;
    RestartConfig restart;
    AdvancedConfig advanced;
    
    // Helper to store the full YAML text for saving
    std::string full_yaml_content;
};

class ConfigParser {
public:
    explicit ConfigParser(const std::string& config_file_path);
    SimulationConfig get_config() const;

private:
    SimulationConfig config_;
    void parse_yaml(const YAML::Node& root);
};

} // namespace HydroPlas
