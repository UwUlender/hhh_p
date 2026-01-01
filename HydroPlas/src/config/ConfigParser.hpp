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
    std::string voltage_expression; // "300 * sin(...)"
    double constant_voltage = 0.0;
};

struct PlasmaConfig {
    double gas_pressure; // Pa
    double background_temp; // K
    std::string background_gas; // e.g., "Ar"
};

struct SpeciesConfig {
    std::string name;
    std::string type; // "electron", "ion", "neutral_excited"
    double charge; // In elementary charges (e.g., -1, 1, 0)
    double mass; // kg
    double diffusion_coeff; // For neutrals (m^2/s)
    std::string mobility_file; // For charged (BOLSIG+)
};

struct ReactionConfig {
    std::string equation; // "e + Ar -> 2e + Ar+"
    std::string type; // "ionization", "excitation", etc.
    std::string rate_type; // "constant", "arrhenius", "table"
    double a; // For Arrhenius or Constant
    double b; // For Arrhenius
    double e_a; // For Arrhenius (Activation Energy)
    std::string table_file; // If table
};

struct OutputConfig {
    std::string format; // "hdf5"
    int frequency_step;
    bool save_rates;
    std::string filename;
};

struct SimulationConfig {
    MeshConfig mesh;
    std::vector<ElectrodeConfig> electrodes;
    PlasmaConfig plasma;
    std::vector<SpeciesConfig> species;
    std::vector<ReactionConfig> reactions;
    OutputConfig output;
    
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
