#pragma once

#include <string>
#include <vector>
#include <nlohmann/json.hpp>

namespace HydroPlas {

struct DomainConfig {
    double Lx;
    double Ly;
    int Nx;
    int Ny;
};

struct TimeConfig {
    double dt;
    double t_end;
    int output_interval;
};

struct ElectrodeConfig {
    std::string name;              // Electrode identifier (e.g., "left", "right", "top", "bottom")
    std::string voltage_type;      // "DC", "RF", "AC", "PULSE"
    double voltage_amplitude;      // V_amplitude (Volts)
    double frequency;              // For RF/AC (Hz)
    double bias;                   // DC bias voltage (Volts)
    double phase;                  // Phase offset (radians)
    double duty_cycle;             // For PULSE type (0-1)
    double gamma_see;              // Secondary Electron Emission coefficient
    bool is_dielectric;            // Whether this electrode has a dielectric layer
    double dielectric_permittivity; // Relative permittivity (ε_r)
    double dielectric_thickness;    // Thickness in meters
};

struct BoundaryConfig {
    // Legacy single-electrode support (for backward compatibility)
    std::string voltage_type; // "DC", "RF"
    double voltage_amplitude; // V_rf or V_dc
    double frequency;         // For RF
    double bias;              // V_bias
    double gamma_see;         // Secondary Electron Emission coefficient
    double dielectric_permittivity;
    double dielectric_thickness;
    
    // New multi-electrode support
    std::vector<ElectrodeConfig> electrodes;
    bool use_multi_electrode; // Flag to use new electrode system
};

struct ExcitedSpeciesConfig {
    std::string name;
    double diffusion_coeff;
    double mass;
    double energy_level;
    double wall_quenching_prob; // gamma_k
    double wall_see_prob;       // gamma_see_k (secondary electron emission)
};

struct ReactionConfig {
    std::string type; // "ionization", "excitation", "stepwise", "penning", "quenching"
    std::string species; // Target species name (e.g., "Ar_m")
    std::string rate_source; // "table" or "constant"
    std::string table_col; // If table
    double rate_constant; // If constant
    std::vector<std::string> reactants;
    std::vector<std::string> products;
};

struct ChemistryConfig {
    std::string mode; // "Pre-calculated" or "Inline BOLSIG+"
    std::string cross_section_file;
    std::string transport_table_file;
    double gas_velocity; // u_gas (m/s)
    double gas_temperature; // T_gas (K)
    std::vector<ExcitedSpeciesConfig> excited_species;
    std::vector<ReactionConfig> reactions;
};

struct SimulationConfig {
    DomainConfig domain;
    TimeConfig time;
    BoundaryConfig boundary;
    ChemistryConfig chemistry;
};

class ConfigParser {
public:
    explicit ConfigParser(const std::string& config_file_path);
    SimulationConfig get_config() const;

private:
    SimulationConfig config_;
    void parse_file(const std::string& path);
};

} // namespace HydroPlas
