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

struct BoundaryConfig {
    std::string voltage_type; // "DC", "RF"
    double voltage_amplitude; // V_rf or V_dc
    double frequency;         // For RF
    double bias;              // V_bias
    double gamma_see;         // Secondary Electron Emission coefficient
    double dielectric_permittivity;
    double dielectric_thickness;
};

struct ChemistryConfig {
    std::string mode; // "Pre-calculated" or "Inline BOLSIG+"
    std::string cross_section_file;
    std::string transport_table_file;
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
