#include "ConfigParser.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace HydroPlas {

using json = nlohmann::json;

ConfigParser::ConfigParser(const std::string& config_file_path) {
    parse_file(config_file_path);
}

SimulationConfig ConfigParser::get_config() const {
    return config_;
}

void ConfigParser::parse_file(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open configuration file: " + path);
    }

    json j;
    try {
        file >> j;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("JSON parse error: " + std::string(e.what()));
    }

    // Domain
    if (j.contains("domain")) {
        auto& d = j["domain"];
        config_.domain.Lx = d.value("Lx", 0.01); // Default 1cm
        config_.domain.Ly = d.value("Ly", 0.0);  // 0 means 1D
        config_.domain.Nx = d.value("Nx", 100);
        config_.domain.Ny = d.value("Ny", 1);
    }

    // Time
    if (j.contains("time")) {
        auto& t = j["time"];
        config_.time.dt = t.value("dt", 1e-12);
        config_.time.t_end = t.value("t_end", 1e-6);
        config_.time.output_interval = t.value("output_interval", 100);
    }

    // Boundary
    if (j.contains("boundary")) {
        auto& b = j["boundary"];
        config_.boundary.voltage_type = b.value("voltage_type", "DC");
        config_.boundary.voltage_amplitude = b.value("voltage_amplitude", 300.0);
        config_.boundary.frequency = b.value("frequency", 13.56e6);
        config_.boundary.bias = b.value("bias", 0.0);
        config_.boundary.gamma_see = b.value("gamma_see", 0.1);
        config_.boundary.dielectric_permittivity = b.value("dielectric_permittivity", 1.0);
        config_.boundary.dielectric_thickness = b.value("dielectric_thickness", 1e-3);
    }

    // Chemistry
    if (j.contains("chemistry")) {
        auto& c = j["chemistry"];
        config_.chemistry.mode = c.value("mode", "Pre-calculated");
        config_.chemistry.cross_section_file = c.value("cross_section_file", "data/cs.dat");
        config_.chemistry.transport_table_file = c.value("transport_table_file", "data/transport.dat");
    }
}

} // namespace HydroPlas
