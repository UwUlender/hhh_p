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
        
        // Check if using new multi-electrode format
        if (b.contains("electrodes") && b["electrodes"].is_array()) {
            config_.boundary.use_multi_electrode = true;
            
            for (const auto& elec : b["electrodes"]) {
                ElectrodeConfig ec;
                ec.name = elec.value("name", "unnamed");
                ec.voltage_type = elec.value("voltage_type", "DC");
                ec.voltage_amplitude = elec.value("voltage_amplitude", 0.0);
                ec.frequency = elec.value("frequency", 13.56e6);
                ec.bias = elec.value("bias", 0.0);
                ec.phase = elec.value("phase", 0.0);
                ec.duty_cycle = elec.value("duty_cycle", 0.5);
                ec.gamma_see = elec.value("gamma_see", 0.1);
                ec.is_dielectric = elec.value("is_dielectric", false);
                ec.dielectric_permittivity = elec.value("dielectric_permittivity", 1.0);
                ec.dielectric_thickness = elec.value("dielectric_thickness", 0.0);
                
                config_.boundary.electrodes.push_back(ec);
            }
        } else {
            // Legacy single-electrode format
            config_.boundary.use_multi_electrode = false;
            config_.boundary.voltage_type = b.value("voltage_type", "DC");
            config_.boundary.voltage_amplitude = b.value("voltage_amplitude", 300.0);
            config_.boundary.frequency = b.value("frequency", 13.56e6);
            config_.boundary.bias = b.value("bias", 0.0);
            config_.boundary.gamma_see = b.value("gamma_see", 0.1);
            config_.boundary.dielectric_permittivity = b.value("dielectric_permittivity", 1.0);
            config_.boundary.dielectric_thickness = b.value("dielectric_thickness", 1e-3);
        }
    }

    // Chemistry
    if (j.contains("chemistry")) {
        auto& c = j["chemistry"];
        config_.chemistry.mode = c.value("mode", "Pre-calculated");
        config_.chemistry.cross_section_file = c.value("cross_section_file", "data/cs.dat");
        config_.chemistry.transport_table_file = c.value("transport_table_file", "data/transport.dat");
        config_.chemistry.gas_velocity = c.value("gas_velocity", 0.0);
        config_.chemistry.gas_temperature = c.value("gas_temperature", 300.0);

        if (c.contains("excited_species") && c["excited_species"].is_array()) {
            for (const auto& es : c["excited_species"]) {
                ExcitedSpeciesConfig esc;
                esc.name = es.value("name", "Unknown");
                esc.diffusion_coeff = es.value("diffusion_coeff", 1.5e-4);
                esc.mass = es.value("mass", 6.63e-26);
                esc.energy_level = es.value("energy_level", 11.5);
                esc.wall_quenching_prob = es.value("wall_quenching_prob", 1.0);
                esc.wall_see_prob = es.value("wall_see_prob", 0.0);
                config_.chemistry.excited_species.push_back(esc);
            }
        }

        if (c.contains("reactions") && c["reactions"].is_array()) {
            for (const auto& r : c["reactions"]) {
                ReactionConfig rc;
                rc.type = r.value("type", "unknown");
                rc.species = r.value("species", "");
                rc.rate_source = r.value("rate_source", "constant");
                rc.table_col = r.value("table_col", "");
                rc.rate_constant = r.value("rate_constant", 0.0);
                if (r.contains("reactants")) rc.reactants = r["reactants"].get<std::vector<std::string>>();
                if (r.contains("products")) rc.products = r["products"].get<std::vector<std::string>>();
                config_.chemistry.reactions.push_back(rc);
            }
        }
    }
}

} // namespace HydroPlas
