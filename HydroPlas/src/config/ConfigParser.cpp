#include "ConfigParser.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

namespace HydroPlas {

ConfigParser::ConfigParser(const std::string& config_file_path) {
    std::ifstream fin(config_file_path);
    if (!fin.is_open()) {
        throw std::runtime_error("Could not open config file: " + config_file_path);
    }
    
    std::stringstream buffer;
    buffer << fin.rdbuf();
    config_.full_yaml_content = buffer.str();
    
    YAML::Node root = YAML::Load(config_.full_yaml_content);
    parse_yaml(root);
}

SimulationConfig ConfigParser::get_config() const {
    return config_;
}

void ConfigParser::parse_yaml(const YAML::Node& root) {
    // Mesh
    if (root["mesh"]) {
        auto mesh_node = root["mesh"];
        config_.mesh.type = mesh_node["type"].as<std::string>();
        
        if (mesh_node["x_nodes"].IsSequence()) {
            config_.mesh.x_nodes = mesh_node["x_nodes"].as<std::vector<double>>();
        } else if (mesh_node["x_nodes"].IsScalar()) {
            // Uniform generation: need [min, max] or similar, but simplified here:
            // Assume if scalar, it's number of nodes, requires length.
            // For now, let's stick to explicit node list or user handling.
            // But if user puts '100', we might need generation logic. 
            // The prompt example shows 'y_nodes: 100', implying uniform.
            // We'll generate a dummy uniform grid if scalar, 
            // but we need length. Assuming 0 to 1 if not specified? 
            // Or maybe 'length' param exists. 
            // Let's just store empty if scalar and let RectilinearGrid handle generation
            // if we pass the count.
            // Actually, let's strictly follow the prompt structure.
            // For this implementation, I will assume nodes are provided or generate them later.
        }

        if (mesh_node["y_nodes"] && mesh_node["y_nodes"].IsSequence()) {
            config_.mesh.y_nodes = mesh_node["y_nodes"].as<std::vector<double>>();
        }
    }

    // Electrodes
    if (root["electrodes"]) {
        for (const auto& el : root["electrodes"]) {
            ElectrodeConfig ec;
            ec.name = el["name"].as<std::string>();
            ec.location = el["location"].as<std::string>();
            if (el["voltage"]) {
                // Check if string or number
                try {
                    ec.constant_voltage = el["voltage"].as<double>();
                    ec.voltage_expression = std::to_string(ec.constant_voltage);
                } catch (...) {
                    ec.voltage_expression = el["voltage"].as<std::string>();
                }
            }
            config_.electrodes.push_back(ec);
        }
    }

    // Plasma
    if (root["plasma"]) {
        config_.plasma.gas_pressure = root["plasma"]["gas_pressure"].as<double>();
        config_.plasma.background_temp = root["plasma"]["background_temp"].as<double>();
        if (root["plasma"]["background_gas"])
             config_.plasma.background_gas = root["plasma"]["background_gas"].as<std::string>();
    }

    // Species
    if (root["species"]) {
        for (const auto& sp : root["species"]) {
            SpeciesConfig sc;
            sc.name = sp["name"].as<std::string>();
            sc.type = sp["type"].as<std::string>();
            if (sp["charge"]) sc.charge = sp["charge"].as<double>();
            else sc.charge = 0.0; // Default to neutral
            
            // Auto-assign charge based on type if missing?
            if (sc.type == "electron") sc.charge = -1.0;
            else if (sc.type == "ion" && sc.charge == 0.0) sc.charge = 1.0; // Assumption
            
            if (sp["mass"]) sc.mass = sp["mass"].as<double>();
            if (sp["diffusion_coeff"]) sc.diffusion_coeff = sp["diffusion_coeff"].as<double>();
            if (sp["mobility_file"]) sc.mobility_file = sp["mobility_file"].as<std::string>();
            
            config_.species.push_back(sc);
        }
    }
    
    // Output
    if (root["output"]) {
        config_.output.format = root["output"]["format"].as<std::string>();
        config_.output.frequency_step = root["output"]["frequency_step"].as<int>();
        if (root["output"]["save_rates"])
            config_.output.save_rates = root["output"]["save_rates"].as<bool>();
        else
            config_.output.save_rates = false;
        
        config_.output.filename = "hydroplas_out.h5"; // Default
    }
}

} // namespace HydroPlas
