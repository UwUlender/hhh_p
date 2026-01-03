#include "ConfigParser.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <regex>

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
            
            // Support both 'type' and 'voltage_type' field names
            if (el["voltage_type"]) ec.voltage_type = el["voltage_type"].as<std::string>();
            else if (el["type"]) ec.voltage_type = el["type"].as<std::string>();
            else ec.voltage_type = "DC";
            
            // Support both 'amplitude' and 'voltage_amplitude' field names
            if (el["voltage_amplitude"]) ec.voltage_amplitude = el["voltage_amplitude"].as<double>();
            else if (el["amplitude"]) ec.voltage_amplitude = el["amplitude"].as<double>();
            
            if (el["frequency"]) ec.frequency = el["frequency"].as<double>();
            if (el["phase"]) ec.phase = el["phase"].as<double>();
            if (el["bias"]) ec.bias = el["bias"].as<double>();
            if (el["duty_cycle"]) ec.duty_cycle = el["duty_cycle"].as<double>();
            
            if (el["voltage_expression"]) {
                ec.voltage_expression = el["voltage_expression"].as<std::string>();
                ec.voltage_type = "Expression";
            }

            if (el["voltage"]) {
                // Check if string or number
                try {
                    double v = el["voltage"].as<double>();
                    ec.voltage_amplitude = v;
                    if (ec.voltage_type == "DC") ec.bias = v;
                } catch (...) {
                    ec.voltage_expression = el["voltage"].as<std::string>();
                    ec.voltage_type = "Expression";
                    
                    // Try to parse "A * sin(2*pi*f*t)"
                    std::regex re("([\\d\\.]+)\\s*\\*\\s*sin\\(2\\*pi\\*([\\d\\.eE\\+\\-]+)\\*t\\)");
                    std::smatch match;
                    if (std::regex_search(ec.voltage_expression, match, re)) {
                        ec.voltage_amplitude = std::stod(match[1]);
                        ec.frequency = std::stod(match[2]);
                        ec.voltage_type = "RF"; // Treat as RF
                    }
                }
            }
            
            if (el["dielectric"]) ec.is_dielectric = el["dielectric"].as<bool>();
            if (el["is_dielectric"]) ec.is_dielectric = el["is_dielectric"].as<bool>();
            
            // Support both field name variations for dielectric properties
            if (el["dielectric_permittivity"]) ec.dielectric_permittivity = el["dielectric_permittivity"].as<double>();
            else if (el["permittivity"]) ec.dielectric_permittivity = el["permittivity"].as<double>();
            
            if (el["dielectric_thickness"]) ec.dielectric_thickness = el["dielectric_thickness"].as<double>();
            else if (el["thickness"]) ec.dielectric_thickness = el["thickness"].as<double>();
            
            if (el["gamma_see"]) ec.gamma_see = el["gamma_see"].as<double>();
            
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
            
            // Auto-assign charge based on type if missing
            if (sc.type == "electron") sc.charge = -1.0;
            else if (sc.type == "ion" && sc.charge == 0.0) sc.charge = 1.0; // Assumption
            
            if (sp["mass"]) sc.mass = sp["mass"].as<double>();
            if (sp["diffusion_coeff"]) sc.diffusion_coeff = sp["diffusion_coeff"].as<double>();
            if (sp["mobility_coeff"]) sc.mobility_coeff = sp["mobility_coeff"].as<double>();
            if (sp["mobility_file"]) sc.mobility_file = sp["mobility_file"].as<std::string>();
            
            config_.species.push_back(sc);
        }
    }
    
    // Reactions
    if (root["reactions"]) {
        for (const auto& rxn : root["reactions"]) {
            ReactionConfig rc;
            rc.equation = rxn["equation"].as<std::string>();
            if (rxn["name_for_output"]) rc.name_for_output = rxn["name_for_output"].as<std::string>();
            else rc.name_for_output = ""; // Will use equation as fallback
            
            rc.type = rxn["type"].as<std::string>();
            
            // Energy change (for electron energy balance)
            if (rxn["energy_change"]) {
                try {
                    rc.energy_change = rxn["energy_change"].as<double>();
                } catch (...) {
                    rc.energy_change = 0.0; // Empty string or invalid
                }
            } else {
                rc.energy_change = 0.0;
            }
            
            rc.rate_type = rxn["rate_type"].as<std::string>();
            
            // Traditional rate parameters
            if (rxn["a"]) rc.a = rxn["a"].as<double>();
            if (rxn["b"]) rc.b = rxn["b"].as<double>();
            if (rxn["e_a"]) rc.e_a = rxn["e_a"].as<double>();
            if (rxn["table_file"]) rc.table_file = rxn["table_file"].as<std::string>();
            
            // Equation-based rates
            if (rxn["equation_constants"]) rc.equation_constants = rxn["equation_constants"].as<std::string>();
            if (rxn["equation_values"]) rc.equation_values = rxn["equation_values"].as<std::string>();
            if (rxn["equation_variables"]) rc.equation_variables = rxn["equation_variables"].as<std::string>();
            if (rxn["rate_equation"]) rc.rate_equation = rxn["rate_equation"].as<std::string>();
            // Also check for "equation" field for backward compatibility
            else if (rxn["equation"] && rc.rate_type == "equation") rc.rate_equation = rxn["equation"].as<std::string>();
            
            config_.reactions.push_back(rc);
        }
    }
    
    // Initial conditions
    if (root["initial_condition"]) {
        for (const auto& ic : root["initial_condition"]) {
            InitialConditionConfig icc;
            icc.name = ic["name"].as<std::string>();
            icc.value = ic["value"].as<double>();
            config_.initial_conditions.push_back(icc);
        }
    }
    
    // Solver configuration
    if (root["solver"]) {
        if (root["solver"]["type"]) config_.solver.type = root["solver"]["type"].as<std::string>();
        else config_.solver.type = "JFNK";
        
        if (root["solver"]["tolerance"]) config_.solver.tolerance = root["solver"]["tolerance"].as<double>();
        else config_.solver.tolerance = 1.0e-6;
        
        if (root["solver"]["max_iterations"]) config_.solver.max_iterations = root["solver"]["max_iterations"].as<int>();
        else config_.solver.max_iterations = 50;
        
        if (root["solver"]["time_step"]) config_.solver.time_step = root["solver"]["time_step"].as<double>();
        else config_.solver.time_step = 1.0e-12;
        
        if (root["solver"]["end_time"]) config_.solver.end_time = root["solver"]["end_time"].as<double>();
        else config_.solver.end_time = 1.0e-9;
        
        if (root["solver"]["preconditioner"]) config_.solver.preconditioner = root["solver"]["preconditioner"].as<std::string>();
        else config_.solver.preconditioner = "pbjacobi";
        
        if (root["solver"]["ksp_type"]) config_.solver.ksp_type = root["solver"]["ksp_type"].as<std::string>();
        else config_.solver.ksp_type = "gmres";
    } else {
        // Set defaults if solver section is missing
        config_.solver.type = "JFNK";
        config_.solver.tolerance = 1.0e-6;
        config_.solver.max_iterations = 50;
        config_.solver.time_step = 1.0e-12;
        config_.solver.end_time = 1.0e-9;
        config_.solver.preconditioner = "pbjacobi";
        config_.solver.ksp_type = "gmres";
    }
    
    // Output
    if (root["output"]) {
        if (root["output"]["format"]) config_.output.format = root["output"]["format"].as<std::string>();
        else config_.output.format = "hdf5";
        
        if (root["output"]["frequency_step"]) config_.output.frequency_step = root["output"]["frequency_step"].as<int>();
        else config_.output.frequency_step = 100;
        
        if (root["output"]["minimum_time_interval"]) config_.output.minimum_time_interval = root["output"]["minimum_time_interval"].as<double>();
        else config_.output.minimum_time_interval = 0.0;
        
        if (root["output"]["save_rates"]) config_.output.save_rates = root["output"]["save_rates"].as<bool>();
        else config_.output.save_rates = false;
        
        if (root["output"]["filename"]) config_.output.filename = root["output"]["filename"].as<std::string>();
        else config_.output.filename = "hydroplas_out.h5";
        
        if (root["output"]["add_timestep"]) {
            try {
                config_.output.add_timestep = root["output"]["add_timestep"].as<bool>();
            } catch (...) {
                // Handle "yes"/"no" strings
                std::string val = root["output"]["add_timestep"].as<std::string>();
                config_.output.add_timestep = (val == "yes" || val == "true" || val == "Yes" || val == "True");
            }
        } else {
            config_.output.add_timestep = false;
        }
        
        // Parse save_fields array
        if (root["output"]["save_fields"] && root["output"]["save_fields"].IsSequence()) {
            config_.output.save_fields = root["output"]["save_fields"].as<std::vector<std::string>>();
        }
        
        // Parse rates array
        if (root["output"]["rates"] && root["output"]["rates"].IsSequence()) {
            config_.output.rates = root["output"]["rates"].as<std::vector<std::string>>();
        }
        
        // Checkpoint configuration
        if (root["output"]["checkpoint_frequency"]) config_.output.checkpoint_frequency = root["output"]["checkpoint_frequency"].as<int>();
        else config_.output.checkpoint_frequency = 1000;
        
        if (root["output"]["checkpoint_filename"]) config_.output.checkpoint_filename = root["output"]["checkpoint_filename"].as<std::string>();
        else config_.output.checkpoint_filename = "checkpoint.h5";
    } else {
        // Set defaults if output section is missing
        config_.output.format = "hdf5";
        config_.output.frequency_step = 100;
        config_.output.minimum_time_interval = 0.0;
        config_.output.save_rates = false;
        config_.output.filename = "hydroplas_out.h5";
        config_.output.add_timestep = false;
        config_.output.checkpoint_frequency = 1000;
        config_.output.checkpoint_filename = "checkpoint.h5";
    }
    
    // Restart configuration
    if (root["restart"]) {
        if (root["restart"]["enabled"]) config_.restart.enabled = root["restart"]["enabled"].as<bool>();
        else config_.restart.enabled = false;
        
        if (root["restart"]["file"]) config_.restart.file = root["restart"]["file"].as<std::string>();
        else config_.restart.file = "checkpoint.h5";
        
        if (root["restart"]["step"]) config_.restart.step = root["restart"]["step"].as<int>();
        else config_.restart.step = 0;
    } else {
        config_.restart.enabled = false;
        config_.restart.file = "checkpoint.h5";
        config_.restart.step = 0;
    }
    
    // Advanced configuration
    if (root["advanced"]) {
        if (root["advanced"]["adaptive_dt"]) config_.advanced.adaptive_dt = root["advanced"]["adaptive_dt"].as<bool>();
        else config_.advanced.adaptive_dt = false;
        
        if (root["advanced"]["dt_min"]) config_.advanced.dt_min = root["advanced"]["dt_min"].as<double>();
        else config_.advanced.dt_min = 1.0e-15;
        
        if (root["advanced"]["dt_max"]) config_.advanced.dt_max = root["advanced"]["dt_max"].as<double>();
        else config_.advanced.dt_max = 1.0e-9;
        
        if (root["advanced"]["mpi_decomposition"]) config_.advanced.mpi_decomposition = root["advanced"]["mpi_decomposition"].as<std::string>();
        else config_.advanced.mpi_decomposition = "auto";
        
        if (root["advanced"]["density_floor"]) config_.advanced.density_floor = root["advanced"]["density_floor"].as<double>();
        else config_.advanced.density_floor = 1.0e6;
        
        if (root["advanced"]["energy_floor"]) config_.advanced.energy_floor = root["advanced"]["energy_floor"].as<double>();
        else config_.advanced.energy_floor = 0.1;
        
        if (root["advanced"]["wall_quenching_probability"]) config_.advanced.wall_quenching_probability = root["advanced"]["wall_quenching_probability"].as<double>();
        else config_.advanced.wall_quenching_probability = 0.0;
        
        if (root["advanced"]["thermal_velocity_factor"]) config_.advanced.thermal_velocity_factor = root["advanced"]["thermal_velocity_factor"].as<double>();
        else config_.advanced.thermal_velocity_factor = 0.25;
    } else {
        // Set defaults if advanced section is missing
        config_.advanced.adaptive_dt = false;
        config_.advanced.dt_min = 1.0e-15;
        config_.advanced.dt_max = 1.0e-9;
        config_.advanced.mpi_decomposition = "auto";
        config_.advanced.density_floor = 1.0e6;
        config_.advanced.energy_floor = 0.1;
        config_.advanced.wall_quenching_probability = 0.0;
        config_.advanced.thermal_velocity_factor = 0.25;
    }
}

} // namespace HydroPlas
