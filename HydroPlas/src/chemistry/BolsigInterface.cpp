#include "BolsigInterface.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>

namespace HydroPlas {

/**
 * @brief Create BOLSIG+ input file from cross-section data
 * 
 * BOLSIG+ expects a specific input format:
 * - Cross-section file path
 * - Gas conditions (pressure, temperature)
 * - E/N range or energy range
 */
void BolsigInterface::create_bolsig_input(const ChemistryConfig& config, 
                                         const std::string& input_file) {
    std::ofstream inp(input_file);
    if (!inp.is_open()) {
        throw std::runtime_error("Cannot create BOLSIG+ input file: " + input_file);
    }
    
    // BOLSIG+ input format (simplified)
    inp << "/ Cross section file" << std::endl;
    inp << config.cross_section_file << std::endl;
    inp << std::endl;
    inp << "/ Gas temperature [K]" << std::endl;
    inp << config.gas_temperature << std::endl;
    inp << std::endl;
    inp << "/ Energy range [eV] and number of points" << std::endl;
    inp << "0.01 100.0 100" << std::endl;  // E_min E_max N_points
    inp << std::endl;
    inp << "/ Output file" << std::endl;
    inp << "bolsig_output.dat" << std::endl;
    
    inp.close();
    std::cout << "Created BOLSIG+ input file: " << input_file << std::endl;
}

/**
 * @brief Parse BOLSIG+ output file and convert to lookup table format
 * 
 * BOLSIG+ output typically contains:
 * - Mean energy [eV]
 * - Mobility [m2/Vs]
 * - Diffusion coefficient [m2/s]
 * - Rate coefficients [m3/s]
 * 
 * Format may vary depending on BOLSIG+ version. This parser handles common formats.
 */
bool BolsigInterface::parse_bolsig_output(const std::string& bolsig_output,
                                          const std::string& table_output) {
    std::ifstream inp(bolsig_output);
    if (!inp.is_open()) {
        std::cerr << "Cannot open BOLSIG+ output: " << bolsig_output << std::endl;
        return false;
    }
    
    std::ofstream out(table_output);
    if (!out.is_open()) {
        std::cerr << "Cannot create table output: " << table_output << std::endl;
        return false;
    }
    
    // Write header
    out << "# HydroPlas transport table generated from BOLSIG+" << std::endl;
    out << "# Energy[eV] Mobility_e[m2/Vs] Diff_e[m2/s] Mobility_i[m2/Vs] Diff_i[m2/s] "
        << "Rate_Ionization[m3/s] Rate_Excitation[m3/s]" << std::endl;
    
    std::string line;
    std::vector<double> energies, mobilities, diffusions;
    std::vector<double> rate_ion, rate_exc;
    
    // Parse BOLSIG+ output (format-specific)
    // Common format: columns of E/N, mean energy, mobility, diffusion, rates...
    while (std::getline(inp, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '/') continue;
        
        std::stringstream ss(line);
        double e_n, mean_e, mu, D, k_iz, k_exc;
        
        // Try to parse common BOLSIG+ format
        if (ss >> e_n >> mean_e >> mu >> D) {
            // Successfully parsed transport coefficients
            energies.push_back(mean_e);
            mobilities.push_back(mu);
            diffusions.push_back(D);
            
            // Try to read rate coefficients (may not always be present)
            if (ss >> k_iz) {
                rate_ion.push_back(k_iz);
            } else {
                rate_ion.push_back(0.0);
            }
            
            if (ss >> k_exc) {
                rate_exc.push_back(k_exc);
            } else {
                rate_exc.push_back(0.0);
            }
        }
    }
    
    inp.close();
    
    if (energies.empty()) {
        std::cerr << "No data parsed from BOLSIG+ output!" << std::endl;
        out.close();
        return false;
    }
    
    // Write parsed data in HydroPlas format
    double mu_i = 1.0e-4;  // Typical ion mobility in Ar [m2/Vs]
    double D_i = 0.026;    // Typical ion diffusion in Ar [m2/s]
    
    out << std::scientific << std::setprecision(6);
    for (size_t i = 0; i < energies.size(); ++i) {
        out << energies[i] << " "
            << mobilities[i] << " "
            << diffusions[i] << " "
            << mu_i << " "
            << D_i << " "
            << rate_ion[i] << " "
            << rate_exc[i] << std::endl;
    }
    
    out.close();
    std::cout << "Parsed " << energies.size() << " data points from BOLSIG+ output" << std::endl;
    return true;
}

void BolsigInterface::run_bolsig(const ChemistryConfig& config, const std::string& output_file) {
    std::cout << "========================================" << std::endl;
    std::cout << "BOLSIG+ Interface: Inline Chemistry Mode" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Check if cross-section file exists
    std::ifstream cs_file(config.cross_section_file);
    if (!cs_file.good()) {
        std::cerr << "WARNING: Cross-section file not found: " << config.cross_section_file << std::endl;
        std::cerr << "Generating fallback transport data..." << std::endl;
        generate_fallback_data(output_file);
        return;
    }
    cs_file.close();
    
    // Step 1: Create BOLSIG+ input file
    std::string input_file = "bolsig_input.txt";
    try {
        create_bolsig_input(config, input_file);
    } catch (const std::exception& e) {
        std::cerr << "Error creating BOLSIG+ input: " << e.what() << std::endl;
        generate_fallback_data(output_file);
        return;
    }
    
    // Step 2: Execute BOLSIG+
    // Try multiple common BOLSIG+ executable names
    std::vector<std::string> bolsig_names = {"bolsigplus", "bolsig+", "bolsig", "BOLSIG+"};
    std::string cmd;
    int ret = -1;
    
    for (const auto& name : bolsig_names) {
        cmd = name + " " + input_file + " > bolsig_raw_output.dat 2>&1";
        std::cout << "Attempting: " << cmd << std::endl;
        ret = std::system(cmd.c_str());
        if (ret == 0) {
            std::cout << "BOLSIG+ execution successful!" << std::endl;
            break;
        }
    }
    
    if (ret != 0) {
        std::cerr << "WARNING: BOLSIG+ execution failed (executable not found or error)" << std::endl;
        std::cerr << "Searched for: ";
        for (const auto& name : bolsig_names) std::cerr << name << " ";
        std::cerr << std::endl;
        std::cerr << "Generating fallback transport data..." << std::endl;
        generate_fallback_data(output_file);
        return;
    }
    
    // Step 3: Parse BOLSIG+ output
    bool success = parse_bolsig_output("bolsig_output.dat", output_file);
    if (!success) {
        std::cerr << "WARNING: Failed to parse BOLSIG+ output" << std::endl;
        std::cerr << "Generating fallback transport data..." << std::endl;
        generate_fallback_data(output_file);
    } else {
        std::cout << "Transport table successfully generated: " << output_file << std::endl;
    }
    
    std::cout << "========================================" << std::endl;
}

void BolsigInterface::generate_fallback_data(const std::string& output_file) {
    std::cout << "Generating fallback Argon transport data (analytical approximations)..." << std::endl;
    
    std::ofstream dummy(output_file);
    if (!dummy.is_open()) {
        throw std::runtime_error("Cannot create fallback data file: " + output_file);
    }
    
    // Write header
    dummy << "# Fallback Argon transport data (analytical approximations)" << std::endl;
    dummy << "# Energy[eV] Mobility_e[m2/Vs] Diff_e[m2/s] Mobility_i[m2/Vs] Diff_i[m2/s] "
         << "Rate_Ionization[m3/s] Rate_Excitation[m3/s]" << std::endl;
    
    dummy << std::scientific << std::setprecision(6);
    
    // Generate realistic Argon data using empirical formulas
    std::vector<double> energies = {0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 
                                   15.0, 20.0, 30.0, 50.0, 100.0};
    
    for (double E : energies) {
        // Electron mobility: mu_e ~ 1 / (1 + E^0.5) [rough approximation]
        double mu_e = 100.0 / (1.0 + std::sqrt(E));
        
        // Diffusion coefficient: D_e ~ mu_e * E (Einstein relation)
        double D_e = mu_e * E * 1.6e-19 / (1.38e-23 * 11600.0);  // [m2/s]
        
        // Ion transport (weakly energy-dependent)
        double mu_i = 1.5e-4;  // [m2/Vs]
        double D_i = 0.026;    // [m2/s]
        
        // Ionization rate: k_iz ~ exp(-E_iz / E) for E > E_threshold
        double E_iz = 15.76;  // Argon ionization potential [eV]
        double k_iz = 0.0;
        if (E > 4.0) {
            k_iz = 1e-13 * std::exp(-E_iz / E);
        }
        
        // Excitation rate: k_exc ~ exp(-E_exc / E)
        double E_exc = 11.55;  // Argon metastable [eV]
        double k_exc = 0.0;
        if (E > 2.0) {
            k_exc = 5e-14 * std::exp(-E_exc / E);
        }
        
        dummy << E << " " << mu_e << " " << D_e << " " << mu_i << " " << D_i << " "
             << k_iz << " " << k_exc << std::endl;
    }
    
    dummy.close();
    std::cout << "Fallback data generated with " << energies.size() << " energy points" << std::endl;
}

} // namespace HydroPlas
