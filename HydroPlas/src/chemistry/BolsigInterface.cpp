#include "BolsigInterface.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>

namespace HydroPlas {

void BolsigInterface::run_bolsig(const ChemistryConfig& config, const std::string& output_file) {
    // Placeholder for BOLSIG+ execution.
    // In a real scenario, we would construct the input file required by BOLSIG+.
    
    std::string cmd = "bolsigplus " + config.cross_section_file + " > " + output_file;
    // Assuming standard output redirection or specific flags.
    
    std::cout << "Running BOLSIG+: " << cmd << std::endl;
    int ret = std::system(cmd.c_str());
    
    if (ret != 0) {
        std::cerr << "Warning: BOLSIG+ execution failed or not found. generating dummy data for testing." << std::endl;
        
        std::ofstream dummy(output_file);
        if (dummy.is_open()) {
            dummy << "Energy Mobility Diffusion Rate" << std::endl;
            dummy << "0.01 100.0 10.0 1e-18" << std::endl;
            dummy << "0.1 80.0 8.0 1e-16" << std::endl;
            dummy << "1.0 50.0 5.0 1e-14" << std::endl;
            dummy << "10.0 10.0 1.0 1e-12" << std::endl;
        }
    }
}

} // namespace HydroPlas
