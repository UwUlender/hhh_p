#pragma once
#include <string>
#include "../config/ConfigParser.hpp"

namespace HydroPlas {

/**
 * @brief Interface to BOLSIG+ for generating electron transport coefficients
 * 
 * This class provides methods to:
 * 1. Create BOLSIG+ input files from cross-section data
 * 2. Execute BOLSIG+ (if available on system)
 * 3. Parse BOLSIG+ output files
 * 4. Generate fallback transport data if BOLSIG+ unavailable
 * 
 * BOLSIG+ is a popular two-term Boltzmann equation solver for electron kinetics.
 * Website: www.bolsig.laplace.univ-tlse.fr
 */
class BolsigInterface {
public:
    /**
     * @brief Main entry point: Generate transport table from cross-sections
     * 
     * Workflow:
     * 1. Check if cross-section file exists
     * 2. Create BOLSIG+ input file
     * 3. Execute BOLSIG+ 
     * 4. Parse output and convert to HydroPlas format
     * 5. Fall back to analytical approximations if BOLSIG+ fails
     * 
     * @param config Chemistry configuration (includes cross-section file path, gas temp)
     * @param output_file Output file path for HydroPlas-format lookup table
     */
    static void run_bolsig(const ChemistryConfig& config, const std::string& output_file);
    
private:
    /**
     * @brief Create BOLSIG+ input file from configuration
     */
    static void create_bolsig_input(const ChemistryConfig& config, const std::string& input_file);
    
    /**
     * @brief Parse BOLSIG+ output and convert to lookup table format
     * @return true if parsing successful, false otherwise
     */
    static bool parse_bolsig_output(const std::string& bolsig_output, const std::string& table_output);
    
    /**
     * @brief Generate fallback transport data using analytical approximations
     * 
     * Used when BOLSIG+ is not available or fails. Generates realistic
     * Argon transport data using empirical formulas and Arrhenius rates.
     */
    static void generate_fallback_data(const std::string& output_file);
};

} // namespace HydroPlas
