#pragma once
#include <string>
#include <vector>
#include "LookupTable.hpp"
#include "../config/ConfigParser.hpp"

namespace HydroPlas {

enum class SpeciesType { Electron, Ion, Neutral };

class Species {
public:
    Species(const SpeciesConfig& config);
    
    std::string name;
    SpeciesType type;
    double charge;
    double mass;
    double diffusion_coeff_const; // For neutrals
    
    LookupTable lookup_table;
    bool has_lookup = false;

    // Get transport coefficients
    // For charged: interpolates from table using mean energy (eV)
    // For neutral: returns constant D, mobility=0
    // Returns: mobility (m2/V/s), diffusion (m2/s)
    void get_transport(double mean_energy, double& mu, double& D) const;
};

class Chemistry; // Forward

} // namespace HydroPlas
