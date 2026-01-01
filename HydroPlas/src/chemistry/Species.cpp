#include "Species.hpp"

namespace HydroPlas {

Species::Species(const SpeciesConfig& config) {
    name = config.name;
    if (config.type == "electron") type = SpeciesType::Electron;
    else if (config.type == "ion") type = SpeciesType::Ion;
    else type = SpeciesType::Neutral;

    charge = config.charge;
    mass = config.mass;
    diffusion_coeff_const = config.diffusion_coeff;

    if (!config.mobility_file.empty()) {
        lookup_table.load(config.mobility_file);
        has_lookup = true;
    }
}

void Species::get_transport(double mean_energy, double& mu, double& D) const {
    if (type == SpeciesType::Neutral) {
        mu = 0.0;
        D = diffusion_coeff_const;
    } else {
        if (has_lookup) {
            mu = lookup_table.get_mobility(mean_energy);
            D = lookup_table.get_diffusion(mean_energy);
        } else {
            // Fallback or constant approximations
            // E.g. Einstein relation if D is not in table?
            // For now, assume table handles it or 0.
            mu = 0.0; 
            D = 0.0;
        }
    }
}

} // namespace HydroPlas
