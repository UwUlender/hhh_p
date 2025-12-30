#include "ReactionHandler.hpp"
#include <cmath>
#include <iostream>

namespace HydroPlas {

ReactionHandler::ReactionHandler(const ChemistryConfig& config, LookupTable* lookup)
    : config_(config), lookup_(lookup) {
}

void ReactionHandler::compute_sources(double ne, double ni, 
                                     const std::vector<double>& n_excited,
                                     double mean_energy, double N_gas,
                                     double& S_ne, double& S_ni, double& S_neps,
                                     std::vector<double>& S_excited) {
    // Initialize all sources to zero
    S_ne = 0.0;
    S_ni = 0.0;
    S_neps = 0.0;
    S_excited.assign(n_excited.size(), 0.0);
    
    // Enforce minimum densities to prevent numerical issues
    ne = std::max(ne, 1e10);
    ni = std::max(ni, 1e10);
    N_gas = std::max(N_gas, 1e20);
    
    // 1. Direct ionization: e + A -> 2e + A+
    compute_direct_ionization(ne, N_gas, mean_energy, S_ne, S_ni, S_neps);
    
    // 2. Excitation: e + A -> e + A*
    compute_excitation(ne, N_gas, mean_energy, S_excited, S_neps);
    
    // 3. Stepwise ionization: e + A* -> 2e + A+
    compute_stepwise_ionization(ne, n_excited, mean_energy, S_ne, S_ni, S_neps, S_excited);
    
    // 4. Penning ionization: A* + A* -> A+ + A + e
    compute_penning_ionization(n_excited, S_ne, S_ni, S_excited);
    
    // 5. Superelastic collisions: e + A* -> e(fast) + A
    compute_superelastic(ne, n_excited, mean_energy, S_neps, S_excited);
    
    // 6. Radiative decay: A* -> A + hv
    compute_radiative_decay(S_excited);
    
    // 7. Collisional quenching: A* + M -> A + M
    compute_quenching(n_excited, N_gas, S_excited);
    
    // 8. Metastable pooling: A* + A* -> A+ + A + e
    compute_pooling(n_excited, S_ne, S_ni, S_excited);
    
    // 9. Three-body recombination: e + A+ + M -> A + M
    compute_recombination(ne, ni, N_gas, S_ne, S_ni, S_neps);
}

void ReactionHandler::compute_direct_ionization(double ne, double N_gas, double mean_energy,
                                               double& S_ne, double& S_ni, double& S_neps) {
    // Direct ionization: e + A -> 2e + A+
    // Rate = k_iz(Te) * ne * N_gas
    // Energy cost: E_iz (typically 15.76 eV for Argon)
    
    double k_iz = 0.0;
    try {
        k_iz = lookup_->interpolate(mean_energy, "Rate_Ionization");
    } catch (...) {
        // If not in table, use Arrhenius form: k = A * exp(-E_iz / mean_energy)
        double E_iz = 15.76; // eV (Argon)
        if (mean_energy > 1.0) {
            k_iz = 1e-14 * std::exp(-E_iz / mean_energy);
        }
    }
    
    double R_iz = k_iz * ne * N_gas;
    
    S_ne += R_iz;   // Gain one electron per ionization
    S_ni += R_iz;   // Gain one ion per ionization
    S_neps -= R_iz * 15.76; // Energy loss (E_iz = 15.76 eV for Ar)
}

void ReactionHandler::compute_excitation(double ne, double N_gas, double mean_energy,
                                        std::vector<double>& S_excited, double& S_neps) {
    // Excitation: e + A -> e + A*
    // Rate = k_exc(Te) * ne * N_gas
    // Energy cost: E_exc (e.g., 11.5 eV for Ar metastable)
    
    for (size_t k = 0; k < config_.excited_species.size(); ++k) {
        double k_exc = 0.0;
        double E_exc = config_.excited_species[k].energy_level; // eV
        
        try {
            // Try to get from lookup table
            std::string col_name = "Rate_Excitation_" + config_.excited_species[k].name;
            k_exc = lookup_->interpolate(mean_energy, col_name);
        } catch (...) {
            // Use generic excitation rate
            try {
                k_exc = lookup_->interpolate(mean_energy, "Rate_Excitation");
            } catch (...) {
                // Fallback Arrhenius form
                if (mean_energy > E_exc / 2.0) {
                    k_exc = 5e-15 * std::exp(-E_exc / mean_energy);
                }
            }
        }
        
        double R_exc = k_exc * ne * N_gas;
        S_excited[k] += R_exc;  // Gain excited species
        S_neps -= R_exc * E_exc; // Energy loss
    }
}

void ReactionHandler::compute_stepwise_ionization(double ne, 
                                                  const std::vector<double>& n_excited,
                                                  double mean_energy,
                                                  double& S_ne, double& S_ni, double& S_neps,
                                                  std::vector<double>& S_excited) {
    // Stepwise ionization: e + A* -> 2e + A+
    // This is the "ladder" effect described in the document
    // The threshold energy is much lower: E_step ~ 4.2 eV for Ar* -> Ar+
    
    for (size_t k = 0; k < config_.excited_species.size(); ++k) {
        double E_state = config_.excited_species[k].energy_level;
        double E_iz = 15.76; // Total ionization potential
        double E_step = E_iz - E_state; // Remaining energy needed
        
        double k_step = 0.0;
        try {
            std::string col_name = "Rate_Stepwise_" + config_.excited_species[k].name;
            k_step = lookup_->interpolate(mean_energy, col_name);
        } catch (...) {
            // Fallback: use modified ionization rate with lower threshold
            if (mean_energy > E_step / 2.0) {
                k_step = 1e-13 * std::exp(-E_step / mean_energy);
            }
        }
        
        double n_k = std::max(n_excited[k], 1e10);
        double R_step = k_step * ne * n_k;
        
        S_ne += R_step;        // Gain electron
        S_ni += R_step;        // Gain ion
        S_excited[k] -= R_step; // Lose excited species
        S_neps -= R_step * E_step; // Energy loss (only the remaining barrier)
    }
}

void ReactionHandler::compute_penning_ionization(const std::vector<double>& n_excited,
                                                 double& S_ne, double& S_ni,
                                                 std::vector<double>& S_excited) {
    // Penning ionization: A* + B -> A + B+ + e
    // For Argon: Ar* + Ar -> Ar + Ar+ + e (Hornbeck-Molnar process)
    // Also relevant in mixtures: He* + N2 -> He + N2+ + e
    
    // Self-Penning for Argon metastables
    for (size_t k = 0; k < config_.excited_species.size(); ++k) {
        const auto& species = config_.excited_species[k];
        
        // Check if species name contains "m" (metastable)
        if (species.name.find("_m") != std::string::npos || 
            species.name.find("metastable") != std::string::npos) {
            
            // Rate coefficient for Penning (typically ~1e-16 to 1e-15 m^3/s)
            double k_penning = 1e-15; // m^3/s for Ar* + Ar
            
            double n_k = std::max(n_excited[k], 1e10);
            double R_penning = k_penning * n_k * n_k; // Second order in metastable
            
            S_ne += R_penning;         // Gain electron
            S_ni += R_penning;         // Gain ion
            S_excited[k] -= 2.0 * R_penning; // Lose two metastables
        }
    }
}

void ReactionHandler::compute_superelastic(double ne, const std::vector<double>& n_excited,
                                          double mean_energy, double& S_neps,
                                          std::vector<double>& S_excited) {
    // Superelastic collision: e + A* -> e(fast) + A
    // This is the reverse of excitation, heating the electron gas
    // Critical in afterglows where E ~ 0
    
    for (size_t k = 0; k < config_.excited_species.size(); ++k) {
        double E_exc = config_.excited_species[k].energy_level;
        
        // Superelastic rate ~ de-excitation rate * Boltzmann factor
        // k_superelastic ~ k_exc * exp(-E_exc / Te) / detailed balance
        double k_superelastic = 0.0;
        if (mean_energy > 0.1) {
            k_superelastic = 1e-14 * std::exp(-E_exc / (2.0 * mean_energy));
        }
        
        double n_k = std::max(n_excited[k], 1e10);
        double R_superelastic = k_superelastic * ne * n_k;
        
        S_excited[k] -= R_superelastic;    // Lose excited species
        S_neps += R_superelastic * E_exc;  // Electron energy GAIN (heating)
    }
}

void ReactionHandler::compute_radiative_decay(std::vector<double>& S_excited) {
    // Radiative decay: A* -> A + hv
    // Rate = A_rad (Einstein A coefficient) * n*
    // For resonant states: A_rad ~ 1e8 s^-1 (ns lifetime)
    // For metastables: A_rad ~ 0 (forbidden transition)
    
    for (size_t k = 0; k < config_.excited_species.size(); ++k) {
        double A_rad = 0.0;
        
        // Check if resonant state (fast radiative decay)
        if (config_.excited_species[k].name.find("_r") != std::string::npos ||
            config_.excited_species[k].name.find("resonant") != std::string::npos) {
            A_rad = 1e8; // s^-1 (1/10ns)
        } else if (config_.excited_species[k].name.find("_m") != std::string::npos ||
                   config_.excited_species[k].name.find("metastable") != std::string::npos) {
            A_rad = 0.0; // Forbidden transition (long lifetime)
        } else {
            A_rad = 1e6; // Default intermediate value
        }
        
        double n_k = std::max(n_excited[k], 0.0); // No lower limit for radiative
        double R_rad = A_rad * n_k;
        
        S_excited[k] -= R_rad; // Lose excited species
        // Note: Photon energy is lost from the system (no energy source term)
    }
}

void ReactionHandler::compute_quenching(const std::vector<double>& n_excited, double N_gas,
                                       std::vector<double>& S_excited) {
    // Collisional quenching: A* + M -> A + M
    // Rate = k_Q * n* * N_gas
    // Transfers energy to gas (not tracked in energy equation)
    
    for (size_t k = 0; k < config_.excited_species.size(); ++k) {
        // Quenching rate coefficient (typically 1e-17 to 1e-15 m^3/s)
        double k_Q = 1e-16; // m^3/s
        
        double n_k = std::max(n_excited[k], 1e10);
        double R_quench = k_Q * n_k * N_gas;
        
        S_excited[k] -= R_quench; // Lose excited species
    }
}

void ReactionHandler::compute_pooling(const std::vector<double>& n_excited,
                                     double& S_ne, double& S_ni,
                                     std::vector<double>& S_excited) {
    // Metastable pooling: A* + A* -> A+ + A + e
    // Similar to Penning but both reactants are the same metastable
    // Important memory effect for breakdown
    
    for (size_t k = 0; k < config_.excited_species.size(); ++k) {
        // Check for metastable species
        if (config_.excited_species[k].name.find("_m") != std::string::npos) {
            double k_pool = 5e-16; // m^3/s
            
            double n_k = std::max(n_excited[k], 1e10);
            double R_pool = k_pool * n_k * n_k;
            
            S_ne += R_pool;
            S_ni += R_pool;
            S_excited[k] -= 2.0 * R_pool; // Two metastables consumed
        }
    }
}

void ReactionHandler::compute_recombination(double ne, double ni, double N_gas,
                                           double& S_ne, double& S_ni, double& S_neps) {
    // Three-body recombination: e + A+ + M -> A + M
    // Rate = k_rec * ne * ni * N_gas
    // Energy gain (recombination releases ionization energy)
    
    double k_rec = 1e-39; // m^6/s (three-body rate)
    
    double R_rec = k_rec * ne * ni * N_gas;
    
    S_ne -= R_rec;  // Lose electron
    S_ni -= R_rec;  // Lose ion
    S_neps += R_rec * 15.76; // Energy gain (exothermic)
}

double ReactionHandler::get_energy_cost(const std::string& reaction_type, int species_idx) {
    if (reaction_type == "ionization") {
        return 15.76; // eV for Argon
    } else if (reaction_type == "excitation" && species_idx >= 0) {
        return config_.excited_species[species_idx].energy_level;
    } else if (reaction_type == "stepwise" && species_idx >= 0) {
        return 15.76 - config_.excited_species[species_idx].energy_level;
    }
    return 0.0;
}

} // namespace HydroPlas
