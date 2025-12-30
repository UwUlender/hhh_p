#pragma once
#include <vector>
#include <string>
#include <map>
#include "../config/ConfigParser.hpp"
#include "LookupTable.hpp"

namespace HydroPlas {

/**
 * @brief ReactionHandler: Computes reaction rates and source terms for plasma chemistry
 * 
 * Implements the following reaction types as described in the theoretical framework:
 * 1. Direct ionization: e + A -> 2e + A+
 * 2. Excitation: e + A -> e + A*
 * 3. Stepwise ionization: e + A* -> 2e + A+ (two-step ladder)
 * 4. Penning ionization: A* + B -> A + B+ + e
 * 5. Superelastic collision: e + A* -> e(fast) + A (electron heating)
 * 6. Radiative decay: A* -> A + hv
 * 7. Collisional quenching: A* + M -> A + M
 * 8. Three-body recombination: e + A+ + M -> A + M
 * 9. Metastable pooling: A* + A* -> A+ + A + e
 */
class ReactionHandler {
public:
    ReactionHandler(const ChemistryConfig& config, LookupTable* lookup);
    
    /**
     * @brief Compute all reaction source terms for a given state
     * @param ne Electron density [m^-3]
     * @param ni Ion density [m^-3]
     * @param n_excited Vector of excited species densities [m^-3]
     * @param mean_energy Mean electron energy [eV]
     * @param N_gas Background gas density [m^-3]
     * @param S_ne Output: electron source term [m^-3 s^-1]
     * @param S_ni Output: ion source term [m^-3 s^-1]
     * @param S_neps Output: energy source term [eV m^-3 s^-1]
     * @param S_excited Output: excited species source terms [m^-3 s^-1]
     */
    void compute_sources(double ne, double ni, const std::vector<double>& n_excited,
                        double mean_energy, double N_gas,
                        double& S_ne, double& S_ni, double& S_neps,
                        std::vector<double>& S_excited);
    
    /**
     * @brief Get the energy loss/gain associated with a reaction
     */
    double get_energy_cost(const std::string& reaction_type, int species_idx = -1);
    
private:
    const ChemistryConfig& config_;
    LookupTable* lookup_;
    
    // Cached rate coefficients
    std::map<std::string, double> rate_cache_;
    
    // Helper functions for specific reaction types
    void compute_direct_ionization(double ne, double N_gas, double mean_energy,
                                   double& S_ne, double& S_ni, double& S_neps);
    
    void compute_excitation(double ne, double N_gas, double mean_energy,
                           std::vector<double>& S_excited, double& S_neps);
    
    void compute_stepwise_ionization(double ne, const std::vector<double>& n_excited,
                                    double mean_energy,
                                    double& S_ne, double& S_ni, double& S_neps,
                                    std::vector<double>& S_excited);
    
    void compute_penning_ionization(const std::vector<double>& n_excited,
                                   double& S_ne, double& S_ni,
                                   std::vector<double>& S_excited);
    
    void compute_superelastic(double ne, const std::vector<double>& n_excited,
                             double mean_energy, double& S_neps,
                             std::vector<double>& S_excited);
    
    void compute_radiative_decay(std::vector<double>& S_excited);
    
    void compute_quenching(const std::vector<double>& n_excited, double N_gas,
                          std::vector<double>& S_excited);
    
    void compute_pooling(const std::vector<double>& n_excited,
                        double& S_ne, double& S_ni,
                        std::vector<double>& S_excited);
    
    void compute_recombination(double ne, double ni, double N_gas,
                              double& S_ne, double& S_ni, double& S_neps);
};

} // namespace HydroPlas
