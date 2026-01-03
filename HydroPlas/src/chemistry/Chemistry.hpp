#pragma once
#include <vector>
#include <string>
#include <map>
#include "Species.hpp"
#include "../config/ConfigParser.hpp"
#include "../numerics/EquationEvaluator.hpp"

namespace HydroPlas {

struct Reaction {
    std::string equation;
    // Map species index to stoichiometric coefficient
    std::map<int, int> reactants;
    std::map<int, int> products;
    
    // Rate parameters
    std::string type; // "constant", "arrhenius", "table", "equation"
    double k_const;
    double A, b, E_a;
    
    // Table lookup for electron-impact reactions
    LookupTable rate_table;
    bool has_rate_table = false;
    
    // Equation based rate
    EquationEvaluator evaluator;
    bool has_equation = false;
    
    // Helper to calculate rate coefficient
    // Updated to accept variable map for equation types
    double get_rate_coeff(double mean_energy, double T_gas, const std::map<std::string, double>& vars = {}) const;

    std::string name; // Name for output
};

class Chemistry {
public:
    Chemistry(const SimulationConfig& config);
    
    // Core method: Compute source terms for all species
    // input: densities of all species (same order as species_list), T_e, T_gas
    // output: S vector (size = num_species)
    void compute_source(const std::vector<double>& densities, double mean_energy, double T_gas, std::vector<double>& sources) const;

    const std::vector<Species>& get_species() const { return species_; }
    int get_num_species() const { return species_.size(); }
    int get_species_index(const std::string& name) const;
    
    // Access to reactions for rate saving
    const std::vector<Reaction>& get_reactions() const { return reactions_; }

private:
    std::vector<Species> species_;
    std::vector<Reaction> reactions_;
    
    void parse_reactions(const std::vector<ReactionConfig>& r_configs);
    void parse_equation(Reaction& rxn, const std::string& eq);
};

} // namespace HydroPlas
