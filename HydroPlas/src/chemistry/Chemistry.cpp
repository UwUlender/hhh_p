#include "Chemistry.hpp"
#include <cmath>
#include <regex>
#include <iostream>

namespace HydroPlas {

Chemistry::Chemistry(const SimulationConfig& config) {
    // Initialize Species
    for (const auto& sc : config.species) {
        species_.emplace_back(sc);
    }
    
    // Initialize Reactions
    parse_reactions(config.reactions);
}

int Chemistry::get_species_index(const std::string& name) const {
    for(size_t i=0; i<species_.size(); ++i) {
        if(species_[i].name == name) return (int)i;
    }
    return -1;
}

void Chemistry::parse_reactions(const std::vector<ReactionConfig>& r_configs) {
    for (const auto& rc : r_configs) {
        Reaction r;
        r.equation = rc.equation;
        r.type = rc.rate_type;
        r.k_const = rc.a; // Using 'a' field for constant too
        r.A = rc.a;
        r.b = rc.b;
        r.E_a = rc.e_a;
        
        // Load rate table if specified
        if (rc.rate_type == "table" && !rc.table_file.empty()) {
            r.rate_table.load_rate(rc.table_file);
            r.has_rate_table = true;
        }
        
        parse_equation(r, rc.equation);
        reactions_.push_back(r);
    }
}

void Chemistry::parse_equation(Reaction& rxn, const std::string& eq) {
    std::regex token_re("([0-9]*)\\s*([A-Za-z0-9_\\+\\*\\-]+)");
    // Split by "->"
    size_t arrow_pos = eq.find("->");
    if (arrow_pos == std::string::npos) return;
    
    std::string lhs = eq.substr(0, arrow_pos);
    std::string rhs = eq.substr(arrow_pos + 2);
    
    // Parse LHS
    // Split by "+"
    std::stringstream ss_l(lhs);
    std::string segment;
    while(std::getline(ss_l, segment, '+')) {
        // Trim
        std::smatch match;
        if (std::regex_search(segment, match, token_re)) {
            std::string coef_str = match[1];
            std::string name = match[2];
            int coef = coef_str.empty() ? 1 : std::stoi(coef_str);
            int idx = get_species_index(name);
            if(idx >= 0) rxn.reactants[idx] += coef;
        }
    }
    
    // Parse RHS
    std::stringstream ss_r(rhs);
    while(std::getline(ss_r, segment, '+')) {
        std::smatch match;
        if (std::regex_search(segment, match, token_re)) {
            std::string coef_str = match[1];
            std::string name = match[2];
            int coef = coef_str.empty() ? 1 : std::stoi(coef_str);
            int idx = get_species_index(name);
            if(idx >= 0) rxn.products[idx] += coef;
        }
    }
}

double Reaction::get_rate_coeff(double mean_energy, double T_gas) const {
    if (type == "constant") {
        return k_const;
    } else if (type == "arrhenius") {
        return A * std::pow(T_gas, b) * std::exp(-E_a / T_gas);
    } else if (type == "table") {
        if (has_rate_table) {
            return rate_table.get_rate(mean_energy);
        }
        return 0.0;
    }
    return 0.0;
}

void Chemistry::compute_source(const std::vector<double>& densities, double mean_energy, double T_gas, std::vector<double>& sources) const {
    sources.assign(species_.size(), 0.0);
    
    for (const auto& r : reactions_) {
        double k = r.get_rate_coeff(mean_energy, T_gas);
        
        // Calculate rate of reaction R = k * prod(n_reactants)
        double R_rate = k;
        for (const auto& pair : r.reactants) {
            int idx = pair.first;
            int stoich = pair.second;
            R_rate *= std::pow(densities[idx], stoich);
        }
        
        // Update sources for all participants
        // S_k += (nu_R - nu_L) * R_rate
        for (const auto& pair : r.reactants) {
            int idx = pair.first;
            sources[idx] -= pair.second * R_rate;
        }
        for (const auto& pair : r.products) {
            int idx = pair.first;
            sources[idx] += pair.second * R_rate;
        }
    }
}

} // namespace HydroPlas
