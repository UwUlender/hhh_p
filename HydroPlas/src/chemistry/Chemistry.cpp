#include "Chemistry.hpp"
#include <cmath>
#include <regex>
#include <iostream>
#include <stdexcept>

namespace HydroPlas {

Chemistry::Chemistry(const SimulationConfig& config) {
    // Calculate neutral gas density N = P / (k_B * T)
    // k_B = 1.380649e-23
    double k_B = 1.380649e-23;
    double N_gas = config.plasma.gas_pressure / (k_B * config.plasma.background_temp);
    
    // Initialize Species
    for (const auto& sc : config.species) {
        species_.emplace_back(sc);
        // Scale mobility/diffusion if loaded from table
        // Table values are multiplied by N, so we divide by N
        if (species_.back().has_lookup) {
            species_.back().lookup_table.scale_transport(1.0 / N_gas);
        }
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
            if (!r.rate_table.load_rate(rc.table_file)) {
                 throw std::runtime_error("Failed to load rate table: " + rc.table_file);
            }
            r.has_rate_table = true;
            
            // Table values are multiplied by Avogadro constant N_A
            // So we divide by N_A to get the actual rate coefficient
            double N_A = 6.02214076e23;
            r.rate_table.scale_rate(1.0 / N_A);
        }
        
        // Parse equation if specified
        if (rc.rate_type == "equation" && !rc.rate_equation.empty()) {
            r.has_equation = true;
            r.evaluator.parse(rc.rate_equation);
            
            // Parse constants
            std::stringstream ss_names(rc.equation_constants);
            std::stringstream ss_vals(rc.equation_values);
            std::string name;
            double val;
            while(ss_names >> name && ss_vals >> val) {
                r.evaluator.set_constant(name, val);
            }
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

double Reaction::get_rate_coeff(double mean_energy, double T_gas, const std::map<std::string, double>& vars) const {
    if (type == "constant") {
        return k_const;
    } else if (type == "arrhenius") {
        return A * std::pow(T_gas, b) * std::exp(-E_a / T_gas);
    } else if (type == "table") {
        if (has_rate_table) {
            return rate_table.get_rate(mean_energy);
        }
        return 0.0;
    } else if (type == "equation") {
        if (has_equation) {
            // Add Tgas if not present? The equation parser uses variables provided.
            std::map<std::string, double> v = vars;
            v["Tgas"] = T_gas;
            v["mean_energy"] = mean_energy;
            return evaluator.evaluate(v);
        }
    }
    return 0.0;
}

void Chemistry::compute_source(const std::vector<double>& densities, double mean_energy, double T_gas, std::vector<double>& sources) const {
    sources.assign(species_.size(), 0.0);
    
    // Construct variable map for equations once per cell
    // This is overhead but necessary for generic equation support
    std::map<std::string, double> vars;
    bool needs_vars = false;
    for(const auto& r : reactions_) {
        if(r.type == "equation") { needs_vars = true; break; }
    }
    
    if(needs_vars) {
        for(size_t i=0; i<species_.size(); ++i) {
            vars[species_[i].name + "_density"] = densities[i];
             // Also support just name?
            vars[species_[i].name] = densities[i];
        }
    }
    
    for (const auto& r : reactions_) {
        double k = r.get_rate_coeff(mean_energy, T_gas, vars);
        
        // Calculate rate of reaction R
        double R_rate = 0.0;
        
        if (r.type == "equation") {
             // For equation type, we assume the equation returns the FULL Rate (R) 
             // OR the rate coefficient (k).
             // Standard practice: if equation uses density, it likely returns R directly 
             // (e.g. radiation).
             // However, if the user defines "Ar1s3 -> Ar1s5", usually we multiply by [Ar1s3].
             // BUT if the equation *already* contains [Ar1s5] (which is product/trapping?),
             // it might be the net rate.
             // Given the complexity of the equation, let's assume it returns the NET rate.
             // Wait, if it returns k, then R = k * [Reactants].
             // If the equation returns a value ~ 1e7, it looks like a rate coefficient (1/s for first order).
             // If it returns ~ 1e20, it's a rate.
             // Looking at the equation: 1.89e7 is in there. 
             // Let's assume the equation returns 'k' (coefficient) effectively, 
             // OR it returns the source term.
             // The config says "rate_type: equation".
             // If I treat it as 'k', I multiply by reactants.
             // If the equation calculates the effective Einstein A coefficient (radiative decay rate),
             // then R = A_eff * [Ar1s3].
             // The equation provided depends on [Ar1s5] (trapping factor?).
             // I will assume it returns the COEFFICIENT (1/s or m3/s).
             // Because for "radia1", the reaction is Ar1s3 -> Ar1s5.
             // So R = k * [Ar1s3].
             
             R_rate = k;
             // Multiply by reactants stoichiometry
             for (const auto& pair : r.reactants) {
                 int idx = pair.first;
                 int stoich = pair.second;
                 R_rate *= std::pow(densities[idx], stoich);
             }
        } else {
             R_rate = k;
             for (const auto& pair : r.reactants) {
                 int idx = pair.first;
                 int stoich = pair.second;
                 R_rate *= std::pow(densities[idx], stoich);
             }
        }
        
        // Update sources...
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
