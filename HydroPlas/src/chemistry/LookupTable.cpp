#include "LookupTable.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace HydroPlas {

void LookupTable::load(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin.is_open()) return; // Warning log?

    // Generic simple parser for 2-column or multicolumn
    // Assuming: Energy | Mobility | Diffusion
    // Skip headers
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || !isdigit(line[0])) continue; // rudimentary skip
        std::stringstream ss(line);
        double e, mu, d;
        ss >> e >> mu >> d;
        energy_grid_.push_back(e);
        mobility_data_.push_back(mu);
        diffusion_data_.push_back(d);
    }
}

void LookupTable::load_rate(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Warning: Cannot open rate table file: " << filename << std::endl;
        return;
    }

    // Format: Energy | Rate
    // Skip comment lines starting with # or non-digit
    std::string line;
    energy_grid_.clear();
    rate_data_.clear();
    
    while (std::getline(fin, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#' || !isdigit(line[0])) continue;
        
        std::stringstream ss(line);
        double e, rate;
        if (ss >> e >> rate) {
            energy_grid_.push_back(e);
            rate_data_.push_back(rate);
        }
    }
    
    if (energy_grid_.empty()) {
        std::cerr << "Warning: No data loaded from rate table: " << filename << std::endl;
    }
}

double LookupTable::interpolate(const std::vector<double>& x, const std::vector<double>& y, double val) const {
    if (x.empty()) return 0.0;
    if (val <= x.front()) return y.front();
    if (val >= x.back()) return y.back();

    auto it = std::lower_bound(x.begin(), x.end(), val);
    size_t i = std::distance(x.begin(), it);
    if (i == 0) i = 1;

    double x1 = x[i-1];
    double x2 = x[i];
    double y1 = y[i-1];
    double y2 = y[i];

    // Log-Log interpolation
    // Avoid log(0)
    if (x1 <= 1e-20 || val <= 1e-20) {
        // Fallback to linear
        double m = (y2 - y1) / (x2 - x1);
        return y1 + m * (val - x1);
    }

    double log_x1 = std::log(x1);
    double log_x2 = std::log(x2);
    double log_val = std::log(val);
    
    double log_y1 = (y1 > 0) ? std::log(y1) : -100.0; // Guard
    double log_y2 = (y2 > 0) ? std::log(y2) : -100.0;

    double m = (log_y2 - log_y1) / (log_x2 - log_x1);
    double log_res = log_y1 + m * (log_val - log_x1);
    
    return std::exp(log_res);
}

double LookupTable::get_mobility(double energy) const {
    return interpolate(energy_grid_, mobility_data_, energy);
}

double LookupTable::get_diffusion(double energy) const {
    return interpolate(energy_grid_, diffusion_data_, energy);
}

double LookupTable::get_rate(double energy) const {
    if (rate_data_.empty()) return 0.0;
    return interpolate(energy_grid_, rate_data_, energy);
}

} // namespace HydroPlas
