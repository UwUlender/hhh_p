#include "LookupTable.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace HydroPlas {

void LookupTable::load_from_file(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open lookup table: " + filepath);
    }

    // Assume first line is header: Energy Mobility Diff ...
    std::string line;
    // Skip comments
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '#') break;
    }
    
    std::stringstream ss(line);
    std::string col;
    std::vector<std::string> headers;
    while (ss >> col) {
        headers.push_back(col);
    }
    
    // Resize vectors
    for (const auto& h : headers) {
        if (h != "Energy") data_columns_[h] = {};
    }

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream dss(line);
        double val;
        int idx = 0;
        double energy_val = 0;
        
        // Read line
        std::vector<double> row_values;
        while (dss >> val) {
            row_values.push_back(val);
        }

        if (row_values.size() != headers.size()) continue;

        energy_grid_.push_back(row_values[0]);
        for (size_t i = 1; i < headers.size(); ++i) {
            data_columns_[headers[i]].push_back(row_values[i]);
        }
    }
}

double LookupTable::interpolate(double energy, const std::string& column_name) const {
    if (energy_grid_.empty()) return 0.0;
    if (data_columns_.find(column_name) == data_columns_.end()) {
        throw std::runtime_error("Column not found: " + column_name);
    }

    const auto& y_vals = data_columns_.at(column_name);

    // Find interval
    auto it = std::lower_bound(energy_grid_.begin(), energy_grid_.end(), energy);
    if (it == energy_grid_.begin()) return y_vals.front();
    if (it == energy_grid_.end()) return y_vals.back();

    int idx = std::distance(energy_grid_.begin(), it) - 1;
    double x0 = energy_grid_[idx];
    double x1 = energy_grid_[idx+1];
    double y0 = y_vals[idx];
    double y1 = y_vals[idx+1];

    // Log-Log interpolation
    if (y0 > 1e-30 && y1 > 1e-30 && x0 > 1e-30 && x1 > 1e-30 && energy > 1e-30) {
        double log_y = std::log(y0) + (std::log(y1) - std::log(y0)) * (std::log(energy) - std::log(x0)) / (std::log(x1) - std::log(x0));
        return std::exp(log_y);
    } else {
        return y0 + (y1 - y0) * (energy - x0) / (x1 - x0);
    }
}

} // namespace HydroPlas
