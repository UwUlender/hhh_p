#pragma once
#include <vector>
#include <string>

namespace HydroPlas {

class LookupTable {
public:
    LookupTable() = default;
    
    // Load from BOLSIG+ text file
    // Expects specific column mapping
    void load(const std::string& filename);

    // Interpolate values based on Mean Energy (epsilon)
    // Returns: mobility, diffusion, temperature(if needed)
    double get_mobility(double energy) const;
    double get_diffusion(double energy) const;
    
    // For reaction rates
    double get_rate(double energy, int reaction_index) const;

private:
    std::vector<double> energy_grid_;
    std::vector<double> mobility_data_;
    std::vector<double> diffusion_data_;
    // Map reaction index to data vector
    // For simplicity, we might need a more complex mapping if multiple reactions use the table
    // or just store raw columns.
    
    double interpolate(const std::vector<double>& x, const std::vector<double>& y, double val) const;
};

} // namespace HydroPlas
