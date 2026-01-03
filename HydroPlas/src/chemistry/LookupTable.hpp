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
    
    // Load single-column rate data (Energy, Rate)
    void load_rate(const std::string& filename);

    // Interpolate values based on Mean Energy (epsilon)
    // Returns: mobility, diffusion, temperature(if needed)
    double get_mobility(double energy) const;
    double get_diffusion(double energy) const;
    
    // For reaction rates
    double get_rate(double energy) const;

    // Scale transport data (mobility and diffusion) by factor
    void scale_transport(double factor);

    // Scale rate data by factor
    void scale_rate(double factor);

private:
    std::vector<double> energy_grid_;
    std::vector<double> mobility_data_;
    std::vector<double> diffusion_data_;
    std::vector<double> rate_data_;
    
    double interpolate(const std::vector<double>& x, const std::vector<double>& y, double val) const;
};

} // namespace HydroPlas
