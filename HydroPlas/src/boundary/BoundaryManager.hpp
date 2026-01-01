#pragma once
#include "../config/ConfigParser.hpp"
#include <petsc.h>
#include <string>
#include <map>
#include <vector>

namespace HydroPlas {

// Must match usage in BoundaryManager.cpp
struct BoundaryConfig {
    bool use_multi_electrode = true;
    std::vector<ElectrodeConfig> electrodes;
    
    // Legacy support fields if needed (or we migrate logic)
    std::string voltage_type = "DC";
    double voltage_amplitude = 0.0;
    double frequency = 0.0;
    double bias = 0.0;
    double gamma_see = 0.0;
    double dielectric_permittivity = 1.0;
    bool is_dielectric = false;
};

class BoundaryManager {
public:
    explicit BoundaryManager(const BoundaryConfig& config);
    
    // Legacy interface (for backward compatibility)
    double get_voltage(double t) const;
    double get_secondary_emission_flux(double ion_flux) const;
    
    // New multi-electrode interface
    double get_electrode_voltage(const std::string& electrode_name, double t) const;
    double get_electrode_voltage_by_index(int electrode_index, double t) const;
    double get_electrode_gamma_see(const std::string& electrode_name) const;
    double get_electrode_gamma_see_by_index(int electrode_index) const;
    bool is_electrode_dielectric(const std::string& electrode_name) const;
    bool is_electrode_dielectric_by_index(int electrode_index) const;
    double get_electrode_dielectric_permittivity(const std::string& electrode_name) const;
    double get_electrode_dielectric_thickness(const std::string& electrode_name) const;
    int get_num_electrodes() const;
    
private:
    BoundaryConfig config_;
    std::map<std::string, int> electrode_name_to_index_;
    
    double compute_voltage(const ElectrodeConfig& electrode, double t) const;
    void build_electrode_map();
};

} // namespace HydroPlas
