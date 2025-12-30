#include "BoundaryManager.hpp"
#include <cmath>
#include <stdexcept>

namespace HydroPlas {

BoundaryManager::BoundaryManager(const BoundaryConfig& config) : config_(config) {
    build_electrode_map();
}

void BoundaryManager::build_electrode_map() {
    for (size_t i = 0; i < config_.electrodes.size(); ++i) {
        electrode_name_to_index_[config_.electrodes[i].name] = i;
    }
}

double BoundaryManager::compute_voltage(const ElectrodeConfig& electrode, double t) const {
    if (electrode.voltage_type == "DC") {
        return electrode.voltage_amplitude + electrode.bias;
    } 
    else if (electrode.voltage_type == "RF" || electrode.voltage_type == "AC") {
        return electrode.voltage_amplitude * std::sin(2.0 * M_PI * electrode.frequency * t + electrode.phase) + electrode.bias;
    }
    else if (electrode.voltage_type == "PULSE") {
        // Square wave with duty cycle
        double period = 1.0 / electrode.frequency;
        double t_mod = std::fmod(t, period);
        double t_on = period * electrode.duty_cycle;
        if (t_mod < t_on) {
            return electrode.voltage_amplitude + electrode.bias;
        } else {
            return electrode.bias;
        }
    }
    return 0.0;
}

// Legacy interface (backward compatibility)
double BoundaryManager::get_voltage(double t) const {
    if (config_.use_multi_electrode && !config_.electrodes.empty()) {
        // Use first electrode for legacy calls
        return compute_voltage(config_.electrodes[0], t);
    }
    
    // Old single-electrode behavior
    if (config_.voltage_type == "DC") {
        return config_.voltage_amplitude;
    } else if (config_.voltage_type == "RF") {
        return config_.voltage_amplitude * std::sin(2.0 * M_PI * config_.frequency * t) + config_.bias;
    }
    return 0.0;
}

double BoundaryManager::get_secondary_emission_flux(double ion_flux) const {
    return config_.gamma_see * ion_flux;
}

// Multi-electrode interface
double BoundaryManager::get_electrode_voltage(const std::string& electrode_name, double t) const {
    auto it = electrode_name_to_index_.find(electrode_name);
    if (it == electrode_name_to_index_.end()) {
        throw std::runtime_error("Electrode not found: " + electrode_name);
    }
    return get_electrode_voltage_by_index(it->second, t);
}

double BoundaryManager::get_electrode_voltage_by_index(int electrode_index, double t) const {
    if (!config_.use_multi_electrode || config_.electrodes.empty()) {
        // Fall back to legacy behavior
        return get_voltage(t);
    }
    
    if (electrode_index < 0 || electrode_index >= static_cast<int>(config_.electrodes.size())) {
        throw std::runtime_error("Electrode index out of range: " + std::to_string(electrode_index));
    }
    
    return compute_voltage(config_.electrodes[electrode_index], t);
}

double BoundaryManager::get_electrode_gamma_see(const std::string& electrode_name) const {
    auto it = electrode_name_to_index_.find(electrode_name);
    if (it == electrode_name_to_index_.end()) {
        throw std::runtime_error("Electrode not found: " + electrode_name);
    }
    return get_electrode_gamma_see_by_index(it->second);
}

double BoundaryManager::get_electrode_gamma_see_by_index(int electrode_index) const {
    if (!config_.use_multi_electrode || config_.electrodes.empty()) {
        return config_.gamma_see;
    }
    
    if (electrode_index < 0 || electrode_index >= static_cast<int>(config_.electrodes.size())) {
        return config_.gamma_see;
    }
    
    return config_.electrodes[electrode_index].gamma_see;
}

bool BoundaryManager::is_electrode_dielectric(const std::string& electrode_name) const {
    auto it = electrode_name_to_index_.find(electrode_name);
    if (it == electrode_name_to_index_.end()) {
        return false;
    }
    return is_electrode_dielectric_by_index(it->second);
}

bool BoundaryManager::is_electrode_dielectric_by_index(int electrode_index) const {
    if (!config_.use_multi_electrode || config_.electrodes.empty()) {
        return config_.dielectric_permittivity > 1.5;
    }
    
    if (electrode_index < 0 || electrode_index >= static_cast<int>(config_.electrodes.size())) {
        return false;
    }
    
    return config_.electrodes[electrode_index].is_dielectric;
}

double BoundaryManager::get_electrode_dielectric_permittivity(const std::string& electrode_name) const {
    auto it = electrode_name_to_index_.find(electrode_name);
    if (it == electrode_name_to_index_.end()) {
        return 1.0;
    }
    
    int idx = it->second;
    if (idx >= 0 && idx < static_cast<int>(config_.electrodes.size())) {
        return config_.electrodes[idx].dielectric_permittivity;
    }
    return 1.0;
}

double BoundaryManager::get_electrode_dielectric_thickness(const std::string& electrode_name) const {
    auto it = electrode_name_to_index_.find(electrode_name);
    if (it == electrode_name_to_index_.end()) {
        return 0.0;
    }
    
    int idx = it->second;
    if (idx >= 0 && idx < static_cast<int>(config_.electrodes.size())) {
        return config_.electrodes[idx].dielectric_thickness;
    }
    return 0.0;
}

int BoundaryManager::get_num_electrodes() const {
    if (config_.use_multi_electrode) {
        return config_.electrodes.size();
    }
    return 1; // Legacy mode
}

} // namespace HydroPlas
