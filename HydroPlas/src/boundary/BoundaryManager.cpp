#include "BoundaryManager.hpp"
#include <cmath>

namespace HydroPlas {

BoundaryManager::BoundaryManager(const BoundaryConfig& config) : config_(config) {}

double BoundaryManager::get_voltage(double t) const {
    if (config_.voltage_type == "DC") {
        return config_.voltage_amplitude;
    } else if (config_.voltage_type == "RF") {
        return config_.voltage_amplitude * std::sin(2.0 * M_PI * config_.frequency * t) + config_.bias;
    }
    return 0.0;
}

double BoundaryManager::get_secondary_emission_flux(double ion_flux) const {
    // If ion flux is towards the wall (positive if flow is out), electron flux is gamma * ion_flux into the plasma (negative relative to wall normal?)
    // Depends on sign convention.
    return config_.gamma_see * ion_flux;
}

} // namespace HydroPlas
