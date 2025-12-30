#pragma once
#include "../config/ConfigParser.hpp"
#include <petsc.h>

namespace HydroPlas {

class BoundaryManager {
public:
    explicit BoundaryManager(const BoundaryConfig& config);
    
    double get_voltage(double t) const;
    double get_secondary_emission_flux(double ion_flux) const;
    
private:
    BoundaryConfig config_;
};

} // namespace HydroPlas
