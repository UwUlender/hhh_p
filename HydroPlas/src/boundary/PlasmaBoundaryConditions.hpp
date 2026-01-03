#pragma once
#include <vector>
#include <string>
#include <memory>
#include <cmath>

namespace HydroPlas {

// Forward declarations
class RectilinearGrid;
class Species;

/**
 * Abstract Base Class for Boundary Conditions
 * Mimics the modularity of Zapdos/MOOSE BC system
 */
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    // Returns the contribution to the residual for a specific cell/face on the boundary
    // x_boundary: density/value at the boundary cell
    // x_internal: density/value at the adjacent internal cell
    // dx: distance between centers or to face
    // area: face area
    // species_idx: index of species (if applicable)
    // variable_idx: index in the state vector
    virtual double compute_flux(double x_boundary, double x_internal, 
                              double dx, double area, 
                              const Species* species, 
                              double E_field_normal,
                              double mean_energy) const = 0;

    virtual std::string get_type() const = 0;
};

/**
 * Hagelaar Electron Boundary Condition
 * Based on: G.J.M. Hagelaar et al., Phys. Rev. E 62, 1452 (2000)
 * 
 * Gamma = (1-r)/(1+r) * [ drift_term + 0.25 * v_th * n ]
 */
class HagelaarElectronBC : public BoundaryCondition {
public:
    HagelaarElectronBC(double reflection_coeff = 0.0) 
        : r_(reflection_coeff) {}

    double compute_flux(double n_wall, double n_internal, 
                      double dx, double area, 
                      const Species* species, 
                      double E_field_normal,
                      double mean_energy) const override {
        
        // Constants
        constexpr double k_B = 1.380649e-23;
        constexpr double m_e = 9.10938356e-31;
        constexpr double pi = 3.14159265359;
        constexpr double q_e = 1.6021766e-19;

        // Thermal velocity: v_th = sqrt(8 * k_B * T_e / (pi * m_e))
        // mean_energy (units often eV or similar in code, assume eV for now)
        // T_e (K) = mean_energy * 2/3 * 11604 
        // Or directly: T_e (eV) = 2/3 * mean_energy
        double T_e_eV = (2.0/3.0) * mean_energy;
        double T_e_K = T_e_eV * 11604.5; // eV to K
        
        double v_th = std::sqrt(8.0 * k_B * T_e_K / (pi * m_e));

        // Drift term: only if electric field pushes electrons to wall
        // E_field_normal points OUT of domain? Need convention.
        // Assuming E_field_normal is E dot n. 
        // Force F = -eE. If E points to wall, Force points away.
        // If E points away from wall (positive), electrons attracted to wall.
        
        double mu = 0.0;
        double D = 0.0;
        if (species) {
            species->get_transport(mean_energy, mu, D);
        }

        double drift_vel = 0.0;
        // If field attracts electrons to wall
        // Electron charge is negative. Force is opposite to field.
        // If E_field_normal > 0 (pointing out), Force is IN.
        // If E_field_normal < 0 (pointing in), Force is OUT (towards wall).
        
        // Logic from Zapdos:
        // if (normals * -1.0 * E > 0.0) _a = 1.0 (drift active)
        // normals points OUT. 
        // If (-E) points OUT, then E points IN.
        // So if E points IN (negative), electrons (negative) are pushed OUT.
        
        double directed_drift = 0.0;
        if (E_field_normal < 0.0) { // Field points INTO domain, e- pushed OUT
             directed_drift = std::abs(mu * E_field_normal);
        }

        // Total directed velocity
        // Gamma = (1-r)/(1+r) * (drift + 0.25*v_th) * n
        double factor = (1.0 - r_) / (1.0 + r_);
        double flux = factor * (directed_drift * n_wall + 0.25 * v_th * n_wall);
        
        return flux * area; // Returns total particle/sec across face
    }

    std::string get_type() const override { return "Hagelaar"; }

private:
    double r_; // Reflection coefficient
};

} // namespace HydroPlas
