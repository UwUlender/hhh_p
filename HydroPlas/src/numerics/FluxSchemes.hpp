#pragma once
#include <cmath>
#include <algorithm>

namespace HydroPlas {

inline double bernoulli(double x) {
    if (std::abs(x) < 1e-4) {
        // Taylor expansion: x / (e^x - 1) approx 1 - x/2 + x^2/12
        return 1.0 - 0.5 * x + x * x / 12.0;
    }
    return x / (std::exp(x) - 1.0);
}

// Computes flux at interface i+1/2
// n_L: density at i
// n_R: density at i+1
// D: diffusion coeff at interface
// mu: signed mobility at interface (sgn(q) * |mu|)
// dphi: phi_R - phi_L
// dx: distance between centers
inline double compute_sg_flux(double n_L, double n_R, double D, double mu, double dphi, double dx) {
    // Pe = mu * dphi / D
    // Avoid division by zero if D is tiny (unlikely in plasma except vacuum?)
    if (D < 1e-20) return 0.0; // Or pure drift?
    
    double Pe = mu * dphi / D;
    
    // Gamma = (D/dx) * [ n_L * B(Pe) - n_R * B(-Pe) ]
    return (D / dx) * (n_L * bernoulli(Pe) - n_R * bernoulli(-Pe));
}

// Computes neutral flux at interface i+1/2
// u_gas: background gas velocity at interface
inline double compute_neutral_flux(double n_L, double n_R, double D, double u_gas, double dx) {
    // Advection: Upwind
    double flux_adv = 0.0;
    if (u_gas > 0) flux_adv = u_gas * n_L;
    else flux_adv = u_gas * n_R;
    
    // Diffusion: Central difference
    double flux_diff = -D * (n_R - n_L) / dx;
    
    return flux_adv + flux_diff;
}

} // namespace HydroPlas
