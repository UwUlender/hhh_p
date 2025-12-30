#pragma once
#include <cmath>
#include <algorithm>

namespace HydroPlas {

// Bernoulli function B(x) = x / (exp(x) - 1)
// Robust implementation for small x
inline double Bernoulli(double x) {
    if (std::abs(x) < 1e-4) {
        double x2 = x * x;
        double x4 = x2 * x2;
        return 1.0 - x / 2.0 + x2 / 12.0 - x4 / 720.0;
    }
    return x / (std::exp(x) - 1.0);
}

// Scharfetter-Gummel Flux from cell i to i+1
// J_{i+1/2} = (D / dx) * [ B(-nu) * n_i - B(nu) * n_{i+1} ]
// where nu = (v_drift * dx) / D
// If using Potential Phi: v_drift = sgn(q) * mu * (-dPhi/dx)
// nu = sgn(q) * mu * (Phi_i - Phi_{i+1}) / D
inline double ScharfetterGummelFlux(double n_i, double n_ip1, double nu, double D, double dx) {
    return (D / dx) * (Bernoulli(-nu) * n_i - Bernoulli(nu) * n_ip1);
}

} // namespace HydroPlas
