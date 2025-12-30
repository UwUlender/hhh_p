#include "Solver.hpp"
#include <iostream>
#include "../numerics/FluxSchemes.hpp"

namespace HydroPlas {

Solver::Solver(DM dm, const SimulationConfig& config) : dm_(dm), config_(config) {
    ctx_.dm = dm;
    ctx_.config = config;
    ctx_.boundary = new BoundaryManager(config.boundary);
}

Solver::~Solver() {
    if (ctx_.boundary) delete ctx_.boundary;
    if (ts_) TSDestroy(&ts_);
}

PetscErrorCode Solver::init() {
    PetscErrorCode ierr;
    
    // 1. Load Chemistry Data
    if (config_.chemistry.mode == "Pre-calculated") {
        try {
            // Check if file exists, if not, create a dummy one
            std::ifstream f(config_.chemistry.transport_table_file);
            if (!f.good()) {
                 std::cout << "Warning: Transport table not found. Creating dummy at " << config_.chemistry.transport_table_file << std::endl;
                 std::ofstream dummy(config_.chemistry.transport_table_file);
                 dummy << "Energy Mobility_e Diff_e Mobility_i Diff_i Rate_Ionization Rate_Excitation" << std::endl;
                 dummy << "0.01 100.0 10.0 1.0 0.1 0.0 0.0" << std::endl;
                 dummy << "1.0 50.0 5.0 0.5 0.05 1e-16 1e-17" << std::endl;
                 dummy << "10.0 10.0 1.0 0.1 0.01 1e-14 1e-15" << std::endl;
                 dummy << "100.0 1.0 0.1 0.01 0.001 1e-13 1e-14" << std::endl;
                 dummy.close();
            }
            ctx_.lookup.load_from_file(config_.chemistry.transport_table_file);
        } catch (const std::exception& e) {
            std::cerr << "Error loading chemistry: " << e.what() << std::endl;
            return -1;
        }
    }

    // 2. Setup DM fields (5 DOFs)
    ierr = DMDASetFieldName(dm_, 0, "ne"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(dm_, 1, "ni"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(dm_, 2, "neps"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(dm_, 3, "phi"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(dm_, 4, "sigma"); CHKERRQ(ierr); // Surface Charge

    // 3. Setup TS
    ierr = TSCreate(PETSC_COMM_WORLD, &ts_); CHKERRQ(ierr);
    ierr = TSSetDM(ts_, dm_); CHKERRQ(ierr);
    ierr = TSSetProblemType(ts_, TS_NONLINEAR); CHKERRQ(ierr);
    ierr = TSSetType(ts_, TSBDF); CHKERRQ(ierr); 
    ierr = TSSetIFunction(ts_, NULL, FormIFunction, &ctx_); CHKERRQ(ierr);
    
    // Create Matrix for Jacobian
    Mat J;
    ierr = DMCreateMatrix(dm_, &J); CHKERRQ(ierr);
    ierr = TSSetIJacobian(ts_, J, J, FormIJacobian, &ctx_); CHKERRQ(ierr);
    ierr = MatDestroy(&J); CHKERRQ(ierr); 
    
    ierr = TSSetTimeStep(ts_, config_.time.dt); CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts_, config_.time.t_end); CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts_, TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
    
    // 4. Setup PCFIELDSPLIT
    // Split 0: Transport+Sigma (0,1,2,4)
    // Split 1: Poisson (3)
    PetscOptionsSetValue(NULL, "-pc_type", "fieldsplit");
    PetscOptionsSetValue(NULL, "-pc_fieldsplit_type", "multiplicative");
    PetscOptionsSetValue(NULL, "-pc_fieldsplit_0_fields", "0,1,2,4");
    PetscOptionsSetValue(NULL, "-pc_fieldsplit_1_fields", "3");
    
    PetscOptionsSetValue(NULL, "-fieldsplit_0_ksp_type", "gmres");
    PetscOptionsSetValue(NULL, "-fieldsplit_0_pc_type", "ilu");
    PetscOptionsSetValue(NULL, "-fieldsplit_1_ksp_type", "preonly");
    PetscOptionsSetValue(NULL, "-fieldsplit_1_pc_type", "lu"); 

    ierr = TSSetFromOptions(ts_); CHKERRQ(ierr);
    
    return 0;
}

PetscErrorCode Solver::solve() {
    PetscErrorCode ierr;
    Vec U;
    ierr = DMCreateGlobalVector(dm_, &U); CHKERRQ(ierr);
    
    // Initial Conditions
    PetscScalar ***u;
    ierr = DMDAVecGetArrayDOF(dm_, U, &u); CHKERRQ(ierr);
    PetscInt xs, ys, xm, ym;
    ierr = DMDAGetCorners(dm_, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
    
    for (PetscInt j = ys; j < ys + ym; j++) {
        for (PetscInt i = xs; i < xs + xm; i++) {
            u[j][i][0] = 1e14; // ne
            u[j][i][1] = 1e14; // ni
            u[j][i][2] = 2.0 * 1e14; // neps
            u[j][i][3] = 0.0;  // phi
            u[j][i][4] = 0.0;  // sigma
        }
    }
    ierr = DMDAVecRestoreArrayDOF(dm_, U, &u); CHKERRQ(ierr);
    
    ierr = TSSolve(ts_, U); CHKERRQ(ierr);
    
    ierr = VecDestroy(&U); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode FormIFunction(TS ts, PetscReal t, Vec U, Vec Udot, Vec F, void* ctx) {
    AppCtx* app = (AppCtx*)ctx;
    PetscErrorCode ierr;
    DM dm = app->dm;
    
    Vec locU;
    ierr = DMGetLocalVector(dm, &locU); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dm, U, INSERT_VALUES, locU); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dm, U, INSERT_VALUES, locU); CHKERRQ(ierr);

    PetscInt xs, ys, xm, ym;
    PetscInt M, N;
    ierr = DMDAGetInfo(dm, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
    ierr = DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
    
    const PetscInt NE = 0, NI = 1, NEPS = 2, PHI = 3, SIGMA = 4;
    
    PetscScalar ***u, ***udot, ***f;
    
    ierr = DMDAVecGetArrayDOF(dm, locU, &u); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(dm, Udot, &udot); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(dm, F, &f); CHKERRQ(ierr);
    
    PetscReal Lx = app->config.domain.Lx;
    PetscReal dx = Lx / (M > 1 ? M-1 : 1);
    
    const double q = 1.602e-19;
    const double eps0 = 8.854e-12;
    const double eps = eps0; 
    
    // Boundary parameters
    bool is_dielectric = (app->config.boundary.dielectric_permittivity > 1.5); // Heuristic check
    double eps_d = app->config.boundary.dielectric_permittivity * eps0;
    double d_diel = app->config.boundary.dielectric_thickness;
    double V_applied = app->boundary->get_voltage(t);
    double gamma_see = app->config.boundary.gamma_see;

    for (PetscInt j = ys; j < ys + ym; j++) {
        for (PetscInt i = xs; i < xs + xm; i++) {
            
            // Time derivatives
            f[j][i][NE] = udot[j][i][NE];
            f[j][i][NI] = udot[j][i][NI];
            f[j][i][NEPS] = udot[j][i][NEPS];
            f[j][i][PHI] = 0.0;
            f[j][i][SIGMA] = udot[j][i][SIGMA]; // Default d(sigma)/dt = 0 in bulk

            // Transport Properties Lookups (Needed for Fluxes)
            double ne_val = u[j][i][NE];
            double neps_val = u[j][i][NEPS];
            double mean_energy = (ne_val > 1e-20) ? (neps_val / ne_val) : 0.0;
            if (mean_energy < 0.01) mean_energy = 0.01;

            double mu_e = app->lookup.interpolate(mean_energy, "Mobility_e");
            double D_e  = app->lookup.interpolate(mean_energy, "Diff_e");
            double mu_i = 1.0; 
            double D_i = 0.026;
            try { mu_i = app->lookup.interpolate(mean_energy, "Mobility_i"); D_i = app->lookup.interpolate(mean_energy, "Diff_i"); } catch(...){}

            // Source Terms
            double R_ion = 0.0;
            try { R_ion = app->lookup.interpolate(mean_energy, "Rate_Ionization"); } catch(...) {}
            double N_gas = 3.22e22; 
            double S_ion = R_ion * N_gas * ne_val;
            
            f[j][i][NE] -= S_ion; 
            f[j][i][NI] -= S_ion; 

             // Energy Loss
            double EnergyLoss = S_ion * 15.76;
            f[j][i][NEPS] += EnergyLoss; // Loss term, so + because F = udot - S + divG. S is negative loss. Wait.
            // Eq: dn/dt + div = S. F = dn/dt + div - S.
            // S_ion is positive source. So - S_ion.
            // Energy eq: ... = -Loss. So S_eps = -Loss.
            // F = ... - (-Loss) = ... + Loss. Correct.
            
            // --- Fluxes ---
            double flux_e_net = 0.0;
            double flux_i_net = 0.0;
            double flux_eps_net = 0.0;
            double heating_net = 0.0;

            // --- Right Face ---
            if (i < M-1) {
                double dphi = u[j][i+1][PHI] - u[j][i][PHI];
                
                // Electron
                double nu_e = mu_e * dphi / D_e;
                double flux_e = ScharfetterGummelFlux(u[j][i][NE], u[j][i+1][NE], nu_e, D_e, dx);
                
                // Ion
                double nu_i = -mu_i * dphi / D_i;
                double flux_i = ScharfetterGummelFlux(u[j][i][NI], u[j][i+1][NI], nu_i, D_i, dx);

                // Energy
                double mu_eps = (5.0/3.0) * mu_e;
                double D_eps = (5.0/3.0) * D_e;
                double nu_eps = mu_eps * dphi / D_eps;
                double flux_eps = ScharfetterGummelFlux(u[j][i][NEPS], u[j][i+1][NEPS], nu_eps, D_eps, dx);
                
                flux_e_net += flux_e;
                flux_i_net += flux_i;
                flux_eps_net += flux_eps;
                
                // Heating contribution (Right half)
                // -e * Gamma_e * E. E = -dphi/dx. -> e * Gamma_e * dphi/dx
                // Heating approx at face: flux_e * dphi/dx
                heating_net += flux_e * (dphi/dx); // Contribution to integral
            } else {
                 // Right Boundary (Wall)
                 // Fluxes out to wall (i -> Wall)
                 
                 // E_wall approx
                 double dphi_wall = u[j][i][PHI] - 0.0; // Grounded Right Wall?
                 // Or use boundary condition logic.
                 
                 // If Right is Grounded Metal
                 // E ~ phi[i]/(dx/2). 
                 double E_wall = u[j][i][PHI] / (0.5*dx);
                 
                 // Drift Fluxes
                 // Ion out: mu_i * n * E_wall.
                 double flux_i_out = mu_i * u[j][i][NI] * E_wall; 
                 // If E_wall > 0 (Field points Right), Ions go Right. 
                 if (flux_i_out < 0) flux_i_out = 0; // Ions don't come from wall? 
                 
                 // Electron out
                 double flux_e_out = -mu_e * u[j][i][NE] * E_wall; // Electrons go Left if E>0.
                 // Actually thermal velocity dominates at wall.
                 double v_th_e = sqrt(8.0 * mean_energy * 1.6e-19 / (3.14 * 9.11e-31)); 
                 flux_e_out += 0.25 * u[j][i][NE] * v_th_e; // Thermal flux
                 
                 flux_e_net += flux_e_out;
                 flux_i_net += flux_i_out;
                 flux_eps_net += flux_e_out * mean_energy; // Convective energy loss
            }

            // --- Left Face ---
            if (i > 0) {
                 // Already handled by neighbor's Right Face in loop?
                 // No, standard loop calculates divergence.
                 // Need flux from i-1 to i.
                 
                 double dphi = u[j][i][PHI] - u[j][i-1][PHI];
                 // Use neighbor's props? Or average? SG uses local props approx.
                 // For consistent SG, use average or upstream.
                 // Simplification: Use local props for coefficient (ok for small gradients).
                 
                 double nu_e = mu_e * dphi / D_e;
                 double flux_e = ScharfetterGummelFlux(u[j][i-1][NE], u[j][i][NE], nu_e, D_e, dx);
                 
                 double nu_i = -mu_i * dphi / D_i;
                 double flux_i = ScharfetterGummelFlux(u[j][i-1][NI], u[j][i][NI], nu_i, D_i, dx);

                 double mu_eps = (5.0/3.0) * mu_e;
                 double D_eps = (5.0/3.0) * D_e;
                 double nu_eps = mu_eps * dphi / D_eps;
                 double flux_eps = ScharfetterGummelFlux(u[j][i-1][NEPS], u[j][i][NEPS], nu_eps, D_eps, dx);
                 
                 flux_e_net -= flux_e;
                 flux_i_net -= flux_i;
                 flux_eps_net -= flux_eps;
                 
                 heating_net += flux_e * (dphi/dx);

            } else {
                 // Left Boundary (i=0) - Wall
                 // Fluxes from Wall to i=0.
                 // Wall is at left. Flux_in is Gamma(LeftWall).
                 
                 // Potential at Wall: V_applied (if Metal) or Phi_surface (if Dielectric)
                 double Phi_wall = V_applied;
                 if (is_dielectric) {
                     // E_d = (V_applied - u[i][PHI]) / (dx/2 + d_diel/eps_r)?
                     // Assume Phi(i) is center.
                 }
                 
                 double dphi_wall = u[j][i][PHI] - Phi_wall;
                 double E_wall_l = -dphi_wall / (0.5*dx); // Field pointing Right
                 
                 // Ion Flux from Wall (should be 0 unless reflected)
                 // Ion Flux TO Wall (going Left).
                 // Flux_i_left = - mu_i * n * E.
                 // If E points Right (positive), Ions go Right (away from wall).
                 // If E points Left (negative), Ions go Left (to wall).
                 
                 double flux_i_boundary = 0.0;
                 if (E_wall_l < 0) { // Field pulls ions to wall
                     flux_i_boundary = mu_i * u[j][i][NI] * E_wall_l; // Negative value
                 }
                 
                 // Electron Flux
                 // SEE: Gamma_e = - gamma * Gamma_i
                 // Gamma_e_boundary = -gamma * flux_i_boundary (positive value)
                 double flux_e_boundary = -gamma_see * flux_i_boundary;
                 
                 // Add Thermal component if E pushes electrons to wall (E > 0)
                 // Gamma_thermal = -0.25 * n * v_th
                 if (E_wall_l > 0) {
                      double v_th_e = sqrt(8.0 * mean_energy * 1.6e-19 / (3.14 * 9.11e-31));
                      flux_e_boundary -= 0.25 * u[j][i][NE] * v_th_e;
                 }
                 
                 // Net flux leaving cell i to left (negative of flux entering)
                 // Flux_net += Flux_right - Flux_left.
                 // Flux_left is flux at interface i-1/2.
                 // Here it is flux_e_boundary.
                 
                 flux_e_net -= flux_e_boundary;
                 flux_i_net -= flux_i_boundary;
                 flux_eps_net -= flux_e_boundary * mean_energy; // Approx energy of emitted electrons? Usually low.
                 
                 // Surface Charge Evolution (Left Wall)
                 if (is_dielectric) {
                     // d(sigma)/dt = J_plasma
                     // J_plasma = q * (Gamma_i - Gamma_e)
                     // Current TO wall (going left) = -q * (flux_i_boundary - flux_e_boundary)
                     // Sigma is charge on dielectric surface facing plasma.
                     // Accumulation: J_in.
                     double J_wall = -q * (flux_i_boundary - flux_e_boundary);
                     f[j][i][SIGMA] = udot[j][i][SIGMA] - J_wall;
                 } else {
                     f[j][i][SIGMA] = u[j][i][SIGMA]; // Dummy fix to 0
                 }
            }
            
            f[j][i][NE] += flux_e_net / dx;
            f[j][i][NI] += flux_i_net / dx;
            f[j][i][NEPS] += flux_eps_net / dx;
            
            // Add Heating
            // Heating term should be -e * Gamma * E.
            // approximated above.
            // Need to match units.
            // Heating term is power density W/m3?
            // If neps is eV/m3. Power is eV/m3s.
            // JouleHeating above was Flux * V/m.
            // Flux (1/m2s) * V/m = V/m3s = eV/m3s.
            // Sign: Gamma_e * GradPhi.
            f[j][i][NEPS] -= heating_net / dx; // Divergence-like sum

            // --- Poisson ---
            // Interior: -eps * d2phi - rho = 0
            // Boundary:
            // Left (i=0):
            // If Metal: Phi = V(t) -> f = Phi - V
            // If Dielectric:
            //   (eps_d * (V - Phi_s)/d - eps_0 * (Phi_s - Phi_1)/dx) = sigma
            //   Approximation using Ghost Point?
            //   Let's use the Algebraic equation at node 0 for Phi.
            
            if (i == 0) {
                if (is_dielectric) {
                    double phi_s = u[j][i][PHI];
                    double phi_1 = u[j][i+1][PHI];
                    double E_p = -(phi_1 - phi_s) / dx;
                    double E_d = (V_applied - phi_s) / d_diel;
                    
                    // Gauss: eps_d * E_d - eps_0 * E_p = sigma
                    double sigma = u[j][i][SIGMA];
                    f[j][i][PHI] = eps_d * E_d - eps0 * E_p - sigma;
                } else {
                    f[j][i][PHI] = u[j][i][PHI] - V_applied;
                }
            } else if (i == M-1) {
                f[j][i][PHI] = u[j][i][PHI]; // Ground
            } else {
                double rho = q * (u[j][i][NI] - u[j][i][NE]);
                double d2phi_dx2 = (u[j][i+1][PHI] - 2*u[j][i][PHI] + u[j][i-1][PHI]) / (dx*dx);
                f[j][i][PHI] = - eps * d2phi_dx2 - rho;
            }
        }
    }
    
    ierr = DMDAVecRestoreArrayDOF(dm, locU, &u); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(dm, Udot, &udot); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(dm, F, &f); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &locU); CHKERRQ(ierr);
    
    return 0;
}

PetscErrorCode FormIJacobian(TS ts, PetscReal t, Vec U, Vec Udot, PetscReal shift, Mat J, Mat P, void* ctx) {
    PetscErrorCode ierr;
    ierr = MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    if (J != P) {
        ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    }
    return 0;
}

} // namespace HydroPlas
