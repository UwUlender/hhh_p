#include "Solver.hpp"
#include <fstream>
#include <iostream>
#include "../numerics/FluxSchemes.hpp"
#include <cmath>

namespace HydroPlas {

Solver::Solver(DM dm, const SimulationConfig& config) : dm_(dm), config_(config) {
    ctx_.dm = dm;
    ctx_.config = config;
    ctx_.boundary = new BoundaryManager(config.boundary);
    ctx_.reactions = nullptr; // Initialized after lookup table is loaded
    ctx_.output = nullptr;    // Initialized in init()
    
    // Initialize Mappings
    ctx_.idx_ne = 0;
    ctx_.idx_ni = 1;
    ctx_.idx_neps = 2;
    ctx_.idx_phi = 3;
    ctx_.idx_sigma = 4;
    ctx_.idx_excited_start = 5;
    ctx_.num_excited = config.chemistry.excited_species.size();
    
    for(int k=0; k<ctx_.num_excited; ++k) {
        ctx_.species_map[config.chemistry.excited_species[k].name] = ctx_.idx_excited_start + k;
    }
}

Solver::~Solver() {
    if (ctx_.boundary) delete ctx_.boundary;
    if (ctx_.reactions) delete ctx_.reactions;
    if (ctx_.output) delete ctx_.output;
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
    
    // Initialize ReactionHandler
    ctx_.reactions = new ReactionHandler(config_.chemistry, &ctx_.lookup);
    
    // Initialize OutputManager
    ctx_.output = new OutputManager(config_, dm_, ctx_.num_excited);
    std::vector<std::string> field_names = {"ne", "ni", "neps", "phi", "sigma"};
    for(int k=0; k<ctx_.num_excited; ++k) {
        field_names.push_back(config_.chemistry.excited_species[k].name);
    }
    ctx_.output->set_field_names(field_names);

    // 2. Setup DM fields
    ierr = DMDASetFieldName(dm_, ctx_.idx_ne, "ne"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(dm_, ctx_.idx_ni, "ni"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(dm_, ctx_.idx_neps, "neps"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(dm_, ctx_.idx_phi, "phi"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(dm_, ctx_.idx_sigma, "sigma"); CHKERRQ(ierr);
    
    for(int k=0; k<ctx_.num_excited; ++k) {
        ierr = DMDASetFieldName(dm_, ctx_.idx_excited_start + k, config_.chemistry.excited_species[k].name.c_str()); CHKERRQ(ierr);
    }

    // 3. Setup TS
    ierr = TSCreate(PETSC_COMM_WORLD, &ts_); CHKERRQ(ierr);
    ierr = TSSetDM(ts_, dm_); CHKERRQ(ierr);
    ierr = TSSetProblemType(ts_, TS_NONLINEAR); CHKERRQ(ierr);
    ierr = TSSetType(ts_, TSBDF); CHKERRQ(ierr); 
    ierr = TSSetIFunction(ts_, NULL, FormIFunction, &ctx_); CHKERRQ(ierr);
    
    // Create Matrix for Jacobian
    Mat J;
    ierr = DMCreateMatrix(dm_, &J); CHKERRQ(ierr);
    
    // Use Finite Difference Jacobian (Coloring)
    ierr = TSSetIJacobian(ts_, J, J, TSComputeIJacobianDefaultColor, &ctx_); CHKERRQ(ierr);

    ierr = MatDestroy(&J); CHKERRQ(ierr); 
    
    // Automatic time step control for RF discharges
    double dt_initial = config_.time.dt;
    if (config_.boundary.voltage_type == "RF" && config_.boundary.frequency > 0.0) {
        double T_rf = 1.0 / config_.boundary.frequency;  // RF period [s]
        double dt_rf = T_rf / 100.0;  // Resolve RF cycle with ~100 points
        if (dt_initial > dt_rf) {
            std::cout << "WARNING: Initial dt = " << dt_initial 
                     << " s is too large for RF frequency " << config_.boundary.frequency 
                     << " Hz (period = " << T_rf << " s)" << std::endl;
            std::cout << "Automatically reducing dt to " << dt_rf << " s (T/100)" << std::endl;
            dt_initial = dt_rf;
        }
    }
    
    ierr = TSSetTimeStep(ts_, dt_initial); CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts_, config_.time.t_end); CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts_, TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
    
    // Enable adaptive timestepping (optional, can be overridden by command line)
    ierr = TSSetMaxSteps(ts_, 1000000); CHKERRQ(ierr);  // Prevent infinite loops
    
    // 4. Setup PCFIELDSPLIT
    std::string split0_fields = "0,1,2,4";
    for(int k=0; k<ctx_.num_excited; ++k) {
        split0_fields += "," + std::to_string(ctx_.idx_excited_start + k);
    }
    std::string split1_fields = "3";
    
    PetscOptionsSetValue(NULL, "-pc_type", "fieldsplit");
    PetscOptionsSetValue(NULL, "-pc_fieldsplit_type", "multiplicative");
    PetscOptionsSetValue(NULL, "-pc_fieldsplit_0_fields", split0_fields.c_str());
    PetscOptionsSetValue(NULL, "-pc_fieldsplit_1_fields", split1_fields.c_str());
    
    PetscOptionsSetValue(NULL, "-fieldsplit_0_ksp_type", "gmres");
    PetscOptionsSetValue(NULL, "-fieldsplit_0_pc_type", "ilu");
    PetscOptionsSetValue(NULL, "-fieldsplit_1_ksp_type", "preonly");
    PetscOptionsSetValue(NULL, "-fieldsplit_1_pc_type", "lu"); 

    // Setup output monitoring
    ierr = TSMonitorSet(ts_, MonitorOutput, &ctx_, NULL); CHKERRQ(ierr);

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
            u[j][i][ctx_.idx_ne] = 1e14; // ne
            u[j][i][ctx_.idx_ni] = 1e14; // ni
            u[j][i][ctx_.idx_neps] = 2.0 * 1e14; // neps
            u[j][i][ctx_.idx_phi] = 0.0;  // phi
            u[j][i][ctx_.idx_sigma] = 0.0;  // sigma
            
            for(int k=0; k<ctx_.num_excited; ++k) {
                u[j][i][ctx_.idx_excited_start + k] = 1e10; // Initial seed for excited species
            }
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
    bool is_dielectric = (app->config.boundary.dielectric_permittivity > 1.5);
    double eps_d = app->config.boundary.dielectric_permittivity * eps0;
    double d_diel = app->config.boundary.dielectric_thickness;
    double V_applied = app->boundary->get_voltage(t);
    double gamma_see_total = app->config.boundary.gamma_see; // Base SEE from ions
    
    double u_gas = app->config.chemistry.gas_velocity;
    double T_gas = app->config.chemistry.gas_temperature;
    
    // Check if mass is valid to avoid div by zero
    double mass = 6.63e-26; 
    if(!app->config.chemistry.excited_species.empty()) {
         mass = app->config.chemistry.excited_species[0].mass;
    }

    for (PetscInt j = ys; j < ys + ym; j++) {
        for (PetscInt i = xs; i < xs + xm; i++) {
            
            // --- Time Derivatives ---
            f[j][i][app->idx_ne] = udot[j][i][app->idx_ne];
            f[j][i][app->idx_ni] = udot[j][i][app->idx_ni];
            f[j][i][app->idx_neps] = udot[j][i][app->idx_neps];
            f[j][i][app->idx_phi] = 0.0;
            f[j][i][app->idx_sigma] = udot[j][i][app->idx_sigma];
            
            for(int k=0; k<app->num_excited; ++k) {
                f[j][i][app->idx_excited_start + k] = udot[j][i][app->idx_excited_start + k];
            }

            // --- Local Properties ---
            double ne_val = u[j][i][app->idx_ne];
            double neps_val = u[j][i][app->idx_neps];
            double mean_energy = (ne_val > 1e-20) ? (neps_val / ne_val) : 0.0;
            if (mean_energy < 0.01) mean_energy = 0.01;

            double mu_e = app->lookup.interpolate(mean_energy, "Mobility_e");
            double D_e  = app->lookup.interpolate(mean_energy, "Diff_e");
            double mu_i = 1.0; 
            double D_i = 0.026;
            try { mu_i = app->lookup.interpolate(mean_energy, "Mobility_i"); D_i = app->lookup.interpolate(mean_energy, "Diff_i"); } catch(...){}

            // --- Source Terms (Chemistry) ---
            double N_gas = 3.22e22; // ~1 atm
            
            // Collect excited species densities
            std::vector<double> n_excited_vals(app->num_excited);
            for(int k=0; k<app->num_excited; ++k) {
                n_excited_vals[k] = u[j][i][app->idx_excited_start + k];
            }
            
            // Compute all reaction sources using ReactionHandler
            double S_ne = 0.0, S_ni = 0.0, S_neps = 0.0;
            std::vector<double> S_excited(app->num_excited, 0.0);
            
            if(app->reactions) {
                app->reactions->compute_sources(ne_val, u[j][i][app->idx_ni], n_excited_vals,
                                               mean_energy, N_gas,
                                               S_ne, S_ni, S_neps, S_excited);
            }

            // Apply source terms to residuals
            f[j][i][app->idx_ne] -= S_ne; 
            f[j][i][app->idx_ni] -= S_ni; 
            f[j][i][app->idx_neps] -= S_neps;
            
            for(int k=0; k<app->num_excited; ++k) {
                f[j][i][app->idx_excited_start + k] -= S_excited[k];
            }

            // --- Fluxes ---
            double flux_e_net = 0.0;
            double flux_i_net = 0.0;
            double flux_eps_net = 0.0;
            double heating_net = 0.0;
            
            std::vector<double> flux_excited_net(app->num_excited, 0.0);

            // --- Right Face ---
            if (i < M-1) {
                double dphi = u[j][i+1][app->idx_phi] - u[j][i][app->idx_phi];
                
                // Charged
                double nu_e = mu_e * dphi / D_e;
                double flux_e = ScharfetterGummelFlux(u[j][i][app->idx_ne], u[j][i+1][app->idx_ne], nu_e, D_e, dx);
                
                double nu_i = -mu_i * dphi / D_i;
                double flux_i = ScharfetterGummelFlux(u[j][i][app->idx_ni], u[j][i+1][app->idx_ni], nu_i, D_i, dx);

                double mu_eps = (5.0/3.0) * mu_e;
                double D_eps = (5.0/3.0) * D_e;
                double nu_eps = mu_eps * dphi / D_eps;
                double flux_eps = ScharfetterGummelFlux(u[j][i][app->idx_neps], u[j][i+1][app->idx_neps], nu_eps, D_eps, dx);
                
                flux_e_net += flux_e;
                flux_i_net += flux_i;
                flux_eps_net += flux_eps;
                heating_net += flux_e * (dphi/dx);
                
                // Neutral
                for(int k=0; k<app->num_excited; ++k) {
                    double D_k = app->config.chemistry.excited_species[k].diffusion_coeff;
                    double nu_k = u_gas * dx / D_k;
                    PetscInt idx = app->idx_excited_start + k;
                    double flux_k = ScharfetterGummelFlux(u[j][i][idx], u[j][i+1][idx], nu_k, D_k, dx);
                    flux_excited_net[k] += flux_k;
                }

            } else {
                 // Right Boundary (Wall)
                 double E_wall = u[j][i][app->idx_phi] / (0.5*dx);
                 
                 double flux_i_out = mu_i * u[j][i][app->idx_ni] * E_wall; 
                 if (flux_i_out < 0) flux_i_out = 0;
                 
                 double v_th_e = sqrt(8.0 * mean_energy * 1.6e-19 / (3.14 * 9.11e-31)); 
                 double flux_e_out = 0.25 * u[j][i][app->idx_ne] * v_th_e;
                 if (-mu_e * E_wall > 0) flux_e_out += -mu_e * u[j][i][app->idx_ne] * E_wall; 
                 
                 flux_e_net += flux_e_out;
                 flux_i_net += flux_i_out;
                 flux_eps_net += flux_e_out * mean_energy;
                 
                 for(int k=0; k<app->num_excited; ++k) {
                     double gamma_k = app->config.chemistry.excited_species[k].wall_quenching_prob;
                     double m_k = app->config.chemistry.excited_species[k].mass;
                     double v_th_k = sqrt(8.0 * 1.38e-23 * T_gas / (3.14 * m_k));
                     PetscInt idx = app->idx_excited_start + k;
                     double flux_k_out = (gamma_k * v_th_k / 4.0) * u[j][i][idx];
                     if(u_gas > 0) flux_k_out += u_gas * u[j][i][idx];
                     flux_excited_net[k] += flux_k_out;
                 }
            }

            // --- Left Face ---
            if (i > 0) {
                 double dphi = u[j][i][app->idx_phi] - u[j][i-1][app->idx_phi];
                 
                 double nu_e = mu_e * dphi / D_e;
                 double flux_e = ScharfetterGummelFlux(u[j][i-1][app->idx_ne], u[j][i][app->idx_ne], nu_e, D_e, dx);
                 
                 double nu_i = -mu_i * dphi / D_i;
                 double flux_i = ScharfetterGummelFlux(u[j][i-1][app->idx_ni], u[j][i][app->idx_ni], nu_i, D_i, dx);

                 double mu_eps = (5.0/3.0) * mu_e;
                 double D_eps = (5.0/3.0) * D_e;
                 double nu_eps = mu_eps * dphi / D_eps;
                 double flux_eps = ScharfetterGummelFlux(u[j][i-1][app->idx_neps], u[j][i][app->idx_neps], nu_eps, D_eps, dx);
                 
                 flux_e_net -= flux_e;
                 flux_i_net -= flux_i;
                 flux_eps_net -= flux_eps;
                 heating_net += flux_e * (dphi/dx);
                 
                 for(int k=0; k<app->num_excited; ++k) {
                    double D_k = app->config.chemistry.excited_species[k].diffusion_coeff;
                    double nu_k = u_gas * dx / D_k;
                    PetscInt idx = app->idx_excited_start + k;
                    double flux_k = ScharfetterGummelFlux(u[j][i-1][idx], u[j][i][idx], nu_k, D_k, dx);
                    flux_excited_net[k] -= flux_k;
                }

            } else {
                 // Left Boundary (i=0)
                 double Phi_wall = V_applied;
                 double dphi_wall = u[j][i][app->idx_phi] - Phi_wall;
                 double E_wall_l = -dphi_wall / (0.5*dx); 
                 
                 double flux_i_boundary = 0.0;
                 if (E_wall_l < 0) { 
                     flux_i_boundary = mu_i * u[j][i][app->idx_ni] * E_wall_l; 
                 }
                 
                 double flux_e_boundary = -gamma_see_total * flux_i_boundary;
                 
                 for(int k=0; k<app->num_excited; ++k) {
                     double gamma_see_k = app->config.chemistry.excited_species[k].wall_see_prob;
                     if (gamma_see_k > 0) {
                        double gamma_k = app->config.chemistry.excited_species[k].wall_quenching_prob;
                        double m_k = app->config.chemistry.excited_species[k].mass;
                        double v_th_k = sqrt(8.0 * 1.38e-23 * T_gas / (3.14 * m_k));
                        PetscInt idx = app->idx_excited_start + k;
                        double n_val = u[j][i][idx];
                        double flux_n_in = (gamma_k * v_th_k / 4.0) * n_val;
                        if(u_gas < 0) flux_n_in += (-u_gas) * n_val;
                        
                        flux_e_boundary += gamma_see_k * flux_n_in;
                     }
                 }

                 if (E_wall_l > 0) {
                      double v_th_e = sqrt(8.0 * mean_energy * 1.6e-19 / (3.14 * 9.11e-31));
                      flux_e_boundary -= 0.25 * u[j][i][app->idx_ne] * v_th_e;
                 }
                 
                 flux_e_net -= flux_e_boundary;
                 flux_i_net -= flux_i_boundary;
                 flux_eps_net -= flux_e_boundary * mean_energy; 
                 
                 double J_wall = -q * (flux_i_boundary - flux_e_boundary);
                 if (is_dielectric) {
                     f[j][i][app->idx_sigma] = udot[j][i][app->idx_sigma] - J_wall;
                 } else {
                     f[j][i][app->idx_sigma] = u[j][i][app->idx_sigma];
                 }
                 
                 for(int k=0; k<app->num_excited; ++k) {
                     double gamma_k = app->config.chemistry.excited_species[k].wall_quenching_prob;
                     double m_k = app->config.chemistry.excited_species[k].mass;
                     double v_th_k = sqrt(8.0 * 1.38e-23 * T_gas / (3.14 * m_k));
                     PetscInt idx = app->idx_excited_start + k;
                     double flux_k_boundary = -(gamma_k * v_th_k / 4.0) * u[j][i][idx];
                     if(u_gas < 0) flux_k_boundary += u_gas * u[j][i][idx];
                     
                     flux_excited_net[k] -= flux_k_boundary;
                 }
            }
            
            f[j][i][app->idx_ne] += flux_e_net / dx;
            f[j][i][app->idx_ni] += flux_i_net / dx;
            f[j][i][app->idx_neps] += flux_eps_net / dx;
            f[j][i][app->idx_neps] -= heating_net / dx;
            
            for(int k=0; k<app->num_excited; ++k) {
                f[j][i][app->idx_excited_start + k] += flux_excited_net[k] / dx;
            }

            // --- Poisson ---
            if (i == 0) {
                if (is_dielectric) {
                    double phi_s = u[j][i][app->idx_phi];
                    double phi_1 = u[j][i+1][app->idx_phi];
                    double E_p = -(phi_1 - phi_s) / dx;
                    double E_d = (V_applied - phi_s) / d_diel;
                    double sigma = u[j][i][app->idx_sigma];
                    f[j][i][app->idx_phi] = eps0 * E_p - eps_d * E_d - sigma;
                } else {
                    f[j][i][app->idx_phi] = u[j][i][app->idx_phi] - V_applied;
                }
            } else if (i == M-1) {
                f[j][i][app->idx_phi] = u[j][i][app->idx_phi]; 
            } else {
                double rho = q * (u[j][i][app->idx_ni] - u[j][i][app->idx_ne]);
                double d2phi_dx2 = (u[j][i+1][app->idx_phi] - 2*u[j][i][app->idx_phi] + u[j][i-1][app->idx_phi]) / (dx*dx);
                f[j][i][app->idx_phi] = - eps * d2phi_dx2 - rho;
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

PetscErrorCode MonitorOutput(TS ts, PetscInt step, PetscReal time, Vec U, void* ctx) {
    AppCtx* app = (AppCtx*)ctx;
    PetscErrorCode ierr;
    
    // Write output at specified intervals
    if (step % app->config.time.output_interval == 0) {
        if (app->output) {
            ierr = app->output->write_output(U, time, step); CHKERRQ(ierr);
        }
        PetscPrintf(PETSC_COMM_WORLD, "Step %d, Time %.3e s\n", step, time);
    }
    
    return 0;
}

} // namespace HydroPlas
