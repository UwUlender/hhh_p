#include "PlasmaSolver.hpp"
#include "../numerics/FluxSchemes.hpp"
#include "../io/OutputManager.hpp" 
#include <iostream>
#include <cmath>
#include <mpi.h>

namespace HydroPlas {

PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void* ctx_void);
PetscErrorCode FormJacobian(SNES snes, Vec X, Mat J, Mat P, void* ctx_void);

PlasmaSolver::PlasmaSolver(RectilinearGrid& grid, Chemistry& chemistry, BoundaryManager& boundary, SimulationConfig& config)
    : grid_(grid), chemistry_(chemistry), boundary_(boundary), config_(config) {
    ctx_.grid = &grid;
    ctx_.chemistry = &chemistry;
    ctx_.boundary = &boundary;
    ctx_.config = &config;
    ctx_.idx_n_eps = chemistry.get_num_species();
    ctx_.idx_phi = chemistry.get_num_species() + 1;
    ctx_.X_prev = nullptr;
    ctx_.dt = 1e-12;
    ctx_.time = 0.0;
    X_ = nullptr;
    F_ = nullptr;
    snes_ = nullptr;
}

PlasmaSolver::~PlasmaSolver() {
    if (ctx_.X_prev) VecDestroy(&ctx_.X_prev);
    if (snes_) SNESDestroy(&snes_);
    if (X_) VecDestroy(&X_);
    if (F_) VecDestroy(&F_);
}

void PlasmaSolver::initialize() {
    setup_dofs();
    setup_solver();
}

void PlasmaSolver::setup_dofs() {
    int num_species = chemistry_.get_num_species();
    int dofs = num_species + 2; // Species + Energy + Potential
    grid_.initialize(dofs);
    
    DMCreateGlobalVector(grid_.get_dm(), &X_);
    VecDuplicate(X_, &F_);
    VecDuplicate(X_, &ctx_.X_prev);
    
    // Initialize X
    // Set initial densities to small value, potential to 0
    PetscScalar *x;
    VecGetArray(X_, &x);
    int localsize;
    VecGetLocalSize(X_, &localsize);
    for (int i=0; i<localsize; ++i) {
        if (i % dofs == ctx_.idx_phi) x[i] = 0.0;
        else x[i] = 1e14; // Background plasma
    }
    VecRestoreArray(X_, &x);
}

void PlasmaSolver::setup_solver() {
    SNESCreate(PETSC_COMM_WORLD, &snes_);
    SNESSetDM(snes_, grid_.get_dm());
    
    SNESSetFunction(snes_, F_, FormFunction, &ctx_);
    
    Mat J;
    MatCreateSNESMF(snes_, &J);
    DMCreateMatrix(grid_.get_dm(), &P_poisson_); 
    SNESSetJacobian(snes_, J, P_poisson_, FormJacobian, &ctx_);
    MatDestroy(&J);
    
    // Configure JFNK and Preconditioner
    KSP ksp;
    SNESGetKSP(snes_, &ksp);
    
    // Set KSP type from config
    KSPSetType(ksp, config_.solver.ksp_type.c_str());
    
    // Set PC type from config
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, config_.solver.preconditioner.c_str());
    
    // Set tolerances
    SNESSetTolerances(snes_, config_.solver.tolerance, PETSC_DEFAULT, PETSC_DEFAULT, 
                      config_.solver.max_iterations, PETSC_DEFAULT);
    
    SNESSetFromOptions(snes_);
}

void PlasmaSolver::solve_step(double dt, double time) {
    ctx_.dt = dt;
    ctx_.time = time;
    
    // Store previous state
    VecCopy(X_, ctx_.X_prev);
    
    // Solve
    SNESSolve(snes_, NULL, X_);
}

void PlasmaSolver::save_state(const std::string& filename, int step, double time) {
    // Use binary viewer to save X_
    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer);
    VecView(X_, viewer);
    PetscViewerDestroy(&viewer);
}

void PlasmaSolver::save_rates(OutputManager& output, int step) {
    // Compute rates
    int nx = grid_.get_nx();
    int ny = grid_.get_ny();
    int num_reactions = chemistry_.get_reactions().size();
    
    // Prepare data structure: rates[reaction_idx][cell_idx]
    std::vector<std::vector<double>> rates_data(num_reactions, std::vector<double>(nx*ny, 0.0));
    
    // Get local vector
    Vec Xloc;
    DMGetLocalVector(grid_.get_dm(), &Xloc);
    DMGlobalToLocalBegin(grid_.get_dm(), X_, INSERT_VALUES, Xloc);
    DMGlobalToLocalEnd(grid_.get_dm(), X_, INSERT_VALUES, Xloc);
    
    PetscScalar ***x;
    DMDAVecGetArrayDOF(grid_.get_dm(), Xloc, &x);
    
    int xs, ys, xm, ym;
    DMDAGetCorners(grid_.get_dm(), &xs, &ys, NULL, &xm, &ym, NULL);
    
    int num_species = chemistry_.get_num_species();
    int idx_eps = ctx_.idx_n_eps;
    
    const auto& reactions = chemistry_.get_reactions();
    const auto& species = chemistry_.get_species();
    
    // Compute rates for local cells
    for (int j = ys; j < ys + ym; ++j) {
        for (int i = xs; i < xs + xm; ++i) {
            // Gather densities
            std::vector<double> densities(num_species);
            for (int k = 0; k < num_species; ++k) {
                densities[k] = x[j][i][k];
            }
            
            // Calculate mean energy
            double n_e = 0.0;
            for (int k = 0; k < num_species; ++k) {
                if (species[k].type == SpeciesType::Electron) {
                    n_e = densities[k];
                    break;
                }
            }
            
            double n_eps = x[j][i][idx_eps];
            double mean_en = (n_e > 1e6) ? n_eps / n_e : 0.0;
            mean_en = std::max(0.0, std::min(mean_en, 100.0));
            
            // Calculate rate for each reaction
            for (size_t r = 0; r < reactions.size(); ++r) {
                double k = reactions[r].get_rate_coeff(mean_en, config_.plasma.background_temp);
                
                // R_rate = k * prod(n_reactants^stoich)
                double R_rate = k;
                for (const auto& pair : reactions[r].reactants) {
                    int idx = pair.first;
                    int stoich = pair.second;
                    R_rate *= std::pow(densities[idx], stoich);
                }
                
                // Store in global index
                int global_idx = j * nx + i;
                rates_data[r][global_idx] = R_rate;
            }
        }
    }
    
    DMDAVecRestoreArrayDOF(grid_.get_dm(), Xloc, &x);
    DMRestoreLocalVector(grid_.get_dm(), &Xloc);
    
    // Gather rates from all processes
    // For simplicity with MPI, use MPI_Allreduce with MPI_SUM
    // (cells not owned by a rank have 0, so sum gives correct result)
    for (size_t r = 0; r < num_reactions; ++r) {
        std::vector<double> global_rates(nx * ny, 0.0);
        MPI_Allreduce(rates_data[r].data(), global_rates.data(), nx * ny, 
                      MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        rates_data[r] = global_rates;
    }
    
    output.write_rates(step, rates_data, nx, ny);
}

// --------------------------------------------------------------------------
// Physics Implementation
// --------------------------------------------------------------------------

PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void* ctx_void) {
    SolverContext* ctx = (SolverContext*)ctx_void;
    RectilinearGrid* grid = ctx->grid;
    Chemistry* chem = ctx->chemistry;
    BoundaryManager* boundary = ctx->boundary;
    
    DM dm = grid->get_dm();
    PetscErrorCode ierr;
    
    Vec Xloc;
    ierr = DMGetLocalVector(dm, &Xloc); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dm, X, INSERT_VALUES, Xloc); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dm, X, INSERT_VALUES, Xloc); CHKERRQ(ierr);
    
    // Access arrays
    PetscScalar ***x, ***x_prev, ***f;
    ierr = DMDAVecGetArrayDOF(dm, Xloc, &x); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(dm, ctx->X_prev, &x_prev); CHKERRQ(ierr); // X_prev is global, but we need local part...
    // X_prev is stored as global vector. To get ghost values, we need a local version too.
    // For time term (x - x_prev)/dt, we only need local values (no gradients of x_prev).
    // So DMDAVecGetArrayDOF on global X_prev works for local part.
    // Wait, DMDAVecGetArrayDOF on global vector gives access to local portion. 
    // Ghosts are invalid unless we scatter. But we don't need ghosts of x_prev for time term.
    
    ierr = DMDAVecGetArrayDOF(dm, F, &f); CHKERRQ(ierr);
    
    int xs, ys, xm, ym;
    ierr = DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
    int nx = grid->get_nx();
    int ny = grid->get_ny();
    
    int num_species = chem->get_num_species();
    int idx_eps = ctx->idx_n_eps;
    int idx_phi = ctx->idx_phi;
    
    double dt = ctx->dt;
    double epsilon_0 = 8.8541878e-12;
    double q_e = 1.6021766e-19; // Elementary charge magnitude
    
    // Reset residual
    for (int j=ys; j<ys+ym; ++j) {
        for (int i=xs; i<xs+xm; ++i) {
            for (int k=0; k<=idx_phi; ++k) f[j][i][k] = 0.0;
        }
    }

    // --- 1. Sources & Time Term ---
    for (int j=ys; j<ys+ym; ++j) {
        for (int i=xs; i<xs+xm; ++i) {
            double vol = grid->get_cell_volume(i, j);
            
            // Gather densities
            std::vector<double> densities(num_species);
            double density_floor = ctx->config->advanced.density_floor;
            
            for (int k=0; k<num_species; ++k) {
                // Apply floor for physics calculations (rates, transport)
                // But use actual x for time derivatives
                densities[k] = std::max(x[j][i][k], density_floor);
            }
            
            double n_e = 0.0;
            // Find electron
            const auto& species = chem->get_species();
            for(int k=0; k<num_species; ++k) {
                if(species[k].type == SpeciesType::Electron) {
                     n_e = densities[k];
                     break;
                }
            }
            double n_eps = x[j][i][idx_eps];
            // Apply energy floor
            double energy_floor = ctx->config->advanced.energy_floor;
            double mean_en = (n_e > density_floor) ? n_eps / n_e : 0.0;
            mean_en = std::max(energy_floor, std::min(mean_en, 100.0));
            
            std::vector<double> sources;
            chem->compute_source(densities, mean_en, ctx->config->plasma.background_temp, sources);
            
            // Species Residuals
            for (int k=0; k<num_species; ++k) {
                // Time term: (n - n_old)/dt * vol
                // Note: x_prev access assumes it matches local layout of global vec, which it does for owned part.
                f[j][i][k] += (x[j][i][k] - x_prev[j][i][k]) / dt * vol;
                // Source term: -S * vol
                f[j][i][k] -= sources[k] * vol;
            }
            
            // Energy Time Term
            f[j][i][idx_eps] += (x[j][i][idx_eps] - x_prev[j][i][idx_eps]) / dt * vol;
            // Energy Source: approx -Loss. (Simplified: assume -CollisionalLoss proportional to density)
            // Need detailed energy loss from Chemistry. Assuming sources has it? 
            // The current Chemistry::compute_source only returns particle sources. 
            // We need energy loss. For now, neglect or add placeholder.
            double energy_loss = 0.0; // Implement later
             f[j][i][idx_eps] -= energy_loss * vol;
        }
    }
    
    // --- 2. Fluxes (X) ---
    for (int j=ys; j<ys+ym; ++j) {
        for (int i=xs; i<xs+xm+1; ++i) { // Interfaces
            // Skip boundaries for now (handled separately or implicitly via ghost values which are BCs)
            // We iterate interfaces i=0 to nx. 
            // i corresponds to face between i-1 and i.
            // Wait, DMDA bounds xs..xs+xm are CELLS. 
            // Faces are i and i+1 for cell i.
            // We loop cells and add Flux_Right to i, subtract Flux_Right from i+1.
            // But this requires access to i+1 which might be remote.
            // Better: Compute flux at i+1/2 (Right face of cell i).
            // Add to f[j][i] (+Flux) and f[j][i+1] (-Flux). 
            // Caution with parallel: Only update OWNED cells.
            
            if (i >= nx) continue; // Boundary right handled differently?
            
            // Indices for Left and Right cells of face i+1/2
            int iL = i;
            int iR = i+1;
            
            // Get Geometry
            double dx = grid->get_dx_face(i);
            double area = grid->get_face_area_x(i, j);
            
            // Variables
            double phi_L = x[j][iL][idx_phi];
            double phi_R = x[j][iR][idx_phi];
            double dphi = phi_R - phi_L;
            double E_field = -dphi / dx;
            
            // Species Fluxes
            const auto& species = chem->get_species();
            for (int k=0; k<num_species; ++k) {
                const auto& sp = species[k];
                double n_L = x[j][iL][k];
                double n_R = x[j][iR][k];
                
                double flux = 0.0;
                
                if (sp.type == SpeciesType::Neutral) {
                     flux = compute_neutral_flux(n_L, n_R, sp.diffusion_coeff_const, 0.0, dx);
                } else {
                     // Charged
                     double mean_en = 0.0; // Should interpolate mean energy at interface
                     // Simple average
                     double ne_L = 1.0, ne_R = 1.0, eps_L = 0.0, eps_R = 0.0; 
                     // Need to find electron index again or cache it
                     // Assuming k=0 is electron for energy lookup (simplified)
                     // Proper way: pass mean_energy to get_transport
                     double mu, D;
                     sp.get_transport(2.0, mu, D); // Placeholder Energy 2.0eV
                     
                     // Charge sign for mobility direction
                     double signed_mu = mu * (sp.charge > 0 ? 1.0 : -1.0);
                     flux = compute_sg_flux(n_L, n_R, D, signed_mu, dphi, dx);
                }
                
                // Add to residuals
                // Net Flux out of cell i = Flux_R - Flux_L
                // f[i] += Flux_R * Area
                // f[i+1] -= Flux_R * Area
                
                // Note: In PETSc DMDA, we usually only update f[j][i] (owned).
                // But Flux(i+1/2) affects i (leaving) and i+1 (entering).
                // We add contribution to i. 
                // We rely on neighbor i+1 to compute Flux(i+1/2) as its Flux_L and subtract it?
                // Yes, consistent flux calculation.
                // Flux_R(i) = Flux(i+1/2). Contribution: +Flux_R(i) * Area.
                
                // Wait. Divergence term is Div(Flux). In FV: Sum(Flux_out * Area).
                // So for cell i: +Flux_Right*Area_R - Flux_Left*Area_L.
                
                f[j][iL][k] += flux * area;
                if (iR < xs+xm) { // If iR is also local (or ghost reachable but we don't write to ghost F)
                    // We shouldn't write to ghost F.
                    // Actually, just update iL. The neighbor will update itself using Flux_Left which is this Flux.
                    // BUT we must ensure Flux_Left(i+1) == Flux_Right(i).
                    // So we compute Flux at i+1/2 and add to i.
                    // And we compute Flux at i-1/2 and subtract from i.
                }
            }
        }
    }
    
    // RE-LOOP for Divergence to avoid double counting or ghost issues
    // Loop over cells i,j
    for (int j=ys; j<ys+ym; ++j) {
        for (int i=xs; i<xs+xm; ++i) {
             // 1. Right Face (i+1/2)
             // ... Compute Flux_R ...
             // f[j][i] += Flux_R * Area_R
             
             // 2. Left Face (i-1/2)
             // ... Compute Flux_L ...
             // f[j][i] -= Flux_L * Area_L
             
             // This requires computing flux twice per face or caching. 
             // Computing twice is fine for now.
             
             // POISSON: -Div(eps Grad phi) = rho
             // - [ eps*(phi_R-phi_C)/dx_R * A_R - eps*(phi_C-phi_L)/dx_L * A_L ] = rho * Vol
             // Residual: -Flux_Phi_Net - rho*Vol
             
             // Get Rho
             double rho = 0.0;
             const auto& species = chem->get_species();
             for (int k=0; k<num_species; ++k) {
                 rho += x[j][i][k] * species[k].charge * q_e;
             }
             
             // Calculate Phi Fluxes
             // Right
             double phi_C = x[j][i][idx_phi];
             double phi_R = x[j][i+1][idx_phi]; // Ghost
             double dx_R = grid->get_dx_face(i);
             double flux_phi_R = epsilon_0 * (phi_R - phi_C) / dx_R;
             
             // Left
             double phi_L = x[j][i-1][idx_phi]; // Ghost
             double dx_L = grid->get_dx_face(i-1);
             double flux_phi_L = epsilon_0 * (phi_C - phi_L) / dx_L;
             
             double div_phi = (flux_phi_R * grid->get_face_area_x(i, j) - flux_phi_L * grid->get_face_area_x(i-1, j));
             double vol = grid->get_cell_volume(i, j);
             
             f[j][i][idx_phi] -= div_phi;
             f[j][i][idx_phi] -= rho * vol;
        }
    }
    
    // --- Boundary Conditions ---
    // Apply Dirichlet for Potential at Walls
    // Loop over boundary nodes (global indices check)
    // If i == 0 (Left Wall)
    if (xs == 0) {
        for (int j=ys; j<ys+ym; ++j) {
             double t = ctx->time;
             
             // Find electrode at x_min
             std::string elec_name = "";
             bool found = false;
             for(const auto& e : ctx->config->electrodes) {
                 if (e.location == "x_min" || e.location == "left" || e.location == "Left") {
                     elec_name = e.name;
                     found = true;
                     break;
                 }
             }

             double V = 0.0;
             double gamma_see = 0.0;
             if (found) {
                 try {
                    V = boundary->get_electrode_voltage(elec_name, t);
                    gamma_see = boundary->get_electrode_gamma_see(elec_name);
                 } catch(...) {
                    // Fallback
                 }
             }
             
             // Hard Dirichlet: f = phi - V
             // Overwrite residual
             f[j][0][idx_phi] = x[j][0][idx_phi] - V;
             
             // Species BC at Wall with SEE
             const auto& species = chem->get_species();
             
             // Calculate ion flux to wall (sum of all positive ions)
             double ion_flux = 0.0;
             int electron_idx = -1;
             
             for (int k=0; k<num_species; ++k) {
                 if (species[k].type == SpeciesType::Electron) {
                     electron_idx = k;
                 } else if (species[k].charge > 0) {
                     // Ion flux approximation (drift term dominates at cathode)
                     double n_wall = x[j][0][k];
                     double E_wall = -(x[j][1][idx_phi] - x[j][0][idx_phi]) / grid->get_dx(0);
                     double mu, D;
                     species[k].get_transport(0.0, mu, D);
                     double v_ion = mu * E_wall;
                     ion_flux += n_wall * std::abs(v_ion);
                 }
             }
             
             // Apply SEE: add electron source at cathode
             if (electron_idx >= 0 && gamma_see > 0.0) {
                 double area = grid->get_face_area_x(0, j);
                 double see_source = gamma_see * ion_flux * area;
                 // Add to electron residual (negative because source term)
                 f[j][0][electron_idx] -= see_source;
             }
             
             // Set zero density BC for ions at walls (absorption)
             for (int k=0; k<num_species; ++k) {
                 if (species[k].charge > 0) {
                     f[j][0][k] = x[j][0][k]; // Enforce n = 0
                 }
             }
             
             // Apply Wall Quenching for Neutrals
             double prob = ctx->config->advanced.wall_quenching_probability;
             if (prob > 0.0) {
                 for (int k=0; k<num_species; ++k) {
                     if (species[k].type == SpeciesType::Neutral) {
                         double n_wall = x[j][0][k];
                         double m = species[k].mass;
                         double T = ctx->config->plasma.background_temp;
                         // Thermal velocity
                         double v_th = std::sqrt(8.0 * 1.380649e-23 * T / (3.14159265359 * m));
                         double gamma = 0.25 * n_wall * v_th * prob;
                         
                         double area = grid->get_face_area_x(0, j);
                         // Flux leaving domain adds to residual (loss)
                         f[j][0][k] += gamma * area;
                     }
                 }
             }
        }
    }
    // Similar for Right Wall (nx-1)
    if (xs + xm == nx) {
        for (int j=ys; j<ys+ym; ++j) {
             double t = ctx->time;
             
             std::string elec_name = "";
             bool found = false;
             for(const auto& e : ctx->config->electrodes) {
                 if (e.location == "x_max" || e.location == "right" || e.location == "Right") {
                     elec_name = e.name;
                     found = true;
                     break;
                 }
             }

             double V = 0.0;
             double gamma_see = 0.0;
             if (found) {
                 try {
                    V = boundary->get_electrode_voltage(elec_name, t);
                    gamma_see = boundary->get_electrode_gamma_see(elec_name);
                 } catch(...) {}
             }
             
             f[j][nx-1][idx_phi] = x[j][nx-1][idx_phi] - V;
             
             // Species BC with SEE
             const auto& species = chem->get_species();
             
             double ion_flux = 0.0;
             int electron_idx = -1;
             
             for (int k=0; k<num_species; ++k) {
                 if (species[k].type == SpeciesType::Electron) {
                     electron_idx = k;
                 } else if (species[k].charge > 0) {
                     double n_wall = x[j][nx-1][k];
                     double E_wall = -(x[j][nx-1][idx_phi] - x[j][nx-2][idx_phi]) / grid->get_dx(nx-1);
                     double mu, D;
                     species[k].get_transport(0.0, mu, D);
                     double v_ion = mu * E_wall;
                     ion_flux += n_wall * std::abs(v_ion);
                 }
             }
             
             if (electron_idx >= 0 && gamma_see > 0.0) {
                 double area = grid->get_face_area_x(nx-1, j);
                 double see_source = gamma_see * ion_flux * area;
                 f[j][nx-1][electron_idx] -= see_source;
             }
             
             for (int k=0; k<num_species; ++k) {
                 if (species[k].charge > 0) {
                     f[j][nx-1][k] = x[j][nx-1][k];
                 }
             }
             
             double prob = ctx->config->advanced.wall_quenching_probability;
             if (prob > 0.0) {
                 for (int k=0; k<num_species; ++k) {
                     if (species[k].type == SpeciesType::Neutral) {
                         double n_wall = x[j][nx-1][k];
                         double m = species[k].mass;
                         double T = ctx->config->plasma.background_temp;
                         double v_th = std::sqrt(8.0 * 1.380649e-23 * T / (3.14159265359 * m));
                         double gamma = 0.25 * n_wall * v_th * prob;
                         
                         double area = grid->get_face_area_x(nx-1, j);
                         f[j][nx-1][k] += gamma * area;
                     }
                 }
             }
        }
    }

    ierr = DMDAVecRestoreArrayDOF(dm, Xloc, &x); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(dm, ctx->X_prev, &x_prev); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(dm, F, &f); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &Xloc); CHKERRQ(ierr);
    
    return 0;
}

PetscErrorCode FormJacobian(SNES snes, Vec X, Mat J, Mat P, void* ctx_void) {
    SolverContext* ctx = (SolverContext*)ctx_void;
    RectilinearGrid* grid = ctx->grid;
    DM dm = grid->get_dm();
    PetscErrorCode ierr;
    
    // Identity for Species (Mass Matrix approx) + Laplacian for Poisson
    // Just putting 1s on diagonal to start, or building Laplacian for P_poisson
    
    // ... Implementation of Preconditioner ...
    // For JFNK, P just needs to be a good approx.
    // Poisson block: Laplacian
    // Transport blocks: Identity or convection-diffusion operator.
    
    MatZeroEntries(P);
    
    // Loop and set stencil
    int xs, ys, xm, ym;
    DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL);
    
    int idx_phi = ctx->idx_phi;
    
    for (int j=ys; j<ys+ym; ++j) {
        for (int i=xs; i<xs+xm; ++i) {
             MatStencil row;
             row.i = i; row.j = j; 
             
             // Diagonal Species
             for (int k=0; k<idx_phi; ++k) {
                 row.c = k;
                 MatSetValuesStencil(P, 1, &row, 1, &row, &ctx->dt, INSERT_VALUES); // 1/dt approx
             }
             
             // Poisson: Laplacian
             row.c = idx_phi;
             double v_diag = 2.0; // Simplified
             MatSetValuesStencil(P, 1, &row, 1, &row, &v_diag, INSERT_VALUES);
             // Neighbors...
        }
    }
    
    MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);
    return 0;
}

} // namespace HydroPlas
