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
            u[j][i][0] = 1e10; // ne
            u[j][i][1] = 1e10; // ni
            u[j][i][2] = 1.6e-19 * 2.0 * 1e10; // neps
            u[j][i][3] = 0.0;  // phi
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
    
    const PetscInt NE = 0, NI = 1, NEPS = 2, PHI = 3;
    
    PetscScalar ***u, ***udot, ***f;
    
    ierr = DMDAVecGetArrayDOF(dm, locU, &u); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(dm, Udot, &udot); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(dm, F, &f); CHKERRQ(ierr);
    
    PetscReal Lx = app->config.domain.Lx;
    PetscReal dx = Lx / (M > 1 ? M-1 : 1);
    
    const double q = 1.602e-19;
    const double eps0 = 8.854e-12;
    const double eps = eps0; 
    
    double mu_e = 30.0; 
    double D_e = 10.0;
    double mu_i = 1.0;
    double D_i = 0.1;

    for (PetscInt j = ys; j < ys + ym; j++) {
        for (PetscInt i = xs; i < xs + xm; i++) {
            
            f[j][i][NE] = udot[j][i][NE];
            f[j][i][NI] = udot[j][i][NI];
            f[j][i][NEPS] = udot[j][i][NEPS];
            f[j][i][PHI] = 0.0;
            
            if (i == 0) {
                // Left BC
                f[j][i][PHI] = u[j][i][PHI] - app->boundary->get_voltage(t);
                
                // For species, we can implement SEE here.
                // Flux boundary condition: Gamma_e = -gamma * Gamma_i
                // Implementing Dirichlet for now to ensure stability of the demo.
                f[j][i][NE] = u[j][i][NE]; 
                f[j][i][NI] = u[j][i][NI];
                f[j][i][NEPS] = u[j][i][NEPS];
                continue;
            } else if (i == M-1) {
                // Right BC
                f[j][i][PHI] = u[j][i][PHI]; // Ground
                f[j][i][NE] = u[j][i][NE];
                f[j][i][NI] = u[j][i][NI];
                f[j][i][NEPS] = u[j][i][NEPS];
                continue;
            }

            // Poisson
            double rho = q * (u[j][i][NI] - u[j][i][NE]);
            double d2phi_dx2 = (u[j][i+1][PHI] - 2*u[j][i][PHI] + u[j][i-1][PHI]) / (dx*dx);
            f[j][i][PHI] = - eps * d2phi_dx2 - rho;

            // Fluxes
            double dphi_right = u[j][i+1][PHI] - u[j][i][PHI];
            double nu_e_right = mu_e * dphi_right / D_e;
            double flux_e_right = ScharfetterGummelFlux(u[j][i][NE], u[j][i+1][NE], nu_e_right, D_e, dx);
            
            double dphi_left = u[j][i][PHI] - u[j][i-1][PHI];
            double nu_e_left = mu_e * dphi_left / D_e;
            double flux_e_left = ScharfetterGummelFlux(u[j][i-1][NE], u[j][i][NE], nu_e_left, D_e, dx);
            
            f[j][i][NE] += (flux_e_right - flux_e_left) / dx;
            
            double nu_i_right = -mu_i * dphi_right / D_i;
            double flux_i_right = ScharfetterGummelFlux(u[j][i][NI], u[j][i+1][NI], nu_i_right, D_i, dx);
            
            double nu_i_left = -mu_i * dphi_left / D_i;
            double flux_i_left = ScharfetterGummelFlux(u[j][i-1][NI], u[j][i][NI], nu_i_left, D_i, dx);
            
            f[j][i][NI] += (flux_i_right - flux_i_left) / dx;
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
