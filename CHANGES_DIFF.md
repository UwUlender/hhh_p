# Exact Code Changes for Segmentation Fault Fix

## File: src/solver/PlasmaSolver.cpp

### Change 1: setup_dofs() function - Complete Replacement

**Location:** Lines 43-62 (old) → Lines 43-117 (new)

**OLD CODE:**
```cpp
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
```

**NEW CODE:**
```cpp
void PlasmaSolver::setup_dofs() {
    int num_species = chemistry_.get_num_species();
    int dofs = num_species + 2; // Species + Energy + Potential
    grid_.initialize(dofs);
    
    DMCreateGlobalVector(grid_.get_dm(), &X_);
    VecDuplicate(X_, &F_);
    VecDuplicate(X_, &ctx_.X_prev);
    
    // Initialize X with default values
    VecSet(X_, 0.0);
    
    // Apply initial conditions from config
    DM dm = grid_.get_dm();
    Vec Xloc;
    DMGetLocalVector(dm, &Xloc);
    DMGlobalToLocalBegin(dm, X_, INSERT_VALUES, Xloc);
    DMGlobalToLocalEnd(dm, X_, INSERT_VALUES, Xloc);
    
    PetscScalar ***x;
    DMDAVecGetArrayDOF(dm, Xloc, &x);
    
    int xs, ys, xm, ym;
    DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym);
    
    // Set default background values
    for (int j = ys; j < ys + ym; ++j) {
        for (int i = xs; i < xs + xm; ++i) {
            // Default values
            for (int k = 0; k < num_species; ++k) {
                x[j][i][k] = config_.advanced.density_floor; // Use density floor as default
            }
            x[j][i][ctx_.idx_n_eps] = 1.5 * config_.advanced.density_floor; // Default energy density
            x[j][i][ctx_.idx_phi] = 0.0; // Potential
        }
    }
    
    // Apply initial conditions from config
    const auto& species = chemistry_.get_species();
    for (const auto& ic : config_.initial_conditions) {
        if (ic.name == "e_energy") {
            // Find electron species to scale energy density
            int e_idx = -1;
            for (int k = 0; k < num_species; ++k) {
                if (species[k].type == SpeciesType::Electron) {
                    e_idx = k;
                    break;
                }
            }
            
            if (e_idx >= 0) {
                for (int j = ys; j < ys + ym; ++j) {
                    for (int i = xs; i < xs + xm; ++i) {
                        double n_e = x[j][i][e_idx];
                        x[j][i][ctx_.idx_n_eps] = ic.value * n_e; // n_eps = mean_energy * n_e
                    }
                }
            }
        } else {
            // Species density
            int species_idx = chemistry_.get_species_index(ic.name);
            if (species_idx >= 0) {
                for (int j = ys; j < ys + ym; ++j) {
                    for (int i = xs; i < xs + xm; ++i) {
                        x[j][i][species_idx] = ic.value;
                    }
                }
            }
        }
    }
    
    DMDAVecRestoreArrayDOF(dm, Xloc, &x);
    DMLocalToGlobalBegin(dm, Xloc, INSERT_VALUES, X_);
    DMLocalToGlobalEnd(dm, Xloc, INSERT_VALUES, X_);
    DMRestoreLocalVector(dm, &Xloc);
}
```

---

### Change 2: FormFunction() - X_prev Array Access

**Location:** Lines 234-243 (old) → Lines 289-300 (new)

**OLD CODE:**
```cpp
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
```

**NEW CODE:**
```cpp
    // Access arrays
    PetscScalar ***x, ***x_prev, ***f;
    ierr = DMDAVecGetArrayDOF(dm, Xloc, &x); CHKERRQ(ierr);
    
    // Need local version of X_prev for accessing ghosted values
    Vec Xprev_loc;
    ierr = DMGetLocalVector(dm, &Xprev_loc); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(dm, ctx->X_prev, INSERT_VALUES, Xprev_loc); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dm, ctx->X_prev, INSERT_VALUES, Xprev_loc); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(dm, Xprev_loc, &x_prev); CHKERRQ(ierr);
    
    ierr = DMDAVecGetArrayDOF(dm, F, &f); CHKERRQ(ierr);
```

---

### Change 3: FormFunction() - Poisson Equation Boundary Check

**Location:** Lines 391-435 (old) → Lines 447-482 (new)

**OLD CODE:**
```cpp
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
```

**NEW CODE:**
```cpp
    // --- 3. Poisson Equation: -Div(eps Grad phi) = rho ---
    // Loop over interior cells only, boundaries handled separately
    for (int j=ys; j<ys+ym; ++j) {
        for (int i=xs; i<xs+xm; ++i) {
             // Skip boundary cells - they get Dirichlet BC
             if (i == 0 || i == nx-1) continue;
             
             // POISSON: -Div(eps Grad phi) = rho
             // - [ eps*(phi_R-phi_C)/dx_R * A_R - eps*(phi_C-phi_L)/dx_L * A_L ] = rho * Vol
             // Residual: -Flux_Phi_Net - rho*Vol
             
             // Get Rho (charge density)
             double rho = 0.0;
             const auto& species = chem->get_species();
             for (int k=0; k<num_species; ++k) {
                 rho += x[j][i][k] * species[k].charge * q_e;
             }
             
             // Calculate Phi Fluxes
             double phi_C = x[j][i][idx_phi];
             
             // Right face
             double phi_R = x[j][i+1][idx_phi];
             double dx_R = grid->get_dx_face(i);
             double flux_phi_R = epsilon_0 * (phi_R - phi_C) / dx_R;
             
             // Left face
             double phi_L = x[j][i-1][idx_phi];
             double dx_L = grid->get_dx_face(i-1);
             double flux_phi_L = epsilon_0 * (phi_C - phi_L) / dx_L;
             
             double div_phi = (flux_phi_R * grid->get_face_area_x(i, j) - flux_phi_L * grid->get_face_area_x(i-1, j));
             double vol = grid->get_cell_volume(i, j);
             
             f[j][i][idx_phi] -= div_phi;
             f[j][i][idx_phi] -= rho * vol;
        }
    }
```

---

### Change 4: FormFunction() - Cleanup Section

**Location:** Lines 602-605 (old) → Lines 652-656 (new)

**OLD CODE:**
```cpp
    ierr = DMDAVecRestoreArrayDOF(dm, Xloc, &x); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(dm, ctx->X_prev, &x_prev); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(dm, F, &f); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &Xloc); CHKERRQ(ierr);
```

**NEW CODE:**
```cpp
    ierr = DMDAVecRestoreArrayDOF(dm, Xloc, &x); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(dm, Xprev_loc, &x_prev); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(dm, F, &f); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &Xloc); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &Xprev_loc); CHKERRQ(ierr);
```

---

## Summary of Changes

1. **setup_dofs()**: Complete rewrite to apply initial conditions from config
2. **FormFunction()**: Added proper local vector handling for X_prev
3. **FormFunction()**: Added boundary skip condition in Poisson loop
4. **FormFunction()**: Added cleanup for Xprev_loc

These changes fix the segmentation fault and ensure correct simulation initialization.
