#pragma once
#include <petsc.h>
#include <petscdmda.h>
#include <petscts.h>
#include "../config/ConfigParser.hpp"
#include "../chemistry/LookupTable.hpp"
#include "../boundary/BoundaryManager.hpp"
#include <vector>
#include <map>

namespace HydroPlas {

struct AppCtx {
    DM dm;
    SimulationConfig config;
    LookupTable lookup;
    BoundaryManager* boundary;
    
    // Mappings
    std::map<std::string, int> species_map; // Name -> DOF Index
    int num_excited;
    int idx_ne, idx_ni, idx_neps, idx_phi, idx_sigma;
    int idx_excited_start;
};

class Solver {
public:
    Solver(DM dm, const SimulationConfig& config);
    ~Solver();

    PetscErrorCode init();
    PetscErrorCode solve();

    // Accessor for AppCtx if needed by callbacks
    AppCtx& get_ctx() { return ctx_; }

private:
    DM dm_;
    SimulationConfig config_;
    TS ts_;
    AppCtx ctx_;
};

// PETSc callbacks
PetscErrorCode FormIFunction(TS ts, PetscReal t, Vec U, Vec Udot, Vec F, void* ctx);
PetscErrorCode FormIJacobian(TS ts, PetscReal t, Vec U, Vec Udot, PetscReal shift, Mat J, Mat P, void* ctx);

} // namespace HydroPlas
