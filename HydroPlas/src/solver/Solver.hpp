#pragma once
#include <petsc.h>
#include <petscdmda.h>
#include <petscts.h>
#include "../config/ConfigParser.hpp"
#include "../chemistry/LookupTable.hpp"
#include "../boundary/BoundaryManager.hpp"

namespace HydroPlas {

struct AppCtx {
    DM dm;
    SimulationConfig config;
    LookupTable lookup;
    BoundaryManager* boundary;
};

class Solver {
public:
    Solver(DM dm, const SimulationConfig& config);
    ~Solver();

    PetscErrorCode init();
    PetscErrorCode solve();

    // Accessor for AppCtx if needed by callbacks (they get it via void*)
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
