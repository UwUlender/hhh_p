#pragma once
#include <petscsnes.h>
#include <petscdmda.h>
#include "../mesh/RectilinearGrid.hpp"
#include "../chemistry/Chemistry.hpp"
#include "../config/ConfigParser.hpp"

namespace HydroPlas {

struct SolverContext {
    RectilinearGrid* grid;
    Chemistry* chemistry;
    SimulationConfig* config;
    
    int idx_phi;
    int idx_n_eps;
};

class PlasmaSolver {
public:
    PlasmaSolver(RectilinearGrid& grid, Chemistry& chemistry, SimulationConfig& config);
    ~PlasmaSolver();

    void initialize();
    
    // Solve for one time step
    void solve_step(double dt, double time);
    
    void save_state(const std::string& filename, int step, double time);

    SolverContext* get_context() { return &ctx_; }
    
    // Accessor for Output
    Vec get_solution() const { return X_; }

private:
    RectilinearGrid& grid_;
    Chemistry& chemistry_;
    SimulationConfig& config_;
    
    SNES snes_;
    Vec X_; 
    Vec F_; 
    
    SolverContext ctx_;
    
    double current_dt_;
    double current_time_;
    
    PC pc_;
    Mat P_poisson_; 
    
    void setup_dofs();
    void setup_solver();
};

PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void* ctx);
PetscErrorCode FormJacobian(SNES snes, Vec X, Mat J, Mat P, void* ctx);

} // namespace HydroPlas
