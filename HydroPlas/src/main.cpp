#include <petsc.h>
#include <iostream>
#include "config/ConfigParser.hpp"
#include "mesh/RectilinearGrid.hpp"
#include "chemistry/Chemistry.hpp"
#include "solver/PlasmaSolver.hpp"
#include "io/OutputManager.hpp"

static char help[] = "Hydrodynamic Plasma Simulation Code\n\n";

using namespace HydroPlas;

int main(int argc, char **argv) {
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv, (char*)0, help); if (ierr) return ierr;

    PetscMPIInt rank, size;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    PetscPrintf(PETSC_COMM_WORLD, "Running on %d MPI processes\n", size);

    std::string config_file = "config/default_config.yaml";
    char conf[256];
    PetscBool flg;
    PetscOptionsGetString(NULL, NULL, "-config", conf, sizeof(conf), &flg);
    if (!flg) PetscOptionsGetString(NULL, NULL, "--config", conf, sizeof(conf), &flg);
    if (flg) config_file = conf;

    try {
        PetscPrintf(PETSC_COMM_WORLD, "Reading configuration from %s\n", config_file.c_str());
        // Note: For actual run, ensure default_config.yaml exists or pass valid one
        ConfigParser parser(config_file);
        SimulationConfig config = parser.get_config();
        
        PetscPrintf(PETSC_COMM_WORLD, "Initializing Grid...\n");
        RectilinearGrid grid(config.mesh);
        
        PetscPrintf(PETSC_COMM_WORLD, "Initializing Chemistry...\n");
        Chemistry chemistry(config);

        PetscPrintf(PETSC_COMM_WORLD, "Initializing Boundary...\n");
        BoundaryConfig b_conf;
        b_conf.electrodes = config.electrodes;
        BoundaryManager boundary(b_conf);
        
        PetscPrintf(PETSC_COMM_WORLD, "Initializing Solver...\n");
        PlasmaSolver solver(grid, chemistry, boundary, config);
        solver.initialize(); 
        
        PetscPrintf(PETSC_COMM_WORLD, "Initializing Output...\n");
        OutputManager output(config.output, grid, chemistry);
        output.write_mesh();
        
        double t = 0.0;
        int step = 0;
        
        // Restart Logic (command-line overrides config)
        char restart_file[256] = "";
        PetscBool has_restart;
        PetscOptionsGetString(NULL, NULL, "-restart", restart_file, sizeof(restart_file), &has_restart);
        if (!has_restart) PetscOptionsGetString(NULL, NULL, "--restart", restart_file, sizeof(restart_file), &has_restart);
        
        bool should_restart = has_restart || config.restart.enabled;
        std::string restart_filename = has_restart ? std::string(restart_file) : config.restart.file;
        int restart_step = config.restart.step;
        
        if (should_restart) {
             PetscPrintf(PETSC_COMM_WORLD, "Restarting from %s...\n", restart_filename.c_str());
             PetscOptionsGetInt(NULL, NULL, "-restart_step", &restart_step, NULL);
             PetscOptionsGetInt(NULL, NULL, "--restart_step", &restart_step, NULL);
             output.read_state(restart_filename, restart_step, solver.get_solution());
             step = restart_step;
             // t should be read from file too. 
             // For now assume t is updated or manual.
        }
        
        // Use solver parameters from config
        double dt = config.solver.time_step;
        double t_end = config.solver.end_time;
        
        // Allow command-line override of time parameters
        PetscReal dt_override, tend_override;
        PetscBool dt_set, tend_set;
        PetscOptionsGetReal(NULL, NULL, "-dt", &dt_override, &dt_set);
        if (!dt_set) PetscOptionsGetReal(NULL, NULL, "--dt", &dt_override, &dt_set);
        PetscOptionsGetReal(NULL, NULL, "-tend", &tend_override, &tend_set);
        if (!tend_set) PetscOptionsGetReal(NULL, NULL, "--tend", &tend_override, &tend_set);
        if (dt_set) dt = dt_override;
        if (tend_set) t_end = tend_override;
        
        int output_freq = (config.output.frequency_step > 0) ? config.output.frequency_step : 100;
        
        PetscPrintf(PETSC_COMM_WORLD, "Solver Configuration:\n");
        PetscPrintf(PETSC_COMM_WORLD, "  Type: %s\n", config.solver.type.c_str());
        PetscPrintf(PETSC_COMM_WORLD, "  Time step: %g s\n", dt);
        PetscPrintf(PETSC_COMM_WORLD, "  End time: %g s\n", t_end);
        PetscPrintf(PETSC_COMM_WORLD, "  Tolerance: %g\n", config.solver.tolerance);
        PetscPrintf(PETSC_COMM_WORLD, "  Max iterations: %d\n", config.solver.max_iterations);
        
        PetscPrintf(PETSC_COMM_WORLD, "Starting Simulation...\n");
        
        while (t < t_end) {
            solver.solve_step(dt, t);
            
            t += dt;
            step++;
            
            if (step % output_freq == 0) {
                 PetscPrintf(PETSC_COMM_WORLD, "Step %d, Time %g\n", step, t);
                 output.write_state(t, step, solver.get_solution(), chemistry.get_num_species());
                 if (config.output.save_rates) {
                     solver.save_rates(output, step);
                 }
            }
        }
        
    } catch (const std::exception& e) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: %s\n", e.what());
    }

    ierr = PetscFinalize();
    return ierr;
}
