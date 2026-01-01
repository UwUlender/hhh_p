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

    std::string config_file = "config/default_config.yaml";
    char conf[256];
    PetscBool flg;
    PetscOptionsGetString(NULL, NULL, "-config", conf, sizeof(conf), &flg);
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
        OutputManager output(config.output, grid);
        output.write_mesh();
        
        double t = 0.0;
        int step = 0;
        
        // Restart Logic
        char restart_file[256] = "";
        PetscBool has_restart;
        PetscOptionsGetString(NULL, NULL, "-restart", restart_file, sizeof(restart_file), &has_restart);
        if (has_restart) {
             PetscPrintf(PETSC_COMM_WORLD, "Restarting from %s...\n", restart_file);
             int r_step = 0;
             PetscOptionsGetInt(NULL, NULL, "-restart_step", &r_step, NULL);
             output.read_state(restart_file, r_step, solver.get_solution());
             step = r_step;
             // t should be read from file too. 
             // For now assume t is updated or manual.
        }
        
        double dt = 1e-12; 
        double t_end = 1e-9; 
        
        int output_freq = (config.output.frequency_step > 0) ? config.output.frequency_step : 100;
        
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
