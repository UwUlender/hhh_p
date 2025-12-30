#include <petscsys.h>
#include <iostream>
#include "config/ConfigParser.hpp"
#include "mesh/MeshGenerator.hpp"
#include "solver/Solver.hpp"
#include "chemistry/BolsigInterface.hpp"

static char help[] = "HydroPlas: Hydrodynamic Plasma Simulation Framework.\n\n";

int main(int argc, char **argv) {
    PetscErrorCode ierr;
    
    // Initialize PETSc
    ierr = PetscInitialize(&argc, &argv, (char *)0, help); if (ierr) return ierr;

    try {
        std::cout << "Initializing HydroPlas..." << std::endl;

        // Load Configuration
        std::string config_file = "config/default_config.json";
        HydroPlas::ConfigParser parser(config_file);
        HydroPlas::SimulationConfig config = parser.get_config();

        std::cout << "Configuration loaded." << std::endl;
        std::cout << "Domain: " << config.domain.Nx << " x " << config.domain.Ny << std::endl;

        if (config.chemistry.mode == "Inline BOLSIG+") {
            HydroPlas::BolsigInterface::run_bolsig(config.chemistry, "bolsig_output.dat");
        }

        // Create Mesh
        HydroPlas::MeshGenerator meshGen(config.domain);
        DM dm;
        ierr = meshGen.create_dm(&dm); CHKERRQ(ierr);

        std::cout << "Mesh generated." << std::endl;

        // Create and Init Solver
        HydroPlas::Solver solver(dm, config);
        ierr = solver.init(); CHKERRQ(ierr);
        
        std::cout << "Solver initialized. Starting simulation..." << std::endl;
        
        // Solve
        ierr = solver.solve(); CHKERRQ(ierr);

        std::cout << "Simulation completed." << std::endl;

        // Cleanup
        ierr = DMDestroy(&dm); CHKERRQ(ierr);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        // Proceed to finalize
    }

    std::cout << "Finalizing..." << std::endl;
    ierr = PetscFinalize();
    return ierr;
}
