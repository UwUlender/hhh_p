#include "MeshGenerator.hpp"
#include <iostream>

namespace HydroPlas {

MeshGenerator::MeshGenerator(const DomainConfig& config) : config_(config) {}

MeshGenerator::~MeshGenerator() {}

PetscErrorCode MeshGenerator::create_dm(DM* dm, int num_dofs) {
    PetscErrorCode ierr;

    // Always use 2D DMDA. For 1D, Ny=1.
    PetscInt Ny = (config_.Ny > 0) ? config_.Ny : 1;
    PetscReal Ly = (config_.Ly > 0.0) ? config_.Ly : 1.0; // Dummy Ly if 1D

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
                        config_.Nx, Ny, PETSC_DECIDE, PETSC_DECIDE,
                        num_dofs, 1, NULL, NULL, dm); CHKERRQ(ierr);


    ierr = DMSetFromOptions(*dm); CHKERRQ(ierr);
    ierr = DMSetUp(*dm); CHKERRQ(ierr);

    // Set coordinates
    ierr = DMDASetUniformCoordinates(*dm, 0.0, config_.Lx, 0.0, config_.Ly, 0.0, 0.0); CHKERRQ(ierr);

    return 0;
}

} // namespace HydroPlas
