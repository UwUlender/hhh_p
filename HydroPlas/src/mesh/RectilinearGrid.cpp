#include "RectilinearGrid.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>

namespace HydroPlas {

RectilinearGrid::RectilinearGrid(const MeshConfig& config) : config_(config), dm_(nullptr) {}

RectilinearGrid::~RectilinearGrid() {
    if (dm_) DMDestroy(&dm_);
}

void RectilinearGrid::initialize(int num_dofs) {
    if (dm_) DMDestroy(&dm_);

    // Process X coordinates
    if (config_.x_nodes.size() > 1) {
        x_faces_ = config_.x_nodes;
    } else {
        int N = 50;
        double L = 0.01;
        x_faces_.resize(N + 1);
        for(int i=0; i<=N; ++i) x_faces_[i] = (double)i / N * L;
    }

    // Process Y coordinates
    if (config_.y_nodes.size() > 1) {
        y_faces_ = config_.y_nodes;
    } else {
        if (config_.type.find("2D") != std::string::npos) {
             int N = 50; 
             double L = 0.01;
             y_faces_.resize(N + 1);
             for(int i=0; i<=N; ++i) y_faces_[i] = (double)i / N * L;
        } else {
             y_faces_ = {0.0, 1.0}; 
        }
    }
    
    nx_ = x_faces_.size() - 1;
    ny_ = y_faces_.size() - 1;

    x_centers_.resize(nx_);
    dx_.resize(nx_);
    for(int i=0; i<nx_; ++i) {
        dx_[i] = x_faces_[i+1] - x_faces_[i];
        x_centers_[i] = 0.5 * (x_faces_[i] + x_faces_[i+1]);
    }

    y_centers_.resize(ny_);
    dy_.resize(ny_);
    for(int j=0; j<ny_; ++j) {
        dy_[j] = y_faces_[j+1] - y_faces_[j];
        y_centers_[j] = 0.5 * (y_faces_[j] + y_faces_[j+1]);
    }

    DMBoundaryType bx = DM_BOUNDARY_NONE;
    DMBoundaryType by = DM_BOUNDARY_NONE;
    DMDAStencilType stype = DMDA_STENCIL_STAR; 
    
    // Create DMDA with requested DOFs
    PetscErrorCode ierr;
    ierr = DMDACreate2d(PETSC_COMM_WORLD, bx, by, stype, 
                        nx_, ny_, PETSC_DECIDE, PETSC_DECIDE, 
                        num_dofs, 1, // stencil width 1
                        nullptr, nullptr, &dm_);
    
    DMSetFromOptions(dm_);
    DMSetUp(dm_);
    
    DMDASetUniformCoordinates(dm_, x_faces_.front(), x_faces_.back(), y_faces_.front(), y_faces_.back(), 0.0, 1.0);
}

double RectilinearGrid::get_cell_volume(int i, int j) const {
    return dx_[i] * dy_[j];
}

double RectilinearGrid::get_face_area_x(int i, int j) const {
    return dy_[j];
}

double RectilinearGrid::get_face_area_y(int i, int j) const {
    return dx_[i];
}

double RectilinearGrid::get_cell_center_x(int i) const {
    return x_centers_[i];
}

double RectilinearGrid::get_cell_center_y(int j) const {
    return y_centers_[j];
}

double RectilinearGrid::get_dx(int i) const { return dx_[i]; }
double RectilinearGrid::get_dy(int j) const { return dy_[j]; }

double RectilinearGrid::get_dx_face(int i) const {
    if (i >= nx_-1) return dx_[nx_-1]; 
    return x_centers_[i+1] - x_centers_[i];
}

double RectilinearGrid::get_dy_face(int j) const {
    if (j >= ny_-1) return dy_[ny_-1];
    return y_centers_[j+1] - y_centers_[j];
}

void RectilinearGrid::get_ownership_range(int& xs, int& ys, int& xm, int& ym) const {
    DMDAGetCorners(dm_, &xs, &ys, nullptr, &xm, &ym, nullptr);
}

} // namespace HydroPlas
