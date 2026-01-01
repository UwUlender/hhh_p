#pragma once

#include <petscdm.h>
#include <petscdmda.h>
#include <vector>
#include <string>
#include "../config/ConfigParser.hpp"

namespace HydroPlas {

class RectilinearGrid {
public:
    explicit RectilinearGrid(const MeshConfig& config);
    ~RectilinearGrid();

    // Initialize the DMDA and metrics
    // num_dofs: number of degrees of freedom per node
    void initialize(int num_dofs = 1);

    DM get_dm() const { return dm_; }

    // Metric getters (local indices)
    double get_cell_volume(int i, int j) const;
    double get_face_area_x(int i, int j) const; 
    double get_face_area_y(int i, int j) const; 
    
    // Coordinate getters
    double get_cell_center_x(int i) const;
    double get_cell_center_y(int j) const;
    
    // Spacings
    double get_dx(int i) const; 
    double get_dy(int j) const; 
    double get_dx_face(int i) const; 
    double get_dy_face(int j) const; 

    // Global dimensions
    int get_nx() const { return nx_; }
    int get_ny() const { return ny_; }

    // Ownership
    void get_ownership_range(int& xs, int& ys, int& xm, int& ym) const;

private:
    MeshConfig config_;
    DM dm_;
    
    int nx_, ny_; 
    
    std::vector<double> x_faces_;
    std::vector<double> y_faces_;
    std::vector<double> x_centers_;
    std::vector<double> y_centers_;
    
    std::vector<double> dx_; 
    std::vector<double> dy_; 
};

} // namespace HydroPlas
