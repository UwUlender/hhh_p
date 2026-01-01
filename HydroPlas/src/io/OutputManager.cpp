#include "OutputManager.hpp"
#include <iostream>

namespace HydroPlas {

OutputManager::OutputManager(const OutputConfig& config, const RectilinearGrid& grid) 
    : config_(config), grid_(grid) {
    
    // Create file
    file_id_ = H5Fcreate(config.filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id_ < 0) {
        std::cerr << "Failed to create HDF5 file" << std::endl;
    }
}

OutputManager::~OutputManager() {
    if (file_id_ >= 0) H5Fclose(file_id_);
}

void OutputManager::write_mesh() {
    if (file_id_ < 0) return;
    
    hid_t group = H5Gcreate(file_id_, "mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Collect coords
    int nx = grid_.get_nx();
    int ny = grid_.get_ny();
    std::vector<double> x(nx), y(ny);
    for(int i=0; i<nx; ++i) x[i] = grid_.get_cell_center_x(i);
    for(int j=0; j<ny; ++j) y[j] = grid_.get_cell_center_y(j);
    
    write_dataset(group, "x_coords", x, {(hsize_t)nx});
    write_dataset(group, "y_coords", y, {(hsize_t)ny});
    
    H5Gclose(group);
}

void OutputManager::write_state(double time, int step, Vec X, int num_species) {
    if (file_id_ < 0) return;
    
    std::string group_name = "data/time_" + std::to_string(time);
    // Create group path hierarchy if needed, for now flat assumption or check
    // HDF5 requires groups to exist.
    // Use H5Gcreate with intermediate? H5Lexists...
    // Simplification: just /data group then subgroup
    
    // Ensure /data exists
    if (!H5Lexists(file_id_, "data", H5P_DEFAULT)) {
        hid_t g = H5Gcreate(file_id_, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(g);
    }
    
    // Create time group
    std::string path = "data/time_" + std::to_string(step); // Use step for uniqueness
    hid_t group = H5Gcreate(file_id_, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Extract data from Vec X
    // Need scatter to zero to get all data on rank 0 for writing (simple parallel I/O)
    // Or simpler: parallel HDF5.
    // For this snippet, assume small scale, gather to rank 0.
    
    VecScatter ctx;
    Vec X_seq;
    VecScatterCreateToAll(X, &ctx, &X_seq);
    VecScatterBegin(ctx, X, X_seq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, X, X_seq, INSERT_VALUES, SCATTER_FORWARD);
    
    // Get array
    const PetscScalar* x_ptr;
    VecGetArrayRead(X_seq, &x_ptr);
    
    // De-interleave
    int nx = grid_.get_nx();
    int ny = grid_.get_ny();
    int dofs = num_species + 2;
    int size = nx * ny;
    
    std::vector<double> buffer(size);
    std::vector<hsize_t> dims = {(hsize_t)ny, (hsize_t)nx}; // 2D layout
    
    // Write each species
    for(int k=0; k<num_species; ++k) {
        for(int j=0; j<ny; ++j) {
            for(int i=0; i<nx; ++i) {
                int idx = (j*nx + i)*dofs + k;
                buffer[j*nx+i] = x_ptr[idx];
            }
        }
        write_dataset(group, "n_" + std::to_string(k), buffer, dims);
    }
    
    // Potential (last)
    for(int j=0; j<ny; ++j) {
        for(int i=0; i<nx; ++i) {
            int idx = (j*nx + i)*dofs + (dofs-1);
            buffer[j*nx+i] = x_ptr[idx];
        }
    }
    write_dataset(group, "phi", buffer, dims);

    VecRestoreArrayRead(X_seq, &x_ptr);
    VecScatterDestroy(&ctx);
    VecDestroy(&X_seq);
    
    H5Gclose(group);
}

void OutputManager::write_dataset(hid_t group_id, const std::string& name, const std::vector<double>& data, const std::vector<hsize_t>& dims) {
    hid_t dataspace = H5Screate_simple(dims.size(), dims.data(), NULL);
    hid_t dataset = H5Dcreate(group_id, name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    H5Dclose(dataset);
    H5Sclose(dataspace);
}

} // namespace HydroPlas
