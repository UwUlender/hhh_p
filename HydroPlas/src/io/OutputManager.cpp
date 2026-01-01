#include "OutputManager.hpp"
#include <iostream>
#include <stdexcept>
#include <mpi.h>

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
    
    // Ensure /data exists
    if (!H5Lexists(file_id_, "data", H5P_DEFAULT)) {
        hid_t g = H5Gcreate(file_id_, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(g);
    }
    
    // Create time group
    std::string path = "data/step_" + std::to_string(step);
    hid_t group = H5Gcreate(file_id_, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Write Time attribute
    hid_t dataspace = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate(group, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &time);
    H5Aclose(attr);
    H5Sclose(dataspace);
    
    // Extract data from Vec X
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
    
    // Electron Energy (second to last)
    for(int j=0; j<ny; ++j) {
        for(int i=0; i<nx; ++i) {
            int idx = (j*nx + i)*dofs + (dofs-2);
            buffer[j*nx+i] = x_ptr[idx];
        }
    }
    write_dataset(group, "n_eps", buffer, dims);

    VecRestoreArrayRead(X_seq, &x_ptr);
    VecScatterDestroy(&ctx);
    VecDestroy(&X_seq);
    
    H5Gclose(group);
}

void OutputManager::write_rates(int step, const std::vector<std::vector<double>>& rates, int nx, int ny) {
     if (file_id_ < 0) return;
     std::string path = "data/step_" + std::to_string(step) + "/rates";
     
     // Ensure /rates group exists
     hid_t step_group = H5Gopen(file_id_, ("data/step_" + std::to_string(step)).c_str(), H5P_DEFAULT);
     if (step_group < 0) return;
     hid_t group = H5Gcreate(step_group, "rates", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     
     std::vector<hsize_t> dims = {(hsize_t)ny, (hsize_t)nx};
     
     for (size_t r=0; r<rates.size(); ++r) {
         write_dataset(group, "reaction_" + std::to_string(r), rates[r], dims);
     }
     
     H5Gclose(group);
     H5Gclose(step_group);
}

void OutputManager::read_state(const std::string& filename, int step, Vec X) {
    hid_t fid = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fid < 0) throw std::runtime_error("Cannot open restart file: " + filename);
    
    std::string path = "data/step_" + std::to_string(step);
    hid_t group = H5Gopen(fid, path.c_str(), H5P_DEFAULT);
    if (group < 0) {
        H5Fclose(fid);
        throw std::runtime_error("Step not found in restart file");
    }
    
    // Read grid dimensions
    int nx = grid_.get_nx();
    int ny = grid_.get_ny();
    int size = nx * ny;
    
    // Get the number of DOFs from the Vec
    PetscInt vec_size;
    VecGetSize(X, &vec_size);
    int dofs = vec_size / size;  // Total DOFs per node
    
    // Create sequential vector for reading on all ranks
    Vec X_seq;
    VecScatter ctx;
    VecScatterCreateToAll(X, &ctx, &X_seq);
    
    // Allocate buffer for reading
    std::vector<double> buffer(size);
    
    // On rank 0, read all datasets
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    PetscScalar* x_ptr;
    VecGetArray(X_seq, &x_ptr);
    
    if (rank == 0) {
        // Read species densities (n_0, n_1, ...)
        int num_species = dofs - 2;  // Subtract phi and n_eps
        for (int k = 0; k < num_species; ++k) {
            std::string dset_name = "n_" + std::to_string(k);
            hid_t dataset = H5Dopen(group, dset_name.c_str(), H5P_DEFAULT);
            if (dataset >= 0) {
                H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
                H5Dclose(dataset);
                
                // Interleave into X_seq
                for (int j = 0; j < ny; ++j) {
                    for (int i = 0; i < nx; ++i) {
                        int idx = (j * nx + i) * dofs + k;
                        x_ptr[idx] = buffer[j * nx + i];
                    }
                }
            }
        }
        
        // Read electron energy (n_eps)
        hid_t dset_eps = H5Dopen(group, "n_eps", H5P_DEFAULT);
        if (dset_eps >= 0) {
            H5Dread(dset_eps, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
            H5Dclose(dset_eps);
            
            int idx_eps = dofs - 2;
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int idx = (j * nx + i) * dofs + idx_eps;
                    x_ptr[idx] = buffer[j * nx + i];
                }
            }
        }
        
        // Read potential (phi)
        hid_t dset_phi = H5Dopen(group, "phi", H5P_DEFAULT);
        if (dset_phi >= 0) {
            H5Dread(dset_phi, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
            H5Dclose(dset_phi);
            
            int idx_phi = dofs - 1;
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int idx = (j * nx + i) * dofs + idx_phi;
                    x_ptr[idx] = buffer[j * nx + i];
                }
            }
        }
    }
    
    VecRestoreArray(X_seq, &x_ptr);
    
    // Scatter from sequential to parallel
    VecScatterBegin(ctx, X_seq, X, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(ctx, X_seq, X, INSERT_VALUES, SCATTER_REVERSE);
    
    VecScatterDestroy(&ctx);
    VecDestroy(&X_seq);
    
    H5Gclose(group);
    H5Fclose(fid);
}

void OutputManager::write_dataset(hid_t group_id, const std::string& name, const std::vector<double>& data, const std::vector<hsize_t>& dims) {
    hid_t dataspace = H5Screate_simple(dims.size(), dims.data(), NULL);
    hid_t dataset = H5Dcreate(group_id, name.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    H5Dclose(dataset);
    H5Sclose(dataspace);
}

} // namespace HydroPlas
