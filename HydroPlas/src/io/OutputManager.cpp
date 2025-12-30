#include "OutputManager.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>

namespace HydroPlas {

OutputManager::OutputManager(const SimulationConfig& config, DM dm, int num_excited)
    : config_(config), dm_(dm), num_excited_(num_excited) {
    
    // Determine output directory
    output_dir_ = "output";
    mkdir(output_dir_.c_str(), 0755);
    
    // Determine number of DOFs
    DMDAGetInfo(dm, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &num_dofs_, 
                NULL, NULL, NULL, NULL, NULL);
    
    // Check if HDF5 is available
#ifdef USE_HDF5
    use_hdf5_ = true;
    std::cout << "OutputManager: HDF5 output enabled (OpenPMD-compatible)" << std::endl;
#else
    use_hdf5_ = false;
    std::cout << "OutputManager: Using text output (HDF5 not available)" << std::endl;
#endif
}

OutputManager::~OutputManager() {
}

void OutputManager::set_field_names(const std::vector<std::string>& names) {
    field_names_ = names;
}

PetscErrorCode OutputManager::write_output(Vec U, PetscReal time, PetscInt step) {
    PetscErrorCode ierr;
    
    // Always write text output for debugging
    ierr = write_text_output(U, time, step); CHKERRQ(ierr);
    
#ifdef USE_HDF5
    if (use_hdf5_) {
        ierr = write_hdf5_output(U, time, step); CHKERRQ(ierr);
    }
#endif
    
    return 0;
}

PetscErrorCode OutputManager::write_text_output(Vec U, PetscReal time, PetscInt step) {
    PetscErrorCode ierr;
    
    // Get array
    PetscScalar ***u;
    ierr = DMDAVecGetArrayDOF(dm_, U, &u); CHKERRQ(ierr);
    
    PetscInt xs, ys, xm, ym, M, N;
    ierr = DMDAGetCorners(dm_, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
    ierr = DMDAGetInfo(dm_, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, 
                       NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
    
    // Write each field to separate file
    for (int dof = 0; dof < num_dofs_; ++dof) {
        std::ostringstream filename;
        filename << output_dir_ << "/";
        if (dof < field_names_.size()) {
            filename << field_names_[dof];
        } else {
            filename << "field_" << dof;
        }
        filename << "_" << std::setfill('0') << std::setw(6) << step << ".txt";
        
        std::ofstream file(filename.str());
        file << std::scientific << std::setprecision(10);
        file << "# Time: " << time << std::endl;
        file << "# Field: " << (dof < field_names_.size() ? field_names_[dof] : "unknown") << std::endl;
        file << "# Grid: " << M << " x " << N << std::endl;
        
        for (PetscInt j = ys; j < ys + ym; j++) {
            for (PetscInt i = xs; i < xs + xm; i++) {
                file << i << " " << j << " " << u[j][i][dof] << std::endl;
            }
        }
        file.close();
    }
    
    ierr = DMDAVecRestoreArrayDOF(dm_, U, &u); CHKERRQ(ierr);
    return 0;
}

#ifdef USE_HDF5

hid_t OutputManager::create_openpmd_file(const std::string& filename) {
    // Create HDF5 file with OpenPMD structure
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    // Add OpenPMD attributes
    hid_t root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
    
    // OpenPMD version
    const char* openpmd_version = "1.1.0";
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(openpmd_version));
    hid_t attr = H5Acreate(root_id, "openPMD", atype, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, atype, openpmd_version);
    H5Aclose(attr);
    H5Tclose(atype);
    
    // OpenPMD extension
    const char* extension = "HydroPlas";
    atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(extension));
    attr = H5Acreate(root_id, "openPMDextension", atype, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, atype, extension);
    H5Aclose(attr);
    H5Tclose(atype);
    
    // Base path
    const char* basePath = "/data/%T/";
    atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(basePath));
    attr = H5Acreate(root_id, "basePath", atype, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, atype, basePath);
    H5Aclose(attr);
    H5Tclose(atype);
    
    // Meshes path
    const char* meshesPath = "meshes/";
    atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(meshesPath));
    attr = H5Acreate(root_id, "meshesPath", atype, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, atype, meshesPath);
    H5Aclose(attr);
    H5Tclose(atype);
    
    H5Gclose(root_id);
    return file_id;
}

void OutputManager::write_openpmd_metadata(hid_t file_id, PetscReal time, PetscInt step) {
    // Create iteration group: /data/{step}
    std::ostringstream group_path;
    group_path << "/data/" << step;
    hid_t iter_group = H5Gcreate(file_id, group_path.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Write time attribute
    hid_t aspace = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate(iter_group, "time", H5T_NATIVE_DOUBLE, aspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &time);
    H5Aclose(attr);
    H5Sclose(aspace);
    
    // Write dt attribute (time step)
    double dt = config_.time.dt;
    aspace = H5Screate(H5S_SCALAR);
    attr = H5Acreate(iter_group, "dt", H5T_NATIVE_DOUBLE, aspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &dt);
    H5Aclose(attr);
    H5Sclose(aspace);
    
    // Write timeUnitSI
    double timeUnitSI = 1.0; // seconds
    aspace = H5Screate(H5S_SCALAR);
    attr = H5Acreate(iter_group, "timeUnitSI", H5T_NATIVE_DOUBLE, aspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &timeUnitSI);
    H5Aclose(attr);
    H5Sclose(aspace);
    
    H5Gclose(iter_group);
}

void OutputManager::write_field_hdf5(hid_t group_id, const std::string& field_name,
                                    PetscScalar** data, PetscInt M, PetscInt N) {
    // Create dataset for field
    hsize_t dims[2] = {(hsize_t)M, (hsize_t)N};
    hid_t dataspace = H5Screate_simple(2, dims, NULL);
    hid_t dataset = H5Dcreate(group_id, field_name.c_str(), H5T_NATIVE_DOUBLE,
                             dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Allocate temporary buffer for 2D data (if needed)
    // Note: This is simplified - in production, use MPI-IO for parallel writes
    std::vector<double> buffer(M * N);
    for (PetscInt j = 0; j < N; j++) {
        for (PetscInt i = 0; i < M; i++) {
            buffer[j * M + i] = data[j][i];
        }
    }
    
    // Write data
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
    
    // Add unitSI attribute (SI units for plasma density: m^-3)
    double unitSI = 1.0; // m^-3 for densities, V for potential
    if (field_name.find("phi") != std::string::npos) {
        unitSI = 1.0; // Volts
    } else if (field_name.find("eps") != std::string::npos) {
        unitSI = 1.602e-19; // eV in Joules
    }
    
    hid_t aspace = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate(dataset, "unitSI", H5T_NATIVE_DOUBLE, aspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &unitSI);
    H5Aclose(attr);
    H5Sclose(aspace);
    
    // Add geometry attribute
    const char* geometry = "cartesian";
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(geometry));
    aspace = H5Screate(H5S_SCALAR);
    attr = H5Acreate(dataset, "geometry", atype, aspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, atype, geometry);
    H5Aclose(attr);
    H5Sclose(aspace);
    H5Tclose(atype);
    
    H5Dclose(dataset);
    H5Sclose(dataspace);
}

PetscErrorCode OutputManager::write_hdf5_output(Vec U, PetscReal time, PetscInt step) {
    PetscErrorCode ierr;
    
    // Create HDF5 filename following OpenPMD convention
    std::ostringstream filename;
    filename << output_dir_ << "/hydroplas_" << std::setfill('0') << std::setw(6) << step << ".h5";
    
    hid_t file_id = create_openpmd_file(filename.str());
    write_openpmd_metadata(file_id, time, step);
    
    // Create meshes group: /data/{step}/meshes
    std::ostringstream meshes_path;
    meshes_path << "/data/" << step << "/meshes";
    hid_t meshes_group = H5Gcreate(file_id, meshes_path.str().c_str(), 
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Get data
    PetscScalar ***u;
    ierr = DMDAVecGetArrayDOF(dm_, U, &u); CHKERRQ(ierr);
    
    PetscInt xs, ys, xm, ym, M, N;
    ierr = DMDAGetCorners(dm_, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
    ierr = DMDAGetInfo(dm_, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
    
    // Write each field
    for (int dof = 0; dof < num_dofs_; ++dof) {
        std::string field_name = (dof < field_names_.size()) ? field_names_[dof] : "field_" + std::to_string(dof);
        
        // Extract 2D slice for this DOF
        std::vector<std::vector<double>> field_data(N, std::vector<double>(M));
        for (PetscInt j = 0; j < N; j++) {
            for (PetscInt i = 0; i < M; i++) {
                field_data[j][i] = u[j][i][dof];
            }
        }
        
        // Convert to pointer array for HDF5
        std::vector<double*> ptrs(N);
        for (PetscInt j = 0; j < N; j++) {
            ptrs[j] = field_data[j].data();
        }
        
        write_field_hdf5(meshes_group, field_name, ptrs.data(), M, N);
    }
    
    ierr = DMDAVecRestoreArrayDOF(dm_, U, &u); CHKERRQ(ierr);
    
    H5Gclose(meshes_group);
    H5Fclose(file_id);
    
    std::cout << "Written HDF5 output: " << filename.str() << std::endl;
    return 0;
}

#endif // USE_HDF5

} // namespace HydroPlas
