#pragma once
#include <petsc.h>
#include <petscdmda.h>
#include <string>
#include <vector>
#include "../config/ConfigParser.hpp"

#ifdef USE_HDF5
#include <hdf5.h>
#endif

namespace HydroPlas {

/**
 * @brief OutputManager: Handles data output in multiple formats
 * 
 * Supports:
 * 1. Text-based output (always available)
 * 2. HDF5/OpenPMD-compatible output (if compiled with HDF5)
 * 
 * The OpenPMD schema organizes data as:
 * /data/{iteration}/meshes/{field_name}
 * /data/{iteration}/meshes/{field_name}/@unitSI
 * /data/{iteration}/@time
 * 
 * This enables visualization with ParaView, VisIt, and other OpenPMD-compatible tools.
 */
class OutputManager {
public:
    OutputManager(const SimulationConfig& config, DM dm, int num_excited);
    ~OutputManager();
    
    /**
     * @brief Write all fields at a given time step
     */
    PetscErrorCode write_output(Vec U, PetscReal time, PetscInt step);
    
    /**
     * @brief Set field names for output
     */
    void set_field_names(const std::vector<std::string>& names);
    
private:
    const SimulationConfig& config_;
    DM dm_;
    int num_dofs_;
    int num_excited_;
    std::vector<std::string> field_names_;
    
    std::string output_dir_;
    bool use_hdf5_;
    
    // Text output
    PetscErrorCode write_text_output(Vec U, PetscReal time, PetscInt step);
    
#ifdef USE_HDF5
    // HDF5/OpenPMD output
    PetscErrorCode write_hdf5_output(Vec U, PetscReal time, PetscInt step);
    hid_t create_openpmd_file(const std::string& filename);
    void write_openpmd_metadata(hid_t file_id, PetscReal time, PetscInt step);
    void write_field_hdf5(hid_t group_id, const std::string& field_name, 
                         PetscScalar** data, PetscInt M, PetscInt N);
#endif
};

} // namespace HydroPlas
