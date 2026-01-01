#pragma once
#include <string>
#include <vector>
#include <hdf5.h>
#include <petscvec.h>
#include "../mesh/RectilinearGrid.hpp"
#include "../config/ConfigParser.hpp"

namespace HydroPlas {

class OutputManager {
public:
    OutputManager(const OutputConfig& config, const RectilinearGrid& grid);
    ~OutputManager();

    void write_mesh();
    void write_state(double time, int step, Vec X, int num_species);
    // write_rates( ... )

private:
    OutputConfig config_;
    const RectilinearGrid& grid_;
    hid_t file_id_;
    
    void write_dataset(hid_t group_id, const std::string& name, const std::vector<double>& data, const std::vector<hsize_t>& dims);
};

} // namespace HydroPlas
