#pragma once
#include <string>
#include <vector>
#include <hdf5.h>
#include <petscvec.h>
#include "../mesh/RectilinearGrid.hpp"
#include "../config/ConfigParser.hpp"
#include "../chemistry/Chemistry.hpp"

namespace HydroPlas {

class OutputManager {
public:
    OutputManager(const OutputConfig& config, const RectilinearGrid& grid, const Chemistry& chemistry);
    ~OutputManager();

    void write_mesh();
    void write_state(double time, int step, Vec X, int num_species);
    void write_rates(int step, const std::vector<std::vector<double>>& rates, int nx, int ny);
    
    // Restart capability
    void read_state(const std::string& filename, int step, Vec X);

private:
    OutputConfig config_;
    const RectilinearGrid& grid_;
    const Chemistry& chemistry_;
    hid_t file_id_;
    
    void write_dataset(hid_t group_id, const std::string& name, const std::vector<double>& data, const std::vector<hsize_t>& dims);
};

} // namespace HydroPlas
