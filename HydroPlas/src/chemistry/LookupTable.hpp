#pragma once
#include <vector>
#include <string>
#include <map>

namespace HydroPlas {

class LookupTable {
public:
    void load_from_file(const std::string& filepath);
    double interpolate(double energy, const std::string& column_name) const;

private:
    std::vector<double> energy_grid_;
    std::map<std::string, std::vector<double>> data_columns_;
};

} // namespace HydroPlas
