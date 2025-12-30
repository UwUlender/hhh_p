#pragma once
#include <string>
#include "../config/ConfigParser.hpp"

namespace HydroPlas {

class BolsigInterface {
public:
    static void run_bolsig(const ChemistryConfig& config, const std::string& output_file);
};

} // namespace HydroPlas
