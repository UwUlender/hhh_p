#pragma once

#include <petscdm.h>
#include <petscdmda.h>
#include "../config/ConfigParser.hpp"

namespace HydroPlas {

class MeshGenerator {
public:
    explicit MeshGenerator(const DomainConfig& config);
    ~MeshGenerator();

    PetscErrorCode create_dm(DM* dm, int num_dofs);

private:
    DomainConfig config_;
};

} // namespace HydroPlas
