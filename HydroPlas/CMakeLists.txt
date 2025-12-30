cmake_minimum_required(VERSION 3.14)
project(HydroPlas LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED PETSc)

find_package(MPI REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)

# Optional: HDF5 for advanced output (OpenPMD-compatible)
find_package(HDF5 COMPONENTS C HL)
if(HDF5_FOUND)
    message(STATUS "HDF5 found - enabling OpenPMD output")
    add_definitions(-DUSE_HDF5)
else()
    message(STATUS "HDF5 not found - using text output only")
endif()

include_directories(${PETSC_INCLUDE_DIRS})
include_directories(${MPI_CXX_INCLUDE_PATH})
link_directories(${PETSC_LIBRARY_DIRS})

add_executable(HydroPlas
    src/main.cpp
    src/config/ConfigParser.cpp
    src/mesh/MeshGenerator.cpp
    src/solver/Solver.cpp
    src/boundary/BoundaryManager.cpp
    src/chemistry/LookupTable.cpp
    src/chemistry/BolsigInterface.cpp
    src/chemistry/ReactionHandler.cpp
    src/io/OutputManager.cpp
)

target_link_libraries(HydroPlas PRIVATE
    ${PETSC_LIBRARIES}
    MPI::MPI_CXX
    nlohmann_json::nlohmann_json
)

if(HDF5_FOUND)
    target_include_directories(HydroPlas PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(HydroPlas PRIVATE ${HDF5_LIBRARIES})
endif()

target_include_directories(HydroPlas PRIVATE
    src
)
