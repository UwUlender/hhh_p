cmake_minimum_required(VERSION 3.14)
project(HydroPlas LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
include(FetchContent)

# YAML-CPP
FetchContent_Declare(
    yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG master
)
FetchContent_MakeAvailable(yaml-cpp)

# PETSc
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED PETSc)

# MPI
find_package(MPI REQUIRED)

# HDF5
find_package(HDF5 COMPONENTS C HL REQUIRED)

# Include Directories
include_directories(${PETSC_INCLUDE_DIRS})
include_directories(${MPI_CXX_INCLUDE_PATH})
link_directories(${PETSC_LIBRARY_DIRS})

# Source Files
set(SOURCES
    src/main.cpp
    src/config/ConfigParser.cpp
    src/mesh/RectilinearGrid.cpp
    src/solver/PlasmaSolver.cpp
    src/chemistry/Species.cpp
    src/chemistry/Chemistry.cpp
    src/chemistry/LookupTable.cpp
    src/io/OutputManager.cpp
    src/boundary/BoundaryManager.cpp
    src/numerics/EquationEvaluator.cpp
)

add_executable(HydroPlas ${SOURCES})

target_link_libraries(HydroPlas PRIVATE
    ${PETSC_LIBRARIES}
    MPI::MPI_CXX
    yaml-cpp
    ${HDF5_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
)

target_include_directories(HydroPlas PRIVATE
    src
    ${HDF5_INCLUDE_DIRS}
)
