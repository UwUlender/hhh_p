cmake_minimum_required(VERSION 3.14)
project(HydroPlas LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED PETSc)

find_package(MPI REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)

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
)

target_link_libraries(HydroPlas PRIVATE
    ${PETSC_LIBRARIES}
    MPI::MPI_CXX
    nlohmann_json::nlohmann_json
)

target_include_directories(HydroPlas PRIVATE
    src
)
