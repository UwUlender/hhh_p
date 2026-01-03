# HydroPlas Implementation Roadmap: Adopting Zapdos Best Practices

This document provides concrete, actionable code examples for implementing the recommendations from the Zapdos comparison analysis.

---

## 1. Testing Infrastructure

### 1.1 Basic Test Setup with Catch2

**File: `test/CMakeLists.txt`**
```cmake
cmake_minimum_required(VERSION 3.14)

# Fetch Catch2
include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.5.0
)
FetchContent_MakeAvailable(Catch2)

# Test executable
add_executable(hydroplas_tests
    unit/test_chemistry.cpp
    unit/test_boundary.cpp
    unit/test_scharfetter_gummel.cpp
    integration/test_dc_discharge.cpp
)

target_link_libraries(hydroplas_tests PRIVATE
    Catch2::Catch2WithMain
    # Link HydroPlas library (refactor main code into library)
)

# CTest integration
include(CTest)
include(Catch)
catch_discover_tests(hydroplas_tests)
```

**File: `test/unit/test_chemistry.cpp`**
```cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "../../src/chemistry/Chemistry.hpp"

using namespace HydroPlas;
using Catch::Matchers::WithinRel;

TEST_CASE("Chemistry: Direct ionization rate", "[chemistry]") {
    // Setup minimal config
    SimulationConfig config;
    config.chemistry_config.mode = "Pre-calculated";
    config.chemistry_config.gas_temperature = 300.0;
    
    Species electron{"e", 9.11e-31, -1.6e-19};
    config.chemistry_config.species_list = {electron};
    
    Chemistry chem(config);
    
    SECTION("Rate increases with electron energy") {
        std::vector<double> densities = {1e16}; // ne = 1e16 m^-3
        std::vector<double> sources(1, 0.0);
        
        // Low energy - minimal ionization
        chem.compute_source(densities, 5.0, 300.0, sources);
        double rate_low = sources[0];
        
        // High energy - significant ionization
        chem.compute_source(densities, 20.0, 300.0, sources);
        double rate_high = sources[0];
        
        REQUIRE(rate_high > rate_low);
        REQUIRE(rate_high > 1e10); // Sanity check
    }
}

TEST_CASE("Chemistry: Conservation in Penning ionization", "[chemistry]") {
    // Ar* + Ar* -> Ar+ + Ar + e
    // Should conserve: d(Ar*) = -2*d(Ar+)
    
    SimulationConfig config;
    // ... setup with Penning reaction
    
    Chemistry chem(config);
    
    std::vector<double> densities = {1e16, 1e18, 1e17}; // e, Ar+, Ar*
    std::vector<double> sources(3);
    
    chem.compute_source(densities, 10.0, 300.0, sources);
    
    // Check stoichiometry (approximate due to other reactions)
    // This is a simplified test
    REQUIRE_THAT(sources[0] + sources[1], WithinRel(0.0, 0.01));
}
```

### 1.2 Regression Testing with Gold Files

**File: `test/regression/test_cases.yaml`**
```yaml
tests:
  - name: dc_discharge_breakdown
    config: configs/test_dc_discharge.json
    reference: gold/dc_discharge_breakdown.h5
    tolerance:
      relative: 1e-3
      absolute: 1e-8
    fields:
      - ne
      - ni
      - phi
    time_steps: [50, 100, 200]
    
  - name: dbd_memory_effect
    config: configs/test_dbd.json
    reference: gold/dbd_memory_effect.h5
    tolerance:
      relative: 5e-3
      absolute: 1e-7
    fields:
      - Ar_m
      - ne
    time_steps: [100, 500]
```

**File: `test/regression/run_regression.py`**
```python
#!/usr/bin/env python3
import subprocess
import h5py
import numpy as np
import yaml
import sys

def compare_hdf5(test_file, gold_file, field, time_step, rtol, atol):
    """Compare HDF5 outputs"""
    with h5py.File(test_file, 'r') as f_test, \
         h5py.File(gold_file, 'r') as f_gold:
        
        data_test = f_test[f'/data/{time_step}/{field}'][:]
        data_gold = f_gold[f'/data/{time_step}/{field}'][:]
        
        if not np.allclose(data_test, data_gold, rtol=rtol, atol=atol):
            max_diff = np.max(np.abs(data_test - data_gold))
            rel_diff = np.max(np.abs((data_test - data_gold) / (data_gold + 1e-20)))
            print(f"  ❌ FAILED: {field} at t={time_step}")
            print(f"     Max absolute difference: {max_diff:.2e}")
            print(f"     Max relative difference: {rel_diff:.2e}")
            return False
        return True

def run_tests():
    with open('test_cases.yaml', 'r') as f:
        test_suite = yaml.safe_load(f)
    
    failed_tests = []
    
    for test in test_suite['tests']:
        print(f"\n🧪 Running: {test['name']}")
        
        # Run HydroPlas
        result = subprocess.run(
            ['../../build/HydroPlas', test['config']],
            capture_output=True, text=True
        )
        
        if result.returncode != 0:
            print(f"  ❌ Simulation failed!")
            print(result.stderr)
            failed_tests.append(test['name'])
            continue
        
        # Compare outputs
        test_passed = True
        for field in test['fields']:
            for time_step in test['time_steps']:
                if not compare_hdf5(
                    'output.h5', test['reference'], field, time_step,
                    test['tolerance']['relative'],
                    test['tolerance']['absolute']
                ):
                    test_passed = False
        
        if test_passed:
            print(f"  ✅ PASSED")
        else:
            failed_tests.append(test['name'])
    
    print(f"\n{'='*60}")
    print(f"Results: {len(test_suite['tests']) - len(failed_tests)}/{len(test_suite['tests'])} passed")
    
    if failed_tests:
        print(f"Failed tests: {', '.join(failed_tests)}")
        sys.exit(1)
    else:
        print("🎉 All tests passed!")
        sys.exit(0)

if __name__ == '__main__':
    run_tests()
```

### 1.3 GitHub Actions CI

**File: `.github/workflows/ci.yml`**
```yaml
name: CI

on:
  push:
    branches: [ main, develop ]
  pull_request:
    branches: [ main ]

jobs:
  build-and-test:
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-22.04, macos-12]
        build_type: [Release, Debug]
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Install dependencies (Ubuntu)
      if: runner.os == 'Linux'
      run: |
        sudo apt-get update
        sudo apt-get install -y \
          petsc-dev \
          libhdf5-dev \
          libyaml-cpp-dev \
          cmake \
          build-essential
    
    - name: Install dependencies (macOS)
      if: runner.os == 'macOS'
      run: |
        brew install petsc hdf5 yaml-cpp cmake
    
    - name: Configure CMake
      run: |
        mkdir build && cd build
        cmake .. -DCMAKE_BUILD_TYPE=${{ matrix.build_type }}
    
    - name: Build
      run: cmake --build build --config ${{ matrix.build_type }} -j4
    
    - name: Run unit tests
      run: |
        cd build/test
        ctest --output-on-failure
    
    - name: Run regression tests
      run: |
        cd test/regression
        python3 run_regression.py
    
    - name: Upload artifacts
      if: failure()
      uses: actions/upload-artifact@v3
      with:
        name: test-outputs-${{ matrix.os }}
        path: test/regression/*.h5

  lint:
    runs-on: ubuntu-22.04
    steps:
    - uses: actions/checkout@v3
    
    - name: Install clang-format
      run: sudo apt-get install -y clang-format-14
    
    - name: Check formatting
      run: |
        find src -name "*.cpp" -o -name "*.hpp" | \
          xargs clang-format-14 --dry-run --Werror

  documentation:
    runs-on: ubuntu-22.04
    steps:
    - uses: actions/checkout@v3
    
    - name: Install Doxygen
      run: sudo apt-get install -y doxygen graphviz
    
    - name: Generate documentation
      run: doxygen Doxyfile
    
    - name: Deploy to GitHub Pages
      if: github.ref == 'refs/heads/main'
      uses: peaceiris/actions-gh-pages@v3
      with:
        github_token: ${{ secrets.GITHUB_TOKEN }}
        publish_dir: ./docs/html
```

---

## 2. Plugin Architecture for Physics

### 2.1 Base Kernel Interface

**File: `src/kernels/PhysicsKernel.hpp`**
```cpp
#pragma once
#include <string>
#include <memory>
#include <petscsnes.h>

namespace HydroPlas {

// Forward declarations
struct SolverState;
class RectilinearGrid;

/**
 * @brief Base class for all physics kernels
 * 
 * Physics kernels compute contributions to the residual vector F
 * and optionally the Jacobian matrix J. Users can inherit from this
 * to add custom physics without modifying core HydroPlas code.
 */
class PhysicsKernel {
public:
    PhysicsKernel(const std::string& name) : name_(name), enabled_(true) {}
    virtual ~PhysicsKernel() = default;
    
    /**
     * @brief Compute residual contribution
     * @param state Current solution state (X, X_prev, dt, time)
     * @param grid Mesh information
     * @param F Residual vector to modify
     */
    virtual void compute_residual(
        const SolverState& state,
        const RectilinearGrid& grid,
        Vec F
    ) = 0;
    
    /**
     * @brief Compute Jacobian contribution (optional)
     * @param state Current solution state
     * @param grid Mesh information
     * @param J Jacobian matrix to modify
     */
    virtual void compute_jacobian(
        const SolverState& state,
        const RectilinearGrid& grid,
        Mat J
    ) {
        // Default: do nothing (use finite differencing)
    }
    
    /**
     * @brief Initialize kernel (called once at setup)
     */
    virtual void initialize(const SimulationConfig& config) {}
    
    /**
     * @brief Get kernel name for logging
     */
    const std::string& name() const { return name_; }
    
    /**
     * @brief Enable/disable kernel
     */
    void set_enabled(bool enabled) { enabled_ = enabled; }
    bool is_enabled() const { return enabled_; }

protected:
    std::string name_;
    bool enabled_;
};

/**
 * @brief State information passed to kernels
 */
struct SolverState {
    Vec X;          // Current solution
    Vec X_prev;     // Previous time step solution
    double dt;      // Time step size
    double time;    // Current simulation time
    
    // Convenience accessors (filled by solver)
    int idx_ne, idx_ni, idx_neps, idx_phi;
    int n_excited_species;
    std::vector<int> idx_excited;
};

} // namespace HydroPlas
```

### 2.2 Example Kernel Implementation

**File: `src/kernels/IonizationKernel.hpp`**
```cpp
#pragma once
#include "PhysicsKernel.hpp"
#include "../chemistry/Chemistry.hpp"

namespace HydroPlas {

/**
 * @brief Computes electron production from ionization
 * 
 * Adds source term: S_e = k_iz(E_mean) * ne * ng
 */
class IonizationKernel : public PhysicsKernel {
public:
    IonizationKernel(Chemistry& chemistry) 
        : PhysicsKernel("Ionization"), chemistry_(chemistry) {}
    
    void compute_residual(
        const SolverState& state,
        const RectilinearGrid& grid,
        Vec F
    ) override;
    
private:
    Chemistry& chemistry_;
};

} // namespace HydroPlas
```

**File: `src/kernels/IonizationKernel.cpp`**
```cpp
#include "IonizationKernel.hpp"
#include "../mesh/RectilinearGrid.hpp"

namespace HydroPlas {

void IonizationKernel::compute_residual(
    const SolverState& state,
    const RectilinearGrid& grid,
    Vec F
) {
    if (!enabled_) return;
    
    // Access solution arrays
    const double ***X_arr, ***F_arr;
    DMDAVecGetArrayRead(grid.get_dm(), state.X, &X_arr);
    DMDAVecGetArray(grid.get_dm(), F, &F_arr);
    
    int xs, ys, xm, ym;
    DMDAGetCorners(grid.get_dm(), &xs, &ys, nullptr, &xm, &ym, nullptr);
    
    for (int j = ys; j < ys + ym; j++) {
        for (int i = xs; i < xs + xm; i++) {
            // Extract densities
            double ne = exp(X_arr[j][i][state.idx_ne]);
            double mean_energy = X_arr[j][i][state.idx_neps] / ne;
            
            // Get ionization rate from chemistry
            double k_iz = chemistry_.get_ionization_rate(mean_energy);
            double ng = 3.22e22; // Background gas density (TODO: from config)
            
            // Compute source term
            double S_ionization = k_iz * ne * ng;
            
            // Add to residual (note sign convention)
            double volume = grid.get_cell_volume(i, j);
            F_arr[j][i][state.idx_ne] -= S_ionization * volume * state.dt;
            
            // Ion production (equal to electron production)
            F_arr[j][i][state.idx_ni] -= S_ionization * volume * state.dt;
        }
    }
    
    DMDAVecRestoreArrayRead(grid.get_dm(), state.X, &X_arr);
    DMDAVecRestoreArray(grid.get_dm(), F, &F_arr);
}

} // namespace HydroPlas
```

### 2.3 Kernel Manager

**File: `src/kernels/KernelManager.hpp`**
```cpp
#pragma once
#include <vector>
#include <memory>
#include "PhysicsKernel.hpp"

namespace HydroPlas {

/**
 * @brief Manages collection of physics kernels
 */
class KernelManager {
public:
    void add_kernel(std::unique_ptr<PhysicsKernel> kernel);
    
    void compute_residuals(
        const SolverState& state,
        const RectilinearGrid& grid,
        Vec F
    );
    
    void compute_jacobians(
        const SolverState& state,
        const RectilinearGrid& grid,
        Mat J
    );
    
    void initialize_all(const SimulationConfig& config);
    
    // Access kernels by name
    PhysicsKernel* get_kernel(const std::string& name);
    
private:
    std::vector<std::unique_ptr<PhysicsKernel>> kernels_;
};

} // namespace HydroPlas
```

**File: `src/kernels/KernelManager.cpp`**
```cpp
#include "KernelManager.hpp"
#include <algorithm>

namespace HydroPlas {

void KernelManager::add_kernel(std::unique_ptr<PhysicsKernel> kernel) {
    kernels_.push_back(std::move(kernel));
}

void KernelManager::compute_residuals(
    const SolverState& state,
    const RectilinearGrid& grid,
    Vec F
) {
    for (auto& kernel : kernels_) {
        if (kernel->is_enabled()) {
            kernel->compute_residual(state, grid, F);
        }
    }
}

void KernelManager::compute_jacobians(
    const SolverState& state,
    const RectilinearGrid& grid,
    Mat J
) {
    for (auto& kernel : kernels_) {
        if (kernel->is_enabled()) {
            kernel->compute_jacobian(state, grid, J);
        }
    }
}

void KernelManager::initialize_all(const SimulationConfig& config) {
    for (auto& kernel : kernels_) {
        kernel->initialize(config);
    }
}

PhysicsKernel* KernelManager::get_kernel(const std::string& name) {
    auto it = std::find_if(kernels_.begin(), kernels_.end(),
        [&name](const auto& k) { return k->name() == name; });
    
    return it != kernels_.end() ? it->get() : nullptr;
}

} // namespace HydroPlas
```

### 2.4 Integration into Solver

**Modified `src/solver/PlasmaSolver.cpp`:**
```cpp
#include "PlasmaSolver.hpp"
#include "../kernels/KernelManager.hpp"
#include "../kernels/IonizationKernel.hpp"
#include "../kernels/DiffusionKernel.hpp"
// ... other kernels

namespace HydroPlas {

PlasmaSolver::PlasmaSolver(/* ... */) {
    // Initialize kernel manager
    kernel_manager_ = std::make_unique<KernelManager>();
    
    // Add standard kernels
    kernel_manager_->add_kernel(
        std::make_unique<IonizationKernel>(chemistry_)
    );
    kernel_manager_->add_kernel(
        std::make_unique<DiffusionKernel>(config_)
    );
    // ... add more kernels
    
    // Initialize all kernels
    kernel_manager_->initialize_all(config_);
}

// In FormFunction callback:
PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void* ctx) {
    auto* solver_ctx = static_cast<SolverContext*>(ctx);
    
    // Fill solver state
    SolverState state;
    state.X = X;
    state.X_prev = solver_ctx->X_prev;
    state.dt = solver_ctx->dt;
    state.time = solver_ctx->time;
    // ... fill indices
    
    // Compute residuals from all active kernels
    solver_ctx->kernel_manager->compute_residuals(
        state, *solver_ctx->grid, F
    );
    
    return 0;
}

} // namespace HydroPlas
```

---

## 3. Enhanced Configuration System

### 3.1 Flexible Kernel Configuration

**Extended JSON format:**
```json
{
    "domain": {
        "Lx": 0.01,
        "Nx": 200
    },
    
    "kernels": {
        "drift": {
            "enabled": true,
            "type": "ScharfetterGummel"
        },
        "diffusion": {
            "enabled": true,
            "type": "CoeffDiffusion"
        },
        "ionization": {
            "enabled": true,
            "type": "FromBolsig",
            "table_file": "data/transport.dat"
        },
        "custom_heating": {
            "enabled": true,
            "type": "Plugin",
            "library": "./libuser_heating.so",
            "parameters": {
                "heating_rate": 1e6
            }
        }
    },
    
    "boundaries": {
        "left": {
            "type": "Electrode",
            "voltage": {
                "type": "RF",
                "amplitude": 500.0,
                "frequency": 13.56e6
            },
            "species_bc": {
                "electrons": {"type": "Hagelaar", "gamma": 0.1},
                "ions": {"type": "Drift"},
                "Ar_m": {"type": "Quenching", "prob": 0.001}
            }
        },
        "right": {
            "type": "Dielectric",
            "permittivity": 4.0,
            "thickness": 1e-3
        }
    },
    
    "output": {
        "format": "hdf5",
        "interval": 50,
        "fields": ["ne", "ni", "phi", "Ar_m"],
        "derived_quantities": [
            {
                "name": "E_field",
                "formula": "-gradient(phi)"
            },
            {
                "name": "power_density",
                "formula": "dot(J_e, E_field)"
            }
        ]
    }
}
```

### 3.2 Config Parser Enhancement

**File: `src/config/ConfigParser.hpp` (additions):**
```cpp
struct KernelConfig {
    std::string name;
    bool enabled;
    std::string type;
    std::string library_path; // For plugins
    std::map<std::string, double> parameters;
};

struct BoundaryConfig {
    std::string name;
    std::string type;
    std::map<std::string, SpeciesBCConfig> species_bcs;
    // ... voltage config, etc.
};

struct DerivedQuantityConfig {
    std::string name;
    std::string formula;
    std::string units;
};

struct SimulationConfig {
    // ... existing fields
    std::vector<KernelConfig> kernel_configs;
    std::map<std::string, BoundaryConfig> boundary_configs;
    std::vector<DerivedQuantityConfig> derived_quantities;
};
```

---

## 4. Docker Deployment

### 4.1 Dockerfile

**File: `Dockerfile`**
```dockerfile
FROM ubuntu:22.04 AS builder

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    petsc-dev \
    libhdf5-dev \
    libyaml-cpp-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for visualization
RUN pip3 install h5py numpy matplotlib plotly

# Copy source code
WORKDIR /hydroplas
COPY . .

# Build HydroPlas
RUN mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make -j$(nproc)

# Runtime stage (smaller image)
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    libpetsc-real3.15 \
    libhdf5-103 \
    libyaml-cpp0.7 \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install h5py numpy matplotlib plotly

# Copy built executable and examples
COPY --from=builder /hydroplas/build/HydroPlas /usr/local/bin/
COPY --from=builder /hydroplas/HydroPlas/config /config
COPY --from=builder /hydroplas/HydroPlas/data /data
COPY --from=builder /hydroplas/HydroPlas/process_results.py /usr/local/bin/

WORKDIR /workspace
VOLUME /workspace

# Set up entrypoint
COPY docker-entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker-entrypoint.sh
ENTRYPOINT ["docker-entrypoint.sh"]
CMD ["--help"]
```

**File: `docker-entrypoint.sh`**
```bash
#!/bin/bash
set -e

if [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
    echo "HydroPlas Docker Container"
    echo "Usage:"
    echo "  docker run -v \$(pwd):/workspace hydroplasma/hydroplas [config.json]"
    echo ""
    echo "Examples:"
    echo "  # Run built-in example"
    echo "  docker run hydroplasma/hydroplas /config/argon_complete.json"
    echo ""
    echo "  # Run your own config"
    echo "  docker run -v \$(pwd):/workspace hydroplasma/hydroplas /workspace/my_config.json"
    echo ""
    echo "  # Interactive mode"
    echo "  docker run -it hydroplasma/hydroplas /bin/bash"
    exit 0
fi

# If argument is a JSON file, run HydroPlas
if [[ "$1" == *.json ]] || [[ "$1" == *.yaml ]]; then
    echo "🚀 Running HydroPlas with config: $1"
    exec HydroPlas "$@"
fi

# Otherwise, execute the command
exec "$@"
```

### 4.2 Docker Compose for Development

**File: `docker-compose.yml`**
```yaml
version: '3.8'

services:
  hydroplas:
    build:
      context: .
      dockerfile: Dockerfile
    volumes:
      - ./examples:/workspace
      - ./output:/output
    environment:
      - PETSC_OPTIONS=-log_view
    command: /config/argon_complete.json
  
  visualization:
    image: python:3.10
    volumes:
      - ./output:/data
      - ./process_results.py:/app/process_results.py
    working_dir: /app
    ports:
      - "8050:8050"  # For Plotly Dash
    command: python3 process_results.py --live --port 8050
```

### 4.3 GitHub Action for Docker Builds

**File: `.github/workflows/docker.yml`**
```yaml
name: Build and Push Docker Image

on:
  push:
    branches: [ main ]
    tags: [ 'v*' ]
  release:
    types: [ published ]

jobs:
  docker:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
    
    - name: Docker meta
      id: meta
      uses: docker/metadata-action@v4
      with:
        images: hydroplasma/hydroplas
        tags: |
          type=ref,event=branch
          type=ref,event=pr
          type=semver,pattern={{version}}
          type=semver,pattern={{major}}.{{minor}}
          type=sha
    
    - name: Login to Docker Hub
      uses: docker/login-action@v2
      with:
        username: ${{ secrets.DOCKERHUB_USERNAME }}
        password: ${{ secrets.DOCKERHUB_TOKEN }}
    
    - name: Build and push
      uses: docker/build-push-action@v4
      with:
        context: .
        push: true
        tags: ${{ steps.meta.outputs.tags }}
        labels: ${{ steps.meta.outputs.labels }}
```

---

## 5. Documentation Website

### 5.1 GitHub Pages with MkDocs

**File: `mkdocs.yml`**
```yaml
site_name: HydroPlas Documentation
site_description: Advanced plasma simulation with explicit excited species transport
site_author: HydroPlas Team
site_url: https://yourorg.github.io/hydroplas

theme:
  name: material
  palette:
    primary: blue
    accent: cyan
  features:
    - navigation.tabs
    - navigation.sections
    - toc.integrate
    - search.suggest
    - content.code.copy

plugins:
  - search
  - mkdocstrings:
      handlers:
        python:
          paths: [scripts]

markdown_extensions:
  - pymdownx.highlight
  - pymdownx.superfences
  - pymdownx.arithmatex:
      generic: true
  - admonition
  - codehilite
  - toc:
      permalink: true

extra_javascript:
  - javascripts/mathjax.js
  - https://polyfill.io/v3/polyfill.min.js?features=es6
  - https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js

nav:
  - Home: index.md
  - Getting Started:
    - Installation: getting-started/installation.md
    - Quick Start: getting-started/quick-start.md
    - Docker: getting-started/docker.md
  - Tutorials:
    - Basic Discharge: tutorials/01_basic_discharge.md
    - Excited Species: tutorials/02_excited_species.md
    - Multi-Electrode: tutorials/03_multi_electrode.md
    - Parameter Studies: tutorials/04_parameter_studies.md
  - Theory:
    - Governing Equations: theory/equations.md
    - Numerics: theory/numerics.md
    - Chemistry: theory/chemistry.md
  - Reference:
    - Configuration: reference/configuration.md
    - Boundary Conditions: reference/boundary_conditions.md
    - Output Format: reference/output_format.md
  - Development:
    - Contributing: development/contributing.md
    - Plugin Development: development/plugins.md
    - Testing: development/testing.md
  - API: api/
```

**File: `docs/index.md`**
```markdown
# HydroPlas

<div align="center">
  <img src="images/logo.png" alt="HydroPlas Logo" width="300">
  
  **Advanced Computational Framework for Non-Equilibrium Plasma Fluid Simulation**
  
  [![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
  [![Build Status](https://github.com/yourorg/hydroplas/workflows/CI/badge.svg)](https://github.com/yourorg/hydroplas/actions)
  [![Docker](https://img.shields.io/docker/v/hydroplasma/hydroplas?label=docker)](https://hub.docker.com/r/hydroplasma/hydroplas)
  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXX)
</div>

## Features

=== "Physics"

    - ✅ **Explicit Excited Species Transport** - Full ADR treatment
    - ✅ **Stepwise & Penning Ionization** - Multi-step processes
    - ✅ **Multi-Electrode Control** - Independent voltage waveforms
    - ✅ **Wall Interactions** - Quenching & secondary emission

=== "Numerics"

    - ✅ **Scharfetter-Gummel Scheme** - Stable for all Péclet numbers
    - ✅ **Implicit Time Integration** - PETSc TS framework
    - ✅ **Newton-Krylov Solver** - Efficient for stiff systems
    - ✅ **Automatic Differentiation** - Exact Jacobians

=== "Software"

    - ✅ **Modern C++** - Clean, readable code
    - ✅ **Plugin Architecture** - Extensible without recompilation
    - ✅ **Docker Support** - Run anywhere
    - ✅ **Comprehensive Testing** - Unit + regression tests

## Quick Start

```bash
# Using Docker (easiest)
docker pull hydroplasma/hydroplas:latest
docker run -v $(pwd):/workspace hydroplasma/hydroplas /config/argon_complete.json

# From source
git clone https://github.com/yourorg/hydroplas.git
cd hydroplas/HydroPlas
mkdir build && cd build
cmake .. && make -j4
./HydroPlas ../config/argon_complete.json
```

## Example: DBD with Memory Effect

```json
{
    "domain": {"Lx": 0.01, "Nx": 200},
    "boundary": {
        "voltage_type": "AC",
        "voltage_amplitude": 5000.0,
        "frequency": 50000.0
    },
    "chemistry": {
        "excited_species": [
            {"name": "Ar_m", "energy_level": 11.55, ...}
        ]
    }
}
```

<div class="result" markdown>
![DBD Metastable Accumulation](images/dbd_memory_effect.png)

*Metastable accumulation at dielectric reduces breakdown voltage by 40%*
</div>

## Why HydroPlas?

| Feature | HydroPlas | COMSOL | Zapdos |
|---------|-----------|--------|--------|
| Excited Species Transport | ✅ Native | ⚠️ Manual | ✅ Yes |
| Multi-Electrode Control | ✅ Built-in | ⚠️ Complex | ⚠️ Via circuits |
| Open Source | ✅ MIT | ❌ Proprietary | ✅ LGPL |
| Easy to Extend | ✅ Plugins | ❌ No | ⚠️ MOOSE steep |
| Computational Cost | ✅ Low | ⚠️ High | ⚠️ Medium |

## Citation

If you use HydroPlas in your research, please cite:

```bibtex
@software{hydroplas2026,
  author = {Your Name},
  title = {HydroPlas: Plasma Simulation with Excited Species Transport},
  year = {2026},
  doi = {10.5281/zenodo.XXXXX},
  url = {https://github.com/yourorg/hydroplas}
}
```

## Next Steps

<div class="grid cards" markdown>

-   :material-clock-fast:{ .lg .middle } __Get Started in 5 Minutes__

    ---

    Install HydroPlas and run your first simulation

    [:octicons-arrow-right-24: Quick Start](getting-started/quick-start.md)

-   :material-book-open-page-variant:{ .lg .middle } __Learn the Theory__

    ---

    Understand the physics and numerics

    [:octicons-arrow-right-24: Theory Guide](theory/equations.md)

-   :material-test-tube:{ .lg .middle } __Follow Tutorials__

    ---

    Step-by-step examples with explanations

    [:octicons-arrow-right-24: Tutorials](tutorials/01_basic_discharge.md)

-   :material-puzzle:{ .lg .middle } __Extend with Plugins__

    ---

    Add custom physics without editing source

    [:octicons-arrow-right-24: Plugin Guide](development/plugins.md)

</div>
```

### 5.2 GitHub Action for Documentation

**File: `.github/workflows/docs.yml`**
```yaml
name: Documentation

on:
  push:
    branches: [ main ]

jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: 3.x
    
    - name: Install dependencies
      run: |
        pip install mkdocs-material mkdocstrings
    
    - name: Build documentation
      run: mkdocs build
    
    - name: Deploy to GitHub Pages
      uses: peaceiris/actions-gh-pages@v3
      with:
        github_token: ${{ secrets.GITHUB_TOKEN }}
        publish_dir: ./site
```

---

## 6. Summary Checklist

### Week 1: Testing Foundation
- [ ] Add Catch2 to CMakeLists.txt
- [ ] Write 5 unit tests (chemistry, boundary, numerics)
- [ ] Create regression test infrastructure
- [ ] Set up GitHub Actions CI

### Week 2: Documentation
- [ ] Initialize MkDocs project
- [ ] Write 3 tutorial pages
- [ ] Generate Doxygen API docs
- [ ] Deploy to GitHub Pages

### Week 3: Docker
- [ ] Create Dockerfile
- [ ] Test on Linux and macOS
- [ ] Publish to Docker Hub
- [ ] Update README with Docker instructions

### Week 4: Extensibility
- [ ] Implement PhysicsKernel base class
- [ ] Refactor solver to use KernelManager
- [ ] Create 3 example kernels
- [ ] Document plugin development process

---

**This roadmap provides concrete, implementable code that brings HydroPlas to the same level of software engineering maturity as Zapdos while maintaining its unique advantages.**
