#!/bin/bash
# Quick rebuild script for PETSc configuration fix
# This script rebuilds only the modified ConfigParser files

set -e  # Exit on error

echo "=================================================="
echo "HydroPlas PETSc Configuration Fix - Rebuild Script"
echo "=================================================="
echo ""

# Check if we're in the HydroPlas directory
if [ ! -f "CMakeLists.txt" ]; then
    echo "ERROR: Please run this script from the HydroPlas root directory"
    exit 1
fi

echo "Cleaning old build..."
rm -rf build
mkdir -p build

echo "Configuring with CMake..."
cd build
cmake ..

echo ""
echo "Building HydroPlas..."
make -j$(nproc)

if [ $? -eq 0 ]; then
    echo ""
    echo "=================================================="
    echo "✓ Build successful!"
    echo "=================================================="
    echo ""
    echo "The PETSc configuration fix has been applied."
    echo ""
    echo "Test the program with:"
    echo "  ./build/HydroPlas --config config/test_petsc_fix.yaml"
    echo ""
    echo "or with the default config:"
    echo "  ./build/HydroPlas --config config/default_config.yaml"
    echo ""
else
    echo ""
    echo "=================================================="
    echo "✗ Build failed!"
    echo "=================================================="
    echo ""
    echo "Please check the error messages above."
    echo "You may need to install missing dependencies."
    exit 1
fi
