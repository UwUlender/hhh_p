#!/bin/bash
# Helper script to run benchmarks

if [ ! -d "build" ]; then
    echo "Build directory not found. Please compile first."
    exit 1
fi

echo "Running Benchmark 1: DC Glow Discharge"
./build/HydroPlas config/benchmark_1_dc.json

echo "Running Benchmark 2: RF Discharge"
./build/HydroPlas config/benchmark_2_rf.json

echo "Running Benchmark 3: DBD"
./build/HydroPlas config/benchmark_3_dbd.json
