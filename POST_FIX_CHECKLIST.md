# Post-Fix Checklist

## Step-by-Step Application Guide

### ☐ Step 1: Backup Your Current Code
```bash
cd ~/src/own_code
cp -r HydroPlas HydroPlas.backup
```

### ☐ Step 2: Apply the Fix

**Option A: Copy the fixed file**
```bash
cp /workspace/HydroPlas/src/solver/PlasmaSolver.cpp \
   ~/src/own_code/HydroPlas/src/solver/
```

**Option B: Use git to apply changes** (if in git repo)
```bash
cd ~/src/own_code/HydroPlas
git checkout cursor/segmentation-violation-error-51cc
```

### ☐ Step 3: Clean Previous Build
```bash
cd ~/src/own_code/HydroPlas/build
make clean
# Or completely remove build directory:
# cd ~/src/own_code/HydroPlas && rm -rf build && mkdir build
```

### ☐ Step 4: Rebuild
```bash
cd ~/src/own_code/HydroPlas/build
cmake ..
make -j$(nproc)
```

**Expected output:**
- No compilation errors
- Executable created: `HydroPlas`

### ☐ Step 5: Verify Executable Exists
```bash
ls -lh ~/src/own_code/HydroPlas/build/HydroPlas
# Or if installed:
ls -lh ~/src/own_code/HydroPlas/HydroPlas
```

### ☐ Step 6: Test Run
```bash
cd ~/src/own_code/HydroPlas
./HydroPlas --config config/main_test.yaml
```

**Expected behavior:**
- ✓ No segmentation fault
- ✓ Prints "Starting Simulation..."
- ✓ Shows time step progress
- ✓ Creates output file (output.h5 or hydroplas_out.h5)

### ☐ Step 7: Verify Initial Conditions
Check the output to ensure initial conditions match your config:
```bash
# If you have h5dump or python with h5py:
h5dump -d "/Step_0/e_density" output.h5 | head -20
# Should show values near 1e14
```

### ☐ Step 8: Short Test Run
Modify config for a very short run to test stability:
```yaml
solver:
  time_step: 1.8436578e-11
  end_time: 1.8436578e-09  # Just 100 steps
```

Then run again to verify it completes without errors.

## Troubleshooting

### Problem: Compilation Errors

**Solution 1: Check file was copied correctly**
```bash
grep "if (i == 0 || i == nx-1) continue;" \
     ~/src/own_code/HydroPlas/src/solver/PlasmaSolver.cpp
# Should find the line (around line 452)
```

**Solution 2: Check for proper includes**
```bash
head -20 ~/src/own_code/HydroPlas/src/solver/PlasmaSolver.cpp
# Should show #include statements
```

**Solution 3: Verify CMake configuration**
```bash
cd ~/src/own_code/HydroPlas/build
cmake .. 2>&1 | tee cmake_log.txt
# Check cmake_log.txt for errors
```

### Problem: Still Getting Segmentation Fault

**Check 1: Verify fix was applied**
```bash
grep -n "Skip boundary cells" \
     ~/src/own_code/HydroPlas/src/solver/PlasmaSolver.cpp
# Should show line with comment "// Skip boundary cells - they get Dirichlet BC"
```

**Check 2: Run with debugger**
```bash
gdb ./HydroPlas
(gdb) run --config config/main_test.yaml
# If it crashes:
(gdb) backtrace
# This will show where the crash occurred
```

**Check 3: Enable PETSc error checking**
```bash
./HydroPlas --config config/main_test.yaml -start_in_debugger
```

### Problem: Wrong Initial Conditions

**Check 1: Verify config file is correct**
```bash
grep -A 10 "initial_condition:" ~/src/own_code/HydroPlas/config/main_test.yaml
```

**Check 2: Add debug output** (temporary)
Add to setup_dofs() after applying initial conditions:
```cpp
PetscPrintf(PETSC_COMM_WORLD, "Applied initial conditions from config\n");
```

### Problem: Compilation Takes Too Long

**Solution: Use parallel compilation**
```bash
make -j$(nproc)  # Uses all CPU cores
# Or specify number:
make -j4         # Uses 4 cores
```

## Verification Checklist

After successful rebuild and run, verify:

- ☐ **No segmentation fault** during initialization
- ☐ **No segmentation fault** during first time step
- ☐ **Initial densities** match config file values
- ☐ **Output file created** (output.h5 or similar)
- ☐ **Time steps progressing** (if running long enough)
- ☐ **Memory usage stable** (use `top` or `htop` to monitor)

## Performance Check

Optional: Verify simulation performance:
```bash
time ./HydroPlas --config config/main_test.yaml
# Note the execution time for future reference
```

## Documentation Reference

For more details, see:
- `README_FIX.md` - Quick start (recommended first read)
- `VISUAL_GUIDE.md` - Visual explanation of fixes
- `COMPLETE_FIX_SUMMARY.md` - Comprehensive technical details
- `CHANGES_DIFF.md` - Exact code differences

## Success Criteria

Your fix is successful if:
1. ✅ Code compiles without errors
2. ✅ Executable runs without segmentation fault
3. ✅ Initial conditions are correctly applied
4. ✅ Simulation progresses through time steps
5. ✅ Output files are created

## If Everything Works

Congratulations! The segmentation fault is fixed. You can now:
- Run your full simulation
- Adjust parameters in the config file
- Analyze results in the output HDF5 file

## If Problems Persist

1. Check that you're using the correct source file
2. Verify PETSc installation is working
3. Try compiling with debug flags: `cmake -DCMAKE_BUILD_TYPE=Debug ..`
4. Contact support with the specific error message

---

**Last Updated:** January 3, 2026  
**Fix Version:** 1.0  
**Status:** Ready for deployment
