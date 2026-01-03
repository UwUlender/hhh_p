# PETSc Configuration Fix - Documentation Index

## 🎯 Start Here

If you're seeing PETSc errors about "GMRES" or "PBP", start with:

1. **[FIX_COMPLETE_README.md](FIX_COMPLETE_README.md)** - Complete overview and what to do next
2. **[QUICK_FIX_REFERENCE.txt](QUICK_FIX_REFERENCE.txt)** - Quick reference card
3. **[CHANGE_DIAGRAM.txt](CHANGE_DIAGRAM.txt)** - Visual diagram of changes

---

## 📚 Documentation Structure

### Quick Start Documents

| Document | Purpose | When to Read |
|----------|---------|--------------|
| **FIX_COMPLETE_README.md** | Complete overview with all info | Start here! |
| **QUICK_FIX_REFERENCE.txt** | One-page reference card | Quick lookup |
| **CHANGE_DIAGRAM.txt** | Visual flow diagrams | Understanding changes |

### Detailed Documentation

| Document | Purpose | When to Read |
|----------|---------|--------------|
| **PETSC_SOLVER_FIX.md** | Detailed fix explanation | Understanding the problem |
| **docs/PETSC_SOLVER_GUIDE.md** | Complete solver configuration guide | Configuring solvers |
| **PETSC_FIX_SUMMARY.md** | Summary of all changes | Review what was done |
| **GIT_COMMIT_INFO.md** | Git commit template and diffs | Committing changes |

### Updated Main Documentation

| Document | Purpose | What's New |
|----------|---------|-----------|
| **README.md** | Main project README | Added troubleshooting section |
| **docs/USER_GUIDE.md** | User guide | Reference to solver guide |

---

## 🛠️ Tools and Scripts

### Rebuild Script
```bash
./rebuild_fix.sh
```
Automated rebuild script that cleans, configures, and builds the project.

### Configuration Validator
```bash
python3 validate_config.py
```
Validates all YAML configuration files to ensure correct PETSc types.

---

## 📁 File Changes Summary

### Modified Files (5)

| File | What Changed |
|------|--------------|
| `src/config/ConfigParser.hpp` | Default solver types fixed (lines 85-86) |
| `src/config/ConfigParser.cpp` | Default solver types fixed (lines 217-221, 222-230) |
| `config/default_config.yaml` | Added solver section with correct types |
| `config/complete_feature_demo.yaml` | Fixed solver types (lines 147-148) |
| `README.md` | Added PETSc troubleshooting section |

### New Files (12)

| File | Purpose |
|------|---------|
| `PETSC_SOLVER_FIX.md` | Detailed fix explanation |
| `docs/PETSC_SOLVER_GUIDE.md` | Complete solver configuration guide |
| `PETSC_FIX_SUMMARY.md` | Summary of all changes |
| `QUICK_FIX_REFERENCE.txt` | Quick reference card |
| `GIT_COMMIT_INFO.md` | Git commit template |
| `FIX_COMPLETE_README.md` | Comprehensive overview |
| `CHANGE_DIAGRAM.txt` | Visual diagrams |
| `INDEX.md` | This file - documentation index |
| `config/test_petsc_fix.yaml` | Test configuration file |
| `rebuild_fix.sh` | Automated rebuild script |
| `validate_config.py` | Configuration validator |

---

## 🚀 Quick Action Plan

### Step 1: Understand the Problem
Read: **FIX_COMPLETE_README.md** (5 minutes)

### Step 2: Rebuild
Run: `./rebuild_fix.sh` (2-5 minutes)

### Step 3: Validate
Run: `python3 validate_config.py` (10 seconds)

### Step 4: Test
Run: `./HydroPlas --config config/test_petsc_fix.yaml` (1-2 minutes)

### Step 5: Verify
Check: No PETSc errors, simulation runs successfully

---

## 📖 Reading Guide by Use Case

### "I just want to fix the error quickly"
1. Read: **QUICK_FIX_REFERENCE.txt**
2. Run: `./rebuild_fix.sh`
3. Test: `./HydroPlas --config config/test_petsc_fix.yaml`

### "I want to understand what went wrong"
1. Read: **FIX_COMPLETE_README.md** (Executive Summary)
2. Read: **PETSC_SOLVER_FIX.md** (Detailed explanation)
3. View: **CHANGE_DIAGRAM.txt** (Visual diagrams)

### "I need to configure solvers for my simulation"
1. Read: **docs/PETSC_SOLVER_GUIDE.md** (Complete guide)
2. Use: **QUICK_FIX_REFERENCE.txt** (Valid types reference)
3. Validate: `python3 validate_config.py`

### "I want to commit these changes to git"
1. Read: **GIT_COMMIT_INFO.md** (Commit template)
2. Review: **PETSC_FIX_SUMMARY.md** (Summary of changes)
3. Commit with provided message template

### "I'm troubleshooting solver convergence"
1. Read: **docs/PETSC_SOLVER_GUIDE.md** (Performance section)
2. Read: **FIX_COMPLETE_README.md** (Troubleshooting section)
3. Try: Different solver configurations from examples

---

## 🔍 Quick Reference Tables

### Valid KSP Types (Must be lowercase!)
```
gmres, fgmres, bcgs, cg, bicg, tfqmr, richardson
```

### Valid PC Types (Parallel-safe)
```
pbjacobi, bjacobi, asm, jacobi, none
```

### Valid PC Types (Serial only)
```
ilu, sor, lu, cholesky, icc
```

### Invalid Types to Avoid
```
❌ GMRES   → Use: gmres
❌ PBP     → Use: pbjacobi
❌ ILU     → Use: ilu
```

---

## 🎓 Learning Path

### Beginner: Just Fix the Error
1. **FIX_COMPLETE_README.md** - Overview
2. `./rebuild_fix.sh` - Rebuild
3. Test and verify

### Intermediate: Understand and Configure
1. **PETSC_SOLVER_FIX.md** - Problem explanation
2. **docs/PETSC_SOLVER_GUIDE.md** - Configuration guide
3. **QUICK_FIX_REFERENCE.txt** - Quick reference
4. Experiment with different solver types

### Advanced: Optimize Performance
1. **docs/PETSC_SOLVER_GUIDE.md** - Performance tips
2. **FIX_COMPLETE_README.md** - Troubleshooting
3. Use PETSc monitoring: `-ksp_monitor -snes_monitor`
4. Profile and optimize for your specific problem

---

## 📊 Validation Status

All configuration files have been validated:

```
✅ config/complete_feature_demo.yaml
✅ config/default_config.yaml
✅ config/test_petsc_fix.yaml
```

Run `python3 validate_config.py` to verify your own configurations.

---

## 🎯 Key Takeaways

1. **PETSc type names are case-sensitive** - Always use lowercase
2. **"PBP" is not valid** - Use `pbjacobi` or other valid types
3. **Parallel vs Serial** - Some preconditioners only work in serial
4. **Default is safe** - `gmres` + `pbjacobi` works for most cases
5. **Validation is easy** - Use `validate_config.py` to check configs

---

## 🆘 Getting Help

### If validation fails:
→ Check **QUICK_FIX_REFERENCE.txt** for valid types

### If build fails:
→ Check **FIX_COMPLETE_README.md** troubleshooting section

### If solver doesn't converge:
→ Read **docs/PETSC_SOLVER_GUIDE.md** performance tips

### If PETSc errors persist:
→ Verify types are lowercase and valid
→ Try runtime override: `-ksp_type gmres -pc_type pbjacobi`

---

## 📝 Checklist

Before considering the fix complete:

- [ ] Read **FIX_COMPLETE_README.md**
- [ ] Run `./rebuild_fix.sh`
- [ ] Run `python3 validate_config.py`
- [ ] Test with `config/test_petsc_fix.yaml`
- [ ] Verify no PETSc errors
- [ ] Simulation runs successfully
- [ ] (Optional) Commit changes using **GIT_COMMIT_INFO.md**

---

## 📞 Status

**Status:** ✅ **COMPLETE**

All files have been updated, validated, and documented.
Ready for rebuild and testing.

**Last Updated:** January 3, 2026

---

## 🔗 External Resources

- [PETSc Documentation](https://petsc.org/release/docs/)
- [PETSc KSP Manual](https://petsc.org/release/docs/manual/ksp/)
- [PETSc PC Manual](https://petsc.org/release/docs/manual/pc/)
- [PETSc FAQ](https://petsc.org/release/faq/)

---

**Navigation:**
- 📖 [Complete Overview](FIX_COMPLETE_README.md)
- 🚀 [Quick Reference](QUICK_FIX_REFERENCE.txt)
- 🔧 [Detailed Fix](PETSC_SOLVER_FIX.md)
- 📚 [Solver Guide](docs/PETSC_SOLVER_GUIDE.md)
- 📊 [Change Summary](PETSC_FIX_SUMMARY.md)
