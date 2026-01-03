# HydroPlas Improvement Checklist

**Based on Zapdos Comparison Analysis**  
**Started:** [DATE]  
**Target Completion:** [DATE + 3 months]

---

## 🔴 Week 1: Foundation (CRITICAL)

### Day 1-2: Docker Deployment
- [ ] Create `Dockerfile` in root directory
- [ ] Create `docker-entrypoint.sh` script
- [ ] Test build locally: `docker build -t hydroplas:test .`
- [ ] Test run: `docker run hydroplas:test /config/argon_complete.json`
- [ ] Create Docker Hub account/organization
- [ ] Push image: `docker push hydroplasma/hydroplas:latest`
- [ ] Update README with Docker instructions
- [ ] Test on clean machine (verify dependencies work)

### Day 3: GitHub Setup
- [ ] Enable GitHub Issues (Settings → Features → Issues)
- [ ] Enable GitHub Discussions (Settings → Features → Discussions)
- [ ] Create issue templates:
  - [ ] Bug report template
  - [ ] Feature request template
  - [ ] Question template
- [ ] Create discussion categories:
  - [ ] General
  - [ ] Help/Q&A
  - [ ] Show and Tell
  - [ ] Ideas
- [ ] Create `CONTRIBUTING.md` file
- [ ] Create `CODE_OF_CONDUCT.md` file

### Day 4: CI/CD
- [ ] Create `.github/workflows/ci.yml`
- [ ] Configure build job (Ubuntu, macOS)
- [ ] Test CI by pushing a commit
- [ ] Fix any build errors that appear
- [ ] Add status badge to README
- [ ] Set up automatic Docker builds on push

### Day 5: Academic Credibility
- [ ] Create `CITATION.cff` file with metadata
- [ ] Register project with Zenodo
- [ ] Link Zenodo to GitHub repository
- [ ] Create first release: `git tag v1.0.0`
- [ ] Get DOI from Zenodo
- [ ] Add DOI badge to README
- [ ] Create `AUTHORS.md` file

**Week 1 Success Criteria:**
✅ Docker image available on Docker Hub  
✅ GitHub Actions CI running  
✅ Issues/Discussions enabled  
✅ DOI obtained from Zenodo

---

## 🔴 Week 2: Testing Infrastructure (CRITICAL)

### Unit Testing Setup
- [ ] Add Catch2 to `CMakeLists.txt`
- [ ] Create `test/` directory structure
- [ ] Create `test/CMakeLists.txt`
- [ ] Write first test: `test/unit/test_main.cpp` (sanity check)
- [ ] Verify tests run: `cd build/test && ctest`

### Unit Tests - Chemistry
- [ ] `test_chemistry_ionization_rate.cpp`
- [ ] `test_chemistry_excitation_rate.cpp`
- [ ] `test_chemistry_conservation.cpp`
- [ ] `test_chemistry_table_lookup.cpp`

### Unit Tests - Numerics
- [ ] `test_scharfetter_gummel_drift.cpp`
- [ ] `test_scharfetter_gummel_advection.cpp`
- [ ] `test_flux_scheme_conservation.cpp`

### Unit Tests - Boundary
- [ ] `test_boundary_voltage_dc.cpp`
- [ ] `test_boundary_voltage_rf.cpp`
- [ ] `test_boundary_dielectric.cpp`

### Regression Testing
- [ ] Create `test/regression/` directory
- [ ] Copy 3 config files to `test/regression/configs/`
- [ ] Run and save as gold files: `test/regression/gold/`
- [ ] Write `test/regression/run_regression.py`
- [ ] Create `test/regression/test_cases.yaml`
- [ ] Add regression tests to CI

### Coverage
- [ ] Add `--coverage` flag to CMake
- [ ] Generate coverage report
- [ ] Upload to Codecov (optional)
- [ ] Add coverage badge to README

**Week 2 Success Criteria:**
✅ 15+ unit tests passing  
✅ 3+ regression tests with gold files  
✅ Tests run automatically in CI  
✅ Coverage > 50%

---

## 🟡 Week 3: Documentation Website (HIGH)

### MkDocs Setup
- [ ] Install mkdocs: `pip install mkdocs-material`
- [ ] Create `mkdocs.yml` configuration
- [ ] Create `docs/` directory structure
- [ ] Write `docs/index.md` (landing page)
- [ ] Add mathematical notation support (MathJax)
- [ ] Test locally: `mkdocs serve`
- [ ] Configure GitHub Pages in repo settings

### Core Documentation
- [ ] `docs/getting-started/installation.md`
- [ ] `docs/getting-started/quick-start.md`
- [ ] `docs/getting-started/docker.md`
- [ ] `docs/theory/equations.md` (copy/adapt from THEORY.md)
- [ ] `docs/theory/numerics.md`
- [ ] `docs/theory/chemistry.md`
- [ ] `docs/reference/configuration.md`
- [ ] `docs/reference/output_format.md`

### Jupyter Tutorials
- [ ] Create `docs/tutorials/` directory
- [ ] Tutorial 1: `01_basic_discharge.ipynb`
  - [ ] Setup and run
  - [ ] Visualize results
  - [ ] Parameter study
- [ ] Tutorial 2: `02_excited_species.ipynb`
  - [ ] Metastable transport
  - [ ] Stepwise ionization
  - [ ] Penning effect
- [ ] Tutorial 3: `03_multi_electrode.ipynb`
  - [ ] Dual-frequency setup
  - [ ] Phase control
  - [ ] Comparison with single electrode
- [ ] Convert notebooks to markdown for web

### API Documentation
- [ ] Create `Doxyfile` configuration
- [ ] Add Doxygen comments to 5 key classes
- [ ] Generate HTML docs: `doxygen Doxyfile`
- [ ] Link API docs to MkDocs site
- [ ] Add "API Reference" section to docs

### Deploy
- [ ] Create `.github/workflows/docs.yml`
- [ ] Push to trigger documentation build
- [ ] Verify site live at https://yourorg.github.io/hydroplas
- [ ] Add website link to GitHub repo description

**Week 3 Success Criteria:**
✅ Professional documentation website live  
✅ 3 interactive Jupyter tutorials  
✅ Getting started guide complete  
✅ API documentation generated

---

## 🟡 Week 4: Extensibility (HIGH)

### Plugin Architecture Design
- [ ] Create `src/kernels/PhysicsKernel.hpp` base class
- [ ] Create `src/kernels/KernelManager.hpp`
- [ ] Implement `KernelManager.cpp`
- [ ] Define `SolverState` struct

### Refactor Existing Code
- [ ] Extract ionization to `IonizationKernel.cpp`
- [ ] Extract diffusion to `DiffusionKernel.cpp`
- [ ] Extract drift to `DriftKernel.cpp`
- [ ] Modify `PlasmaSolver` to use `KernelManager`
- [ ] Test that results unchanged (regression tests!)

### Example Plugins
- [ ] Create `examples/plugins/` directory
- [ ] Write `CustomHeatingKernel.cpp` example
- [ ] Write `UserDefinedReactionKernel.cpp` example
- [ ] Create Makefile for plugins
- [ ] Test dynamic loading

### Configuration Support
- [ ] Extend JSON schema to support kernel config
- [ ] Parse kernel list from config file
- [ ] Load plugins from `library` path
- [ ] Pass parameters to kernel `initialize()`

### Documentation
- [ ] Write `docs/development/plugins.md` guide
- [ ] Document `PhysicsKernel` API
- [ ] Create plugin template repository
- [ ] Add examples to documentation

**Week 4 Success Criteria:**
✅ Plugin system working  
✅ 2 example plugins created  
✅ Existing code refactored without breakage  
✅ Plugin development guide written

---

## 🟢 Month 2: Enhancements (MEDIUM)

### Adaptive Time Stepping
- [ ] Use `TSAdapt` from PETSc
- [ ] Add to config: `"adaptive": true, "rtol": 1e-4`
- [ ] Test on stiff problem
- [ ] Document in user guide

### Boundary Condition Library
- [ ] Create `src/bcs/BoundaryCondition.hpp` base class
- [ ] Implement `HagelaarBC.cpp`
- [ ] Implement `SchottkyEmissionBC.cpp`
- [ ] Add to config: per-boundary BC specification
- [ ] Test with example cases
- [ ] Document new BCs

### Mesh Import
- [ ] Write `src/mesh/GmshReader.cpp`
- [ ] Parse Gmsh ASCII format
- [ ] Convert to PETSc DMDA (or DMPlex)
- [ ] Test with simple 2D geometry
- [ ] Document mesh import workflow

### Enhanced Output
- [ ] Create `src/postprocessors/` directory
- [ ] Implement derived quantity framework
- [ ] Add electric field calculation
- [ ] Add power density calculation
- [ ] Support formulas in config
- [ ] Test and document

### Code Quality
- [ ] Add `.clang-format` file
- [ ] Format all code: `find src -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i`
- [ ] Add clang-tidy checks
- [ ] Fix warnings: compile with `-Wall -Wextra -Werror`
- [ ] Run Valgrind memory check
- [ ] Fix any leaks/errors found

**Month 2 Success Criteria:**
✅ Adaptive timestepping working  
✅ 5+ BC types available  
✅ Mesh import from Gmsh functional  
✅ Code formatted and warning-free

---

## 🟢 Month 3: Community & Publication (MEDIUM)

### Example Library
- [ ] Create 10 well-documented example configs
- [ ] DBD with metastables
- [ ] Plasma jet
- [ ] RF discharge (single & dual frequency)
- [ ] Penning mixture
- [ ] Each with README explaining physics
- [ ] Add to documentation website

### Video Tutorial
- [ ] Record 10-15 minute walkthrough
- [ ] Cover: install, run, visualize
- [ ] Upload to YouTube
- [ ] Add to documentation
- [ ] Post on social media (Twitter/LinkedIn)

### Methods Paper
- [ ] Draft paper outline
- [ ] Write Introduction (motivation, excited species)
- [ ] Write Methods (numerics, architecture)
- [ ] Create 3-5 figures (code structure, validation)
- [ ] Run validation cases (comparison with literature)
- [ ] Write Results section
- [ ] Submit to JOSS or Computer Physics Communications

### Community Building
- [ ] Post announcement on relevant forums
- [ ] Create Twitter/X account for updates
- [ ] Email potential users (plasma research groups)
- [ ] Respond to all GitHub issues/discussions
- [ ] Recruit 3-5 beta testers
- [ ] Incorporate feedback

### Performance Benchmarking
- [ ] Reproduce 3 Zapdos test cases in HydroPlas
- [ ] Compare time-to-solution
- [ ] Compare memory usage
- [ ] Document scaling (1-16 cores)
- [ ] Create benchmarking report
- [ ] Add to documentation

**Month 3 Success Criteria:**
✅ Methods paper submitted  
✅ 10+ example cases documented  
✅ Video tutorial published  
✅ 5+ external users engaged  
✅ Performance benchmarks completed

---

## 🔵 Future Work (LOW PRIORITY)

### Advanced Features
- [ ] DMPlex integration (unstructured meshes)
- [ ] GPU acceleration exploration
- [ ] Multi-physics coupling (thermal, flow)
- [ ] Machine learning for chemistry acceleration
- [ ] Graphical interface (web-based)

### Long-term Community
- [ ] Annual workshop
- [ ] Slack/Discord server
- [ ] Mailing list
- [ ] Conference presentations
- [ ] Collaboration with other plasma codes

---

## Progress Tracking

### Weekly Review
- [ ] Week 1 complete: [DATE]
- [ ] Week 2 complete: [DATE]
- [ ] Week 3 complete: [DATE]
- [ ] Week 4 complete: [DATE]
- [ ] Month 2 complete: [DATE]
- [ ] Month 3 complete: [DATE]

### Metrics Dashboard

| Metric | Target | Current | Status |
|--------|--------|---------|--------|
| Unit Tests | 30+ | __ | ⏳ |
| Regression Tests | 10+ | __ | ⏳ |
| Code Coverage | 70% | __% | ⏳ |
| GitHub Stars | 50+ | __ | ⏳ |
| Docker Pulls | 100+ | __ | ⏳ |
| Documentation Pages | 20+ | __ | ⏳ |
| External Users | 5+ | __ | ⏳ |
| Papers Using HydroPlas | 1+ | __ | ⏳ |

### Blockers / Issues

1. Issue: ___________  
   Status: ⏳ / ✅ / ❌  
   Resolution: ___________

2. Issue: ___________  
   Status: ⏳ / ✅ / ❌  
   Resolution: ___________

---

## Notes and Decisions

### Decision Log

| Date | Decision | Rationale |
|------|----------|-----------|
| [DATE] | Example decision | Example rationale |
| | | |
| | | |

### Deferred Items

Items decided NOT to implement (with reasons):

1. **GUI (Peacock-like):** High effort, low ROI for target users
2. **Exodus output:** HDF5 sufficient, no added value
3. **Thread parallelism:** MPI adequate for current use cases

---

## Resources

### Development Tools
- [Docker Desktop](https://www.docker.com/products/docker-desktop)
- [Catch2](https://github.com/catchorg/Catch2)
- [MkDocs Material](https://squidfunk.github.io/mkdocs-material/)
- [Doxygen](https://www.doxygen.nl/)
- [Zenodo](https://zenodo.org/)

### Learning Resources
- [Zapdos Documentation](https://shannon-lab.github.io/zapdos)
- [PETSc Manual](https://petsc.org/release/docs/manual/)
- [MOOSE Framework](https://mooseframework.inl.gov)
- [Software Carpentry](https://software-carpentry.org/)

### Communication Channels
- GitHub Issues: [link]
- GitHub Discussions: [link]
- Email: [email]
- Slack/Discord: [link if created]

---

## Completion Certificate

When all critical and high-priority items are complete:

```
✅ HydroPlas Modernization Complete!

Achievements:
- Professional documentation website
- Comprehensive testing (70%+ coverage)
- Docker deployment
- Plugin architecture
- Methods paper submitted
- Active community (5+ users)

Date: [DATE]
Team: [NAMES]

Next Phase: [What's next?]
```

---

**Remember:** Progress over perfection. Complete Week 1 before moving to Week 2!

**Last Updated:** [DATE]  
**Maintained by:** [NAME]
