---
name: meep-build-and-install
description: Use this skill for installation and build troubleshooting in this modified Meep repository (includes Lorentz-Bath/FDTD-Bath support for condensed-phase polaritons). Prefer the `pymeep-fdtdbath` Conda package from `tel-research`; use source builds only when package installation cannot satisfy requirements.
---

# meep: Build and Install

Repository note: this repository is a modified Meep distribution with Lorentz-Bath/FDTD-Bath extensions for condensed-phase polaritons.

## High-Signal Playbook
### Route conditions
- Use `meep-getting-started` for first simulations after installation succeeds.
- Use `meep-simulation-workflows` for runtime convergence/output issues rather than install issues.
- Use `meep-api-and-scripting` for interface-level symbol or binding behavior questions.

### Triage questions
- Which platform (Linux/macOS/WSL/HPC cluster) and shell environment?
- Python-only workflow or Scheme interface also required?
- Serial run only or MPI parallel required?
- Can `tel-research` + `conda-forge` channels be used, or is source build mandatory?
- Is this a fresh environment (`conda create`) or an existing one (`conda install`)?
- Is this from a release tarball or a `git clone` checkout?
- Are compiler/MPI/HDF5 stacks from one consistent toolchain?

### Canonical workflow
1. Prefer Conda package install for this modified build: `pymeep-fdtdbath` from `tel-research` + `conda-forge`.
2. Use either a fresh environment (`conda create`) or existing environment (`conda install`) with `--override-channels` and selector `"pymeep-fdtdbath * mpi_mpich_*"`.
3. Activate the environment, verify import/version, and confirm Lorentz-Bath API visibility.
4. Run one known example script before any custom model.
5. If source build is required, install prerequisites and keep compiler/MPI/HDF5 toolchains consistent (`doc/docs/Build_From_Source.md`).
6. For `git clone`, run `autogen.sh`; for release tarball `meep-X.Y.Z.tar.gz`, configure directly.
7. Configure required features (`--with-mpi`, `--with-openmp`, custom `PYTHON=` as needed), then `make`, `make check`, `make install`.
8. Validate serial and MPI execution paths before debugging performance.

### Minimal working example
Conda path (recommended for this modified package):
```bash
# Option A: create a new environment
conda create -n fdtdbath-run --override-channels \
  -c tel-research -c conda-forge \
  "pymeep-fdtdbath * mpi_mpich_*"
conda activate fdtdbath-run

# Option B: install into an existing environment
conda install --override-channels \
  -c tel-research -c conda-forge \
  "pymeep-fdtdbath * mpi_mpich_*"

python -c 'import meep as mp; print(mp.__version__); print(hasattr(mp, "BathLorentzianSusceptibility"))'
python python/examples/straight-waveguide.py
```

Source path skeleton (fallback from `doc/docs/Build_From_Source.md`):
```bash
# For git checkout; release tarballs usually skip this
sh autogen.sh --enable-shared

./configure --with-mpi --with-openmp PYTHON=python3
make -j
make check
sudo make install
```

### Pitfalls and fixes
- Missing `--override-channels` or incorrect channel priority can pull upstream `pymeep` (without this modified Lorentz-Bath extension).
- Package name is `pymeep-fdtdbath` but import stays `import meep as mp`.
- Using GitHub auto-generated `vX.Y.Z.tar.gz` instead of `meep-X.Y.Z.tar.gz` changes build requirements (`doc/docs/Build_From_Source.md`).
- MPI + HDF5 toolchain mismatches cause link/runtime failures; compile/link with consistent MPI compilers (`doc/docs/Build_From_Source.md`).
- Supercomputer installs often fail from mixed compiler families; keep dependencies and Meep on one vendor/compiler stack.
- Parallel Python job hangs after one process fails; run via `python -m mpi4py` per parallel docs (`doc/docs/Parallel_Meep.md`).

### Convergence and validation checks
- Verify Python module load/version and Lorentz-Bath symbol: `python -c 'import meep as mp; print(mp.__version__); print(hasattr(mp, "BathLorentzianSusceptibility"))'`.
- Verify Scheme binary if built: `meep` with a minimal `.ctl` file.
- Verify MPI launch path: `mpirun -np 4 python -m mpi4py <script>.py`.
- Run `make check` (or `make RUNCODE="env OMP_NUM_THREADS=4 mpirun -np 3" check`) for source builds.
- Confirm that one known tutorial script runs end-to-end before custom scripts.

## Scope
- Handle build, installation, compilation, and environment setup for this modified Meep package.
- Keep responses docs-first and reproducible.

## Primary documentation references
- `doc/docs/Installation.md`
- `doc/docs/Build_From_Source.md`
- `doc/docs/Parallel_Meep.md`
- `doc/docs/Download.md`
- `doc/docs/FAQ.md`
- `doc/README.md`

## Workflow
- Start from the `pymeep-fdtdbath` Conda commands in this skill.
- Use installation/build docs for troubleshooting and source-build fallback.
- If details are missing, inspect `skills/meep-build-and-install/references/doc_map.md`.
- Escalate to `skills/meep-build-and-install/references/source_map.md` only for unresolved build-system details.
- Cite exact file paths in responses.

## Tutorials and examples
- `python/examples`
- `scheme/examples`
- `doc/docs/Python_Tutorials`
- `doc/docs/Scheme_Tutorials`

## Test references
- `tests`
- `python/tests`

## Optional deeper inspection
- `libpympb`
- `python`
- `scheme`
- `src`

## Source entry points for unresolved issues (behavior-level)
- `configure.ac` - configure flags, dependency detection, MPI/OpenMP/HDF5/libctl checks.
- `src/Makefile.am` - C++ build targets and linkage.
- `python/Makefile.am` - Python module build/install rules.
- `scheme/Makefile.am` - Scheme interface build/install rules.
- `libpympb/Makefile.am` - MPB Python bridge build path.
- `src/support/Makefile.am` - support library build details.
- `python/meep.i` - SWIG Python binding surface.
- Prefer targeted search first: `rg -n "<symbol_or_keyword>" configure.ac python scheme src libpympb`.
