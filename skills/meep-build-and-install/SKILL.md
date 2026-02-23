---
name: meep-build-and-install
description: Use this skill for Meep installation and build troubleshooting (Conda, source builds, MPI/OpenMP/HDF5, platform caveats, and post-install validation). Keep guidance docs-first and escalate to build scripts/source files only for unresolved details.
---

# meep: Build and Install

## High-Signal Playbook
### Route conditions
- Use `meep-getting-started` for first simulations after installation succeeds.
- Use `meep-simulation-workflows` for runtime convergence/output issues rather than install issues.
- Use `meep-api-and-scripting` for interface-level symbol or binding behavior questions.

### Triage questions
- Which platform (Linux/macOS/WSL/HPC cluster) and shell environment?
- Python-only workflow or Scheme interface also required?
- Serial run only or MPI parallel required?
- Conda installation acceptable, or source build mandatory?
- Is this from a release tarball or a `git clone` checkout?
- Are compiler/MPI/HDF5 stacks from one consistent toolchain?

### Canonical workflow
1. Prefer Conda (`pymeep`) for fastest reliable setup (`doc/docs/Installation.md`).
2. Choose serial (`pymeep`) or parallel (`pymeep=*=mpi_mpich_*`) environment.
3. Verify import and version, then run one known example script.
4. If source build is required, install prerequisites and pick one compiler toolchain consistently (`doc/docs/Build_From_Source.md`).
5. For `git clone`, run `autogen.sh`; for release tarball `meep-X.Y.Z.tar.gz`, configure directly.
6. Configure with required features (`--with-mpi`, `--with-openmp`, custom `PYTHON=` as needed), then `make`, `make check`, `make install`.
7. Validate serial and MPI execution paths before debugging performance.

### Minimal working example
Conda path (recommended in `doc/docs/Installation.md`):
```bash
conda create -n mp -c conda-forge pymeep pymeep-extras
conda activate mp
python -c 'import meep as mp; print(mp.__version__)'
python python/examples/straight-waveguide.py
```

Source path skeleton (from `doc/docs/Build_From_Source.md`):
```bash
# For git checkout; release tarballs usually skip this
sh autogen.sh --enable-shared

./configure --with-mpi --with-openmp PYTHON=python3
make -j
make check
sudo make install
```

### Pitfalls and fixes
- `illegal instruction` on import can be OpenBLAS mismatch; docs note downgrading OpenBLAS (`doc/docs/Installation.md`).
- MKL/OpenBLAS mixing in Conda can trigger segmentation faults; prefer strict `conda-forge` + `no-mkl` guidance (`doc/docs/Installation.md`).
- Using GitHub auto-generated `vX.Y.Z.tar.gz` instead of `meep-X.Y.Z.tar.gz` changes build requirements (`doc/docs/Build_From_Source.md`).
- MPI + HDF5 toolchain mismatches cause link/runtime failures; compile/link with consistent MPI compilers (`doc/docs/Build_From_Source.md`).
- Supercomputer installs often fail from mixed compiler families; keep dependencies and Meep on one vendor/compiler stack.
- Parallel Python job hangs after one process fails; run via `python -m mpi4py` per parallel docs (`doc/docs/Parallel_Meep.md`).

### Convergence and validation checks
- Verify Python module load and version: `python -c 'import meep as mp; print(mp.__version__)'`.
- Verify Scheme binary if built: `meep` with a minimal `.ctl` file.
- Verify MPI launch path: `mpirun -np 4 python -m mpi4py <script>.py`.
- Run `make check` (or `make RUNCODE="env OMP_NUM_THREADS=4 mpirun -np 3" check`) for source builds.
- Confirm that one known tutorial script runs end-to-end before custom scripts.

## Scope
- Handle build, installation, compilation, and environment setup for Meep.
- Keep responses docs-first and reproducible.

## Primary documentation references
- `doc/docs/Installation.md`
- `doc/docs/Build_From_Source.md`
- `doc/docs/Parallel_Meep.md`
- `doc/docs/Download.md`
- `doc/docs/FAQ.md`
- `doc/README.md`

## Workflow
- Start from installation/build docs first.
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
