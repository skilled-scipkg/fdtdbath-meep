# meep source map: Build and Install

Use this map only after reading `skills/meep-build-and-install/references/doc_map.md`.

## Starter commands
- `conda create -n fdtdbath-run --override-channels -c tel-research -c conda-forge "pymeep-fdtdbath * mpi_mpich_*" && conda activate fdtdbath-run`
- `conda install --override-channels -c tel-research -c conda-forge "pymeep-fdtdbath * mpi_mpich_*"`
- `sh autogen.sh --enable-shared`
- `./configure --with-mpi --with-openmp PYTHON=python3`
- `make -j && make check`
- `mpirun -np 2 python -m mpi4py python/examples/straight-waveguide.py`

## Validation checkpoints
- `configure` summary reflects intended features (`MPI`, `OpenMP`, `HDF5`, Python/Scheme interfaces).
- `python -c 'import meep as mp; print(mp.__version__); print(hasattr(mp, "BathLorentzianSusceptibility"))'` succeeds in the target environment.
- At least one serial and one MPI example complete successfully.

## Function-level source entry points
- `configure.ac` - `AC_ARG_WITH(mpi)`, `AC_ARG_ENABLE(openmp)`, `AC_ARG_WITH(hdf5)`, `AC_ARG_WITH(python)`.
- `autogen.sh` - autotools bootstrap for git checkouts.
- `src/Makefile.am` - core `libmeep` build targets and linkage.
- `python/Makefile.am` - SWIG build rules (`meep-python.cxx`) and `_meep.la` packaging.
- `scheme/Makefile.am` - Scheme SWIG generation and interface build rules.
- `libpympb/Makefile.am` - MPB bridge build target.
- `src/support/Makefile.am` - support library linkage details.
- `python/meep.i` - Python binding interface definitions.
- `scheme/meep.i` - Scheme binding interface definitions.

## Fast source navigation
- `rg -n "AC_ARG_WITH\(mpi|AC_ARG_ENABLE\(openmp|AC_ARG_WITH\(hdf5|AC_ARG_WITH\(python" configure.ac`
- `rg -n "SWIG|_meep.la|libmeep_la|libpympb" src/Makefile.am python/Makefile.am scheme/Makefile.am libpympb/Makefile.am`
