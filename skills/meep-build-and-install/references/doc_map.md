# meep documentation map: Build and Install

Use this map from `skills/meep-build-and-install/SKILL.md`.

## Primary install/build docs
- `doc/docs/Installation.md` - recommended Conda paths and troubleshooting.
- `doc/docs/Build_From_Source.md` - source-build flow and dependency expectations.
- `doc/docs/Parallel_Meep.md` - MPI launch patterns and scaling caveats.
- `doc/docs/Download.md` - release artifact selection (`meep-X.Y.Z.tar.gz` vs source snapshots).
- `doc/docs/FAQ.md` - common failure modes (convergence/runtime/import).
- `doc/README.md` - documentation build commands useful when local docs need rebuilding.

## Quick validation references
- `python/examples/straight-waveguide.py`
- `python/examples/bend-flux.py`
- `scheme/examples/straight-waveguide.ctl`

## Build-system files to cross-check with docs
- `autogen.sh`
- `configure.ac`
- `src/Makefile.am`
- `python/Makefile.am`
- `scheme/Makefile.am`
- `libpympb/Makefile.am`

## Escalation
If docs are insufficient, use `skills/meep-build-and-install/references/source_map.md`.
