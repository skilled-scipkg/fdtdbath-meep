# meep documentation map: Inputs and Modeling

Use this map from `skills/meep-inputs-and-modeling/SKILL.md`.

## Core modeling docs
- `doc/docs/Materials.md` - material models, conductivity, dispersion conventions.
- `doc/docs/Subpixel_Smoothing.md` - interface accuracy and convergence behavior.
- `doc/docs/Chunks_and_Symmetry.md` - chunking/symmetry implications.
- `doc/docs/Scheme_User_Interface.md` - full object/material/source schema.

## Modeling tutorials
- `doc/docs/Python_Tutorials/Material_Dispersion.md`
- `doc/docs/Scheme_Tutorials/Material_Dispersion.md`
- `doc/docs/Python_Tutorials/Third_Harmonic_Generation.md`
- `doc/docs/Scheme_Tutorials/Third_Harmonic_Generation.md`
- `doc/docs/Python_Tutorials/GDSII_Import.md`
- `doc/docs/Python_Tutorials/Mode_Decomposition.md`
- `doc/docs/Python_Tutorials/Multilevel_Atomic_Susceptibility.md`

## Runnable modeling examples
- `python/examples/material-dispersion.py`
- `python/examples/3rd-harm-1d.py`
- `python/examples/coupler.py`
- `scheme/examples/material-dispersion.ctl`
- `scheme/examples/3rd-harm-1d.ctl`

## Escalation
If modeling behavior still needs implementation-level checks, use `skills/meep-inputs-and-modeling/references/source_map.md`.
