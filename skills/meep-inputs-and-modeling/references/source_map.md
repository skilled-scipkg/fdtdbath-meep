# meep source map: Inputs and Modeling

Use this map only after reading `skills/meep-inputs-and-modeling/references/doc_map.md`.

## Starter commands
- `python python/examples/material-dispersion.py`
- `python python/examples/3rd-harm-1d.py`
- `python python/examples/coupler.py`
- `meep scheme/examples/material-dispersion.ctl`

## Validation checkpoints
- Material/geometry edits preserve numerical stability (no late-time blow-up).
- Target metrics converge under resolution and PML thickness sweeps.
- For dispersive/nonlinear cases, parameter changes produce smooth, physically plausible trends.

## Function-level source entry points
- `python/simulation.py` - constructor knobs for `geometry`, `boundary_layers`, `symmetries`, `k_point`, `Courant`.
- `python/geom.py` - geometric object and vector abstractions used by Python input files.
- `python/materials.py` - built-in material definitions and dispersive parameter sets.
- `scheme/materials.scm` - Scheme-side material definitions.
- `src/material_data.hpp` - material data structures and susceptibility interfaces.
- `src/material_data.cpp` - material interpolation/storage implementation.
- `src/meepgeom.hpp` - geometry/material API declarations.
- `src/meepgeom.cpp` - `set_materials_from_geometry`, material-grid interpolation, `chi2`/`chi3` paths.
- `src/structure.cpp` - structure assembly and chunk/material assignment.
- `src/fix_boundary_sources.cpp` - source correction near boundaries.
- `src/GDSIIgeom.cpp` - `set_geometry_from_GDSII`, `get_GDSII_prisms`, layer extraction.
- `src/susceptibility.cpp` - Lorentzian/gyrotropic update paths (`update_P`, `chi1`).
- `src/multilevel-atom.cpp` - multilevel susceptibility updates.
- `src/update_pols.cpp` - polarization update orchestration.

## Fast source navigation
- `rg -n "k_point|symmetries|Courant|boundary_layers|geometry" python/simulation.py`
- `rg -n "set_materials_from_geometry|material_grid_val|chi2|chi3|conductivity" src/meepgeom.cpp`
- `rg -n "update_P|chi1|multilevel_susceptibility|update_pols" src/susceptibility.cpp src/multilevel-atom.cpp src/update_pols.cpp`
- `rg -n "GDSII|set_geometry_from_GDSII|get_GDSII_prisms" src/GDSIIgeom.cpp`
