# meep source map: Advanced Topics

Use this map only after reading `skills/meep-advanced-topics/references/doc_map.md`.

## Starter commands
- `python python/examples/refl-angular-kz2d.py`
- `python python/examples/oblique-source.py`
- `meep scheme/examples/refl-angular-kz2d.ctl`

## Validation checkpoints
- For nonzero `k_point`/`kz`, verify complex-field expectations and memory impact.
- Symmetry-enabled runs match non-symmetry baselines within discretization tolerance.
- PML thickness/padding sweeps reduce boundary-reflection artifacts.
- Courant changes preserve stability when using cylindrical or special-field options.

## Function-level source entry points
- `python/simulation.py` - `Simulation.__init__` options: `k_point`, `kz_2d`, `symmetries`, `force_complex_fields`, `Courant`.
- `python/simulation.py` - `run_k_point` and `run_k_points` for Bloch/eigenfrequency workflows.
- `scheme/meep.scm.in` - `define-param k-point`, `define-param symmetries`, `define-param Courant`, `define-param force-complex-fields?`.
- `scheme/meep.scm.in` - `(run-k-point ...)`, `(run-k-points ...)`, Bloch setup logic.
- `src/boundaries.cpp` - `fields::set_boundary`, periodic/Bloch boundary handling.
- `src/update_eh.cpp` - PML/boundary-relevant E/H updates.
- `src/update_pols.cpp` - polarization/nonlinearity update orchestration.
- `src/susceptibility.cpp` - dispersive/gyrotropic susceptibility updates.
- `src/fields.cpp` - Yee-grid component handling and interpolation behavior.
- `src/vec.cpp` - Yee shifts, component transforms, and symmetry transforms.
- `src/meepgeom.cpp` - nonlinear and material-grid handling (`chi2`, `chi3`, conductivity paths).
- `src/step.cpp` - timestep sequencing used by advanced runs.

## Fast source navigation
- `rg -n "k_point|kz_2d|symmetries|force_complex_fields|Courant|run_k_points" python/simulation.py`
- `rg -n "define-param (k-point|symmetries|Courant|force-complex-fields\?)|run-k-point|run-k-points" scheme/meep.scm.in`
- `rg -n "set_boundary|use_bloch|update_eh|update_pols|chi3|phase_shift|fields::step" src`
