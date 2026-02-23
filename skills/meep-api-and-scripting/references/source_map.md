# meep source map: API and Scripting

Use this map only after reading `skills/meep-api-and-scripting/references/doc_map.md`.

## Starter commands
- `python python/examples/straight-waveguide.py`
- `meep scheme/examples/straight-waveguide.ctl`
- `python python/examples/solve-cw.py`
- `python python/examples/mode-decomposition.py`

## Validation checkpoints
- API calls used in scripts map cleanly to documented interfaces.
- `run` step functions execute in the expected order.
- Flux/mode outputs are finite and reproducible under a small runtime increase.

## Function-level source entry points
- `python/simulation.py` - `Simulation.run`, `add_flux`, `add_near2far`, `solve_cw`, `get_eigenmode_coefficients`.
- `python/simulation.py` - `stop_when_fields_decayed`, `at_every`, `at_beginning`, `synchronized_magnetic`.
- `python/source.py` - source object behavior (`ContinuousSource`, `GaussianSource`, `EigenModeSource`).
- `python/meep.i` - Python SWIG symbol exposure and typemap boundaries.
- `python/meep-python.hpp` - Python/C++ bridge declarations.
- `scheme/meep.scm.in` - `(run-until ...)`, `(run-sources ...)`, `(run-k-points ...)`, `(at-every ...)`.
- `scheme/meep.i` - Scheme SWIG interface.
- `scheme/meep.cpp` - Scheme/C++ glue implementation.
- `scheme/meep_op_renames.i` - compatibility renames for Scheme symbols.
- `src/dft.cpp` - `fields::add_dft_flux`, `fields::get_dft_array`, mode-monitor plumbing.
- `src/near2far.cpp` - `fields::add_dft_near2far`, `dft_near2far::farfield`.
- `src/step.cpp` - runtime stepping behavior underlying high-level APIs.

## Fast source navigation
- `rg -n "def (run|add_flux|add_near2far|solve_cw|get_eigenmode_coefficients)|def stop_when_fields_decayed" python/simulation.py`
- `rg -n "define \(run-until|define \(run-sources|define \(run-k-points|define \(at-every" scheme/meep.scm.in`
- `rg -n "add_dft_flux|add_dft_near2far|get_dft_array|farfield" src/dft.cpp src/near2far.cpp`
