# meep source map: Getting Started

Use this map only after reading `skills/meep-getting-started/references/doc_map.md`.

## Starter commands
- `python python/examples/straight-waveguide.py`
- `meep scheme/examples/straight-waveguide.ctl`
- `python python/examples/bend-flux.py`

## Validation checkpoints
- Baseline runs finish without instability/divergence.
- Fields propagate and are absorbed by PML (no strong late-time ringing).
- A basic observable (e.g., flux/transmission printout) remains finite and stable after a modest runtime increase.

## Function-level source entry points
- `python/simulation.py` - `Simulation.__init__`, `init_sim`, `run`, `run_k_point`, `run_k_points`.
- `python/simulation.py` - `stop_when_fields_decayed`, `at_every`, `at_beginning`, `at_end`.
- `python/source.py` - `ContinuousSource`, `GaussianSource`, `CustomSource`, `EigenModeSource`.
- `scheme/meep.scm.in` - `(run-until ...)`, `(run-sources ...)`, `(run-k-points ...)`, `(stop-when-fields-decayed ...)`.
- `src/step.cpp` - `fields::step`, `fields::step_source`, `fields::step_boundaries`.
- `src/sources.cpp` - `fields::add_point_source`, `gaussian_src_time`, `continuous_src_time`.
- `src/fields.cpp` - core field container setup and interpolation behavior.
- `src/meepgeom.cpp` - `set_materials_from_geometry` and geometry-to-grid mapping.

## Fast source navigation
- `rg -n "def (run|run_k_point|run_k_points|init_sim)|stop_when_fields_decayed" python/simulation.py`
- `rg -n "define \(run-until|define \(run-sources|define \(run-k-points|stop-when-fields-decayed" scheme/meep.scm.in`
- `rg -n "fields::step|add_point_source|set_materials_from_geometry" src`
