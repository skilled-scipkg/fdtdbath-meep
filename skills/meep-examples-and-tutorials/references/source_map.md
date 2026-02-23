# meep source map: Examples and Tutorials

Use this map only after reading `skills/meep-examples-and-tutorials/references/doc_map.md`.

## Starter commands
- `python python/examples/straight-waveguide.py`
- `python python/examples/bend-flux.py`
- `python python/examples/antenna-radiation.py`
- `meep scheme/examples/straight-waveguide.ctl`
- `meep scheme/examples/bend-flux.ctl`

## Validation checkpoints
- The unmodified baseline example runs before any adaptation.
- Adaptations are made one parameter block at a time and compared to baseline outputs.
- Metrics (flux/mode/far-field) remain physically consistent after each modification.

## Function-level source entry points
- `python/examples/README.md` - canonical Python example index and grouping.
- `python/examples` - runnable templates by problem type.
- `scheme/examples` - Scheme control-file analogs.
- `python/simulation.py` - `run`, `add_flux`, `add_near2far`, `get_eigenmode_coefficients`.
- `src/dft.cpp` - flux/mode monitor internals (`fields::add_dft_flux`, `fields::add_mode_monitor`).
- `src/near2far.cpp` - near-to-far internals (`fields::add_dft_near2far`, `get_farfields_array`).
- `src/sources.cpp` - source profile behavior used by many tutorials.
- `src/step.cpp` - timestep orchestration used across all examples.
- `scheme/meep.scm.in` - run helpers (`run-sources`, `run-sources+`, `stop-when-fields-decayed`).

## Fast source navigation
- `rg -n "add_flux|add_near2far|get_eigenmode_coefficients|run\(" python/simulation.py`
- `rg -n "add_dft_flux|add_mode_monitor|add_dft_near2far|get_farfields_array" src/dft.cpp src/near2far.cpp`
- `rg -n "define \(run-sources|define \(run-sources\+|stop-when-fields-decayed" scheme/meep.scm.in`
