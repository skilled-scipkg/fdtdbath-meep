# meep source map: FDTD-Bath Implementation

Use this map only after reading `skills/meep-fdtd-bath-implementation/references/doc_map.md`.

## Starter commands
- `rg -n "BathLorentzianSusceptibility|bath_form|bath_dephasing|obtain_independent_lorentzians" python/geom.py`
- `rg -n "bath_lorentzian_susceptibility|new_internal_data|init_internal_data|copy_internal_data|update_P|dump_params" src/susceptibility.cpp src/meep.hpp`
- `rg -n "num_bath > 0|bath_lorentzian_susceptibility" src/meepgeom.cpp python/typemap_utils.cpp`
- `rg -n "get_bath_pol_array|get_bath_pol_array_slice|bath_base" python/simulation.py src/array_slice.cpp`

## Python API implementation
- `python/geom.py`
  - `class BathLorentzianSusceptibility(LorentzianSusceptibility)`
  - `__init__`: direct bath definition (`bath_frequencies`, `bath_gammas`, `bath_couplings`) and indirect forms (`bath_form`, `bath_width`, `bath_dephasing`, `bath_gamma`).
  - `obtain_independent_lorentzians`: state-space conversion and eigendecomposition helper.
  - `bath_potential_energy_distribution`: postprocessing helper for bath energy analysis.

## C++ runtime path
- `src/meepgeom.cpp`
  - In `geom_epsilon::add_susceptibilities`, the `ss->num_bath > 0` branch instantiates `meep::bath_lorentzian_susceptibility`.
- `src/meep.hpp`
  - `class bath_lorentzian_susceptibility`: constructor fields and overridden bath-specific memory/update methods.
- `src/susceptibility.cpp`
  - `new_internal_data`, `init_internal_data`, `copy_internal_data`: layout/allocation for polarization plus bath state.
  - `update_P`: coupled bright-mode plus bath-oscillator timestep kernel.
  - `dump_params`: HDF5 serialization including per-bath parameters.

## Binding and reconstruction path
- `python/meep.i`
  - Exposes `BathLorentzianSusceptibility` from `python/geom.py` in the Python module namespace.
- `python/typemap_utils.cpp`
  - `susceptibility_to_py_obj`: branch mapping C++ susceptibility structs back to Python `BathLorentzianSusceptibility` objects.

## Bath-field extraction path
- `python/simulation.py`
  - `Simulation.get_bath_pol_array`: Python-facing array retrieval for bath polarization fields.
- `src/array_slice.cpp`
  - `fields::get_bath_pol_array_slice`: C++ slice backend for bath arrays.
  - `bath_base`: helper locating bath data block in internal polarization storage.

## Validation checkpoints
- `num_bath > 0` selects bath susceptibility path end-to-end (constructor dispatch in `src/meepgeom.cpp`).
- Retrieved bath array shape is `[nbath] + spatial_dims` in `python/simulation.py`.
- `BathLorentzianSusceptibility` parameters in Python map to C++ vectors (`bath_frequencies`, `bath_couplings`, `bath_gammas`, `bath_anharmonicities`).
