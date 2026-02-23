---
name: meep-fdtd-bath-implementation
description: "Use this skill for implementation-level analysis of the FDTD-Bath (Lorentz-Bath) extension in this modified Meep tree. Route here when requests need traceability from the paper \"FDTD with Auxiliary Bath Fields for Condensed-Phase Polaritonics: Fundamentals and Implementation\" to Python `BathLorentzianSusceptibility` (`python/geom.py`) and C++ `bath_lorentzian_susceptibility` (`src/susceptibility.cpp`), including wrapper/bridge files."
---

# meep: FDTD-Bath Implementation

## High-Signal Playbook
### Route conditions
- Use this skill when the request is about how the Lorentz-Bath model is implemented in this repository.
- Use `paper_tutorial_lorentz_bath` when the request is figure reproduction or benchmark reproduction for the first FDTD-Bath publication.
- Use `meep-inputs-and-modeling` for general material-modeling questions that are not specifically about the bath implementation internals.
- Use `meep-api-and-scripting` for broad API symbol lookup not centered on the Lorentz-Bath path.

### Triage questions
- Is the user asking for theory-to-code mapping, code-level debugging, or publication reproduction?
- Do they need Python API behavior, C++ update internals, or both?
- Is the question about parameter mapping (`num_bath`, `bath_form`, `bath_dephasing`, `bath_gamma`)?
- Is the question about extracting and analyzing bath arrays during simulation?

### Canonical workflow
1. Start from the equations and implementation claims in the paper "FDTD with Auxiliary Bath Fields for Condensed-Phase Polaritonics: Fundamentals and Implementation".
2. Map paper parameters to `BathLorentzianSusceptibility` construction in `python/geom.py`.
3. Trace Python susceptibility objects into C++ objects via `src/meepgeom.cpp`.
4. Confirm class interfaces in `src/meep.hpp`, then inspect update and memory logic in `src/susceptibility.cpp`.
5. If roundtrip serialization or Python object reconstruction matters, inspect `python/meep.i` and `python/typemap_utils.cpp`.
6. If runtime bath-field inspection is requested, inspect `python/simulation.py` and `src/array_slice.cpp`.
7. Route reproduction requests to `skills/paper_tutorial_lorentz_bath/SKILL.md`.

### Fast theory-to-code anchors
- Paper: "FDTD with Auxiliary Bath Fields for Condensed-Phase Polaritonics: Fundamentals and Implementation" - Lorentz-Bath finite-difference updates (bath and polarization equations) and implementation-details section.
- Paper: "FDTD with Auxiliary Bath Fields for Condensed-Phase Polaritonics: Fundamentals and Implementation" - Python usage snippets for explicit and simplified bath initialization.
- `python/geom.py`: `BathLorentzianSusceptibility.__init__`, `obtain_independent_lorentzians`, `bath_potential_energy_distribution`.
- `src/susceptibility.cpp`: `bath_lorentzian_susceptibility::new_internal_data`, `init_internal_data`, `copy_internal_data`, `update_P`, `dump_params`.

## Scope
- Explain and inspect Lorentz-Bath/FDTD-Bath implementation details in this modified Meep codebase.
- Keep answers grounded in both paper equations and concrete source paths.

## Primary documentation references
- `FDTD with Auxiliary Bath Fields for Condensed-Phase Polaritonics: Fundamentals and Implementation` (paper title reference)
- `skills/meep-fdtd-bath-implementation/references/doc_map.md`
- `skills/paper_tutorial_lorentz_bath/SKILL.md`

## Workflow
- Start with `skills/meep-fdtd-bath-implementation/references/doc_map.md`.
- Escalate to `skills/meep-fdtd-bath-implementation/references/source_map.md` for implementation-level details.
- Prefer targeted symbol search before broad file reads:
  - `rg -n "BathLorentzianSusceptibility|bath_form|bath_dephasing|obtain_independent_lorentzians" python/geom.py`
  - `rg -n "bath_lorentzian_susceptibility|update_P|new_internal_data|dump_params" src/susceptibility.cpp src/meep.hpp src/meepgeom.cpp`
  - `rg -n "get_bath_pol_array|get_bath_pol_array_slice|bath_base" python/simulation.py src/array_slice.cpp`

## Publication reproduction handoff
- For reproducing the first FDTD-Bath publication, use:
  - `skills/paper_tutorial_lorentz_bath/SKILL.md`
- Figure playbooks:
  - `skills/paper_tutorial_lorentz_bath/playbooks/fig_001_model_and_geometry.md`
  - `skills/paper_tutorial_lorentz_bath/playbooks/fig_002_1d_spectra_and_dynamics.md`
  - `skills/paper_tutorial_lorentz_bath/playbooks/fig_003_2d_normal_incidence_spectra.md`
  - `skills/paper_tutorial_lorentz_bath/playbooks/fig_004_mpi_scaling.md`

## Source entry points for unresolved issues (behavior-level)
- `python/geom.py` - `BathLorentzianSusceptibility` API surface and parameter initialization rules.
- `src/meepgeom.cpp` - branch selecting `new meep::bath_lorentzian_susceptibility(...)`.
- `src/meep.hpp` - `bath_lorentzian_susceptibility` declaration and bath-state method overrides.
- `src/susceptibility.cpp` - bath state allocation, time stepping, and parameter dump implementation.
- `python/meep.i` - Python export list for `BathLorentzianSusceptibility`.
- `python/typemap_utils.cpp` - susceptibility-to-Python object conversion for bath susceptibilities.
- `python/simulation.py` - `Simulation.get_bath_pol_array`.
- `src/array_slice.cpp` - bath polarization slice extraction via `get_bath_pol_array_slice`.
