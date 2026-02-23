---
name: meep-simulation-workflows
description: Use this skill for end-to-end simulation execution flow (setup, run control, monitors, normalization, parallel run behavior, termination criteria, and output interpretation). Start with docs, escalate to source only for unresolved behavior.
---

# meep: Simulation Workflows

## High-Signal Playbook
### Route conditions
- Use `meep-inputs-and-modeling` for deep material/geometry parameterization choices.
- Use `meep-build-and-install` for environment, dependency, and compiler/MPI setup failures.
- Use `meep-api-and-scripting` for unresolved Python/Scheme API symbol behavior.
- Use `meep-examples-and-tutorials` when the user needs closest runnable references first.

### Triage questions
- What is the observable: fields, flux/transmittance, modes, forces, near-to-far, or LDOS?
- Python or Scheme interface?
- Is this single run, normalization/scattering pair, or parameter sweep?
- What are `resolution`, PML thickness/padding, and source bandwidth?
- What is the run stop condition (`until`, `until_after_sources`, decay condition)?
- Serial or MPI run, and does behavior differ by process count?

### Canonical workflow
1. Define baseline model: `cell`, geometry/materials, source, PML, and resolution.
2. Add the minimal monitor set (`add_flux`, `add_near2far`, field outputs) before running.
3. Run with explicit termination logic (prefer decay-based stop for spectral tasks).
4. For reflectance/transmittance, perform normalization and scattering runs with identical discretization and monitor placement.
5. Use docs-backed monitor constraints (e.g., near2far surfaces in homogeneous isotropic medium, not in PML).
6. Validate outputs for physical consistency (energy flow sign, expected symmetry/mode content).
7. Only then scale to MPI or large parameter sweeps.

### Minimal working example
Python flux workflow (docs-aligned structure):
```python
import meep as mp

cell = mp.Vector3(16, 8, 0)
geometry = [mp.Block(mp.Vector3(mp.inf, 1, mp.inf),
                     center=mp.Vector3(),
                     material=mp.Medium(epsilon=12))]
sources = [mp.Source(mp.GaussianSource(0.15, fwidth=0.1),
                     component=mp.Ez,
                     center=mp.Vector3(-7, 0))]

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(1.0)],
    resolution=10,
)

flux = sim.add_flux(0.15, 0.1, 100,
                    mp.FluxRegion(center=mp.Vector3(5, 0), size=mp.Vector3(0, 6)))

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(5, 0), 1e-6))
print(mp.get_fluxes(flux)[0])
```

MPI run form (from `doc/docs/Parallel_Meep.md`):
```bash
mpirun -np 4 python -m mpi4py your_script.py
```

### Pitfalls and fixes
- Treating `run`/`run-until` as a loop and passing evaluated expressions instead of step functions (`doc/docs/The_Run_Function_Is_Not_A_Loop.md`, `doc/docs/Field_Functions.md`).
- Stopping too early causes unconverged DFT/flux and spectral leakage; use decay-based criteria (`doc/docs/Python_User_Interface.md#stop_when_fields_decayed`).
- Mixing normalization and scattering setups (different monitor position/resolution/source profile) invalidates reflected/transmitted subtraction.
- Near-to-far regions must be in homogeneous isotropic medium and outside PML (`doc/docs/Python_User_Interface.md`, `doc/docs/Scheme_User_Interface.md`).
- Using too many MPI processes for a small job can hurt performance; scale process count with problem size (`doc/docs/Parallel_Meep.md`).

### Convergence and validation checks
- Sweep `resolution` and monitor dominant metrics (flux peaks, Q, far-field pattern).
- Sweep PML thickness and source/monitor offsets from boundaries.
- Tighten decay conditions (`decay_by`) and source `cutoff` for robust DFT convergence.
- Ensure monitor frequency span stays inside source bandwidth.
- Check conservation-style sanity where applicable (e.g., incident ~= reflected + transmitted + loss).

## Scope
- Handle simulation setup, execution flow, run control, and output interpretation.
- Keep responses abstract and workflow-oriented unless details are requested.

## Primary documentation references
- `doc/docs/Python_User_Interface.md`
- `doc/docs/Scheme_User_Interface.md`
- `doc/docs/The_Run_Function_Is_Not_A_Loop.md`
- `doc/docs/Field_Functions.md`
- `doc/docs/Synchronizing_the_Magnetic_and_Electric_Fields.md`
- `doc/docs/Parallel_Meep.md`
- `doc/docs/Python_Tutorials/Near_to_Far_Field_Spectra.md`
- `doc/docs/Scheme_Tutorials/Near_to_Far_Field_Spectra.md`

## Workflow
- Start from primary docs and closest tutorial examples.
- If details are missing, inspect `skills/meep-simulation-workflows/references/doc_map.md`.
- Escalate to `skills/meep-simulation-workflows/references/source_map.md` only for unresolved implementation behavior.
- Cite exact file paths in responses.

## Tutorials and examples
- `python/examples/bend-flux.py`
- `python/examples/antenna-radiation.py`
- `python/examples/cavity-farfield.py`
- `scheme/examples/bend-flux.ctl`
- `scheme/examples/antenna-radiation.ctl`
- `scheme/examples/cavity-farfield.ctl`

## Test references
- `tests`
- `python/tests`

## Optional deeper inspection
- `libpympb`
- `python`
- `scheme`
- `src`

## Source entry points for unresolved issues (behavior-level)
- `python/simulation.py` - `Simulation.run`, monitor APIs, stopping helpers.
- `src/step.cpp` - main stepping loop and run-function orchestration.
- `src/step_generic.cpp` - generic stepping path details.
- `src/time.cpp` - simulation-time bookkeeping.
- `src/fields.cpp` - field update kernels.
- `src/cw_fields.cpp` - CW / frequency-domain-related field paths.
- `src/near2far.cpp` - near-to-far transform implementation.
- `src/h5fields.cpp` - HDF5 output/update behavior.
- `src/loop_in_chunks.cpp` - chunked iteration and parallel workflow internals.
- Prefer targeted search first: `rg -n "<symbol_or_keyword>" python scheme src`.
