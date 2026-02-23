---
name: meep-getting-started
description: Use this skill for first-run onboarding, core concepts, and minimal Python/Scheme simulations in this modified Meep repository (includes Lorentz-Bath/FDTD-Bath support for condensed-phase polaritons). Route users to build/install, modeling, workflow, API, or examples skills when requests go beyond quickstart scope.
---

# meep: Getting Started

Repository note: this repository is a modified Meep distribution with Lorentz-Bath/FDTD-Bath extensions for condensed-phase polaritons.

## High-Signal Playbook
### Route conditions
- Use `meep-build-and-install` for environment setup, dependency, compiler, MPI, or import errors.
- Use `meep-inputs-and-modeling` for advanced materials, boundary conditions, symmetry, and geometry parameterization.
- Use `meep-simulation-workflows` for monitor setup, run control, normalization, and convergence debugging.
- Use `meep-api-and-scripting` for Python/Scheme symbol lookup and interface-specific behavior.
- Use `meep-examples-and-tutorials` when the user mainly needs closest reference examples by problem type.

### Triage questions
- Is the user starting in Python or Scheme?
- Is `import meep` (Python) or `meep <file>.ctl` (Scheme) already working?
- What unit length `a` are they using, and what physical wavelength/frequency should that map to?
- Is this a 2d, 3d, or cylindrical (`m`) problem?
- Do they need fields, flux/transmittance, resonance, or far-field output first?
- Are they asking for a first runnable script or for interpretation of existing results?

### Canonical workflow
1. Confirm installation sanity with a minimal runtime check (`python -c 'import meep as mp; print(mp.__version__)'`) and route package/setup issues to `meep-build-and-install`.
2. Start from the straight-waveguide basics in Python or Scheme (`doc/docs/Python_Tutorials/Basics.md`, `doc/docs/Scheme_Tutorials/Basics.md`).
3. Build a minimal model: `cell`, `geometry`, `sources`, `pml_layers`, `resolution`, then run.
4. Validate that outputs are physically sensible (field propagation direction, PML absorption, no immediate blow-up).
5. Advance to bent-waveguide or bend-flux examples only after the straight case behaves correctly.
6. Escalate to source only if docs and examples do not explain behavior.

### Minimal working example
Python (from `doc/docs/Python_Tutorials/Basics.md`):
```python
import meep as mp

cell = mp.Vector3(16, 8, 0)
geometry = [mp.Block(mp.Vector3(mp.inf, 1, mp.inf),
                     center=mp.Vector3(),
                     material=mp.Medium(epsilon=12))]
sources = [mp.Source(mp.ContinuousSource(frequency=0.15),
                     component=mp.Ez,
                     center=mp.Vector3(-7, 0))]

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(1.0)],
    resolution=10,
)
sim.run(until=200)
```

Scheme (from `doc/docs/Scheme_Tutorials/Basics.md`):
```scheme
(set! geometry-lattice (make lattice (size 16 8 no-size)))
(set! geometry
      (list (make block (center 0 0)
                        (size infinity 1 infinity)
                        (material (make medium (epsilon 12))))))
(set! sources
      (list (make source
                  (src (make continuous-src (frequency 0.15)))
                  (component Ez)
                  (center -7 0))))
(set! pml-layers (list (make pml (thickness 1.0))))
(set! resolution 10)
(run-until 200 (at-beginning output-epsilon) (at-end output-efield-z))
```

### Pitfalls and fixes
- PML is inside the cell and overlaps objects by design; leave physical padding before the region of interest (`doc/docs/Python_Tutorials/Basics.md`).
- Abrupt CW turn-on excites extra frequencies; use finite source width / smooth turn-on from basics tutorials.
- `run`/`run-until` arguments are step functions, not arbitrary statements (`doc/docs/The_Run_Function_Is_Not_A_Loop.md`).
- Frequency units are inverse length in chosen scale (`f = a / lambda_vac`); wrong unit mapping causes wrong physics (`doc/docs/Introduction.md#units-in-meep`).
- Doubling resolution also tightens timestep via Courant relation and can greatly raise runtime (`doc/docs/Introduction.md`).

### Convergence and validation checks
- Sweep `resolution` and confirm key observable changes are monotonic/saturating.
- Sweep PML thickness and source-to-PML padding for reflection sensitivity.
- Use decay-based stopping (`stop_when_fields_decayed` in Python or `(stop-when-fields-decayed ...)` in Scheme) for spectral tasks.
- For reflectance/transmittance, keep normalization and scattering runs on identical discretization and monitor placement.

## Scope
- Handle first-run Meep questions, core concepts, and minimal simulation startup.
- Keep responses concise and docs-first; avoid source deep-dives unless docs are insufficient.

## Primary documentation references
- `doc/docs/Introduction.md`
- `doc/docs/Python_Tutorials/Basics.md`
- `doc/docs/Scheme_Tutorials/Basics.md`
- `doc/docs/C++_Tutorial.md`
- `doc/bfast/fixed_angle_broadband_simulations_in_Meep.md`

## Workflow
- Start from primary docs and runnable basics examples.
- If details are missing, inspect `skills/meep-getting-started/references/doc_map.md`.
- Only if ambiguity remains, inspect `skills/meep-getting-started/references/source_map.md` entry points.
- Cite exact file paths in responses.

## Tutorials and examples
- `python/examples/straight-waveguide.py`
- `python/examples/bent-waveguide.py`
- `scheme/examples/straight-waveguide.ctl`
- `scheme/examples/bent-waveguide.ctl`
- `doc/docs/Python_Tutorials`
- `doc/docs/Scheme_Tutorials`

## Test references
- `tests`
- `python/tests`

## Optional deeper inspection
- `libpympb`
- `python`
- `scheme`
- `src`

## Source entry points for unresolved issues (behavior-level)
- `python/simulation.py` - high-level `Simulation` setup and run semantics.
- `python/source.py` - source-time/source-shape behavior.
- `python/materials.py` - built-in material definitions and helpers.
- `scheme/meep.scm.in` - Scheme input-variable defaults and run-function glue.
- `src/fields.cpp` - Yee-grid field updates.
- `src/sources.cpp` - source current injection paths.
- `src/meepgeom.cpp` - geometry discretization and periodicity helpers.
- `src/material_data.cpp` - material storage/interpolation internals.
- `src/step.cpp` - timestep orchestration and stepping hooks.
- Prefer targeted search first: `rg -n "<symbol_or_keyword>" python scheme src`.
