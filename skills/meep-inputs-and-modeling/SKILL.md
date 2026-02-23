---
name: meep-inputs-and-modeling
description: Use this skill for geometry/material/source/boundary modeling decisions in this modified Meep repository (including Lorentz-Bath/FDTD-Bath support for condensed-phase polaritons, plus dispersion, nonlinearity, symmetry, chunks, and coordinate choices). Keep guidance docs-first and escalate to source only when docs are insufficient.
---

# meep: Inputs and Modeling

## High-Signal Playbook
### Route conditions
- Use `meep-getting-started` for initial onboarding and first minimal run.
- Use `meep-simulation-workflows` for run control, monitors, termination, and output interpretation.
- Use `meep-api-and-scripting` for unresolved Python/Scheme API symbol semantics.
- Use `meep-build-and-install` for environment and dependency failures.

### Triage questions
- What structure class is being modeled (waveguide, cavity, grating, scattering, periodic)?
- Are materials nondispersive, dispersive, nonlinear, gyrotropic, or multilevel atomic?
- What dimensionality is intended (1d/2d/3d/cylindrical), and is symmetry exploitable?
- Are boundaries metallic, Bloch (`k_point`), or PML?
- Is geometry analytic objects, material function, or imported layout (e.g., GDSII)?
- Which outputs define correctness (flux, fields, Q, far-field, LDOS)?

### Canonical workflow
1. Choose unit length `a` and convert desired wavelengths/frequencies to Meep units.
2. Define geometry and materials from simplest faithful representation first.
3. Set boundaries/symmetries (`PML`, mirror/rotation, `k_point`) consistent with physics.
4. Choose source type (`GaussianSource`, `ContinuousSource`, eigenmode source) aligned to objective.
5. Run a baseline with conservative resolution/PML settings.
6. Sweep dominant error knobs (resolution, PML thickness, runtime, source bandwidth/cutoff).
7. Lock a converged setup before adding advanced features (dispersion/nonlinearity/imported CAD).

### Minimal working example
Python baseline for geometry/material/source setup:
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
    resolution=20,
)
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(5, 0), 1e-6))
```

### Pitfalls and fixes
- Discontinuous interfaces often dominate error; raise resolution and use subpixel smoothing where applicable (`doc/docs/Subpixel_Smoothing.md`).
- Subpixel averaging primarily smooths instantaneous ε/μ; dispersive-interface convergence can remain first-order.
- Material functions with nonlinear/dispersion in Scheme require `extra-materials` hints (`doc/docs/Scheme_User_Interface.md`).
- Nonzero `k_point` (or cylindrical `m`) implies complex fields and higher memory cost.
- PML is inside the cell; insufficient padding can contaminate near fields and spectra.
- Overlapping geometry order matters; later objects in the list take precedence.

### Convergence and validation checks
- Resolution sweep against target metrics (not just visual fields).
- PML thickness/padding sweep for boundary-reflection sensitivity.
- Source bandwidth/cutoff sweep for spectral robustness.
- Symmetry on/off comparison for one validation run to catch parity mistakes.
- For imported geometries (GDSII/pixel/material functions), verify that model changes continuously under small parameter perturbations.

## Scope
- Handle inputs, model construction, and physical parameterization.
- Keep responses concise and numerically grounded.

## Primary documentation references
- `doc/docs/Materials.md`
- `doc/docs/Subpixel_Smoothing.md`
- `doc/docs/Scheme_User_Interface.md`
- `doc/docs/Chunks_and_Symmetry.md`
- `doc/docs/Python_Tutorials/Material_Dispersion.md`
- `doc/docs/Scheme_Tutorials/Material_Dispersion.md`
- `doc/docs/Python_Tutorials/Third_Harmonic_Generation.md`
- `doc/docs/Python_Tutorials/GDSII_Import.md`
- `doc/docs/Python_Tutorials/Mode_Decomposition.md`

## Workflow
- Start from primary docs and closest tutorial examples.
- If details are missing, inspect `skills/meep-inputs-and-modeling/references/doc_map.md`.
- Escalate to `skills/meep-inputs-and-modeling/references/source_map.md` only for unresolved implementation behavior.
- Cite exact file paths in responses.

## Tutorials and examples
- `python/examples/material-dispersion.py`
- `python/examples/3rd-harm-1d.py`
- `python/examples/coupler.py`
- `scheme/examples/material-dispersion.ctl`
- `scheme/examples/3rd-harm-1d.ctl`
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
- `python/materials.py` - Python material library and helper definitions.
- `src/material_data.hpp` - material data structures and interfaces.
- `src/material_data.cpp` - material interpolation/evaluation internals.
- `src/structure.cpp` - structure assembly on the Yee grid.
- `src/structure_dump.cpp` - structure export/debugging paths.
- `src/meepgeom.hpp` - geometry representation APIs.
- `src/meepgeom.cpp` - geometry discretization implementation.
- `src/fix_boundary_sources.cpp` - source handling near boundaries/PML.
- `src/GDSIIgeom.cpp` - GDSII geometry import support.
- `scheme/materials.scm` - Scheme material definitions.
- `scheme/structure.cpp` - Scheme-side structure bridge.
- Prefer targeted search first: `rg -n "<symbol_or_keyword>" python scheme src`.
