---
name: meep-examples-and-tutorials
description: Use this skill to route users to the closest runnable examples/tutorials (Python and Scheme) in this modified Meep repository (includes Lorentz-Bath/FDTD-Bath support for condensed-phase polaritons), with minimal adaptation steps and docs-first references before source escalation.
---

# meep: Examples and Tutorials

## High-Signal Playbook
### Route conditions
- Use `meep-getting-started` for brand-new users needing first concepts/units before choosing examples.
- Use `meep-simulation-workflows` for monitor placement, run control, and convergence troubleshooting.
- Use `meep-api-and-scripting` for unresolved API symbols while adapting examples.
- Use `meep-inputs-and-modeling` for advanced material/geometry design choices.

### Triage questions
- What physics outcome is needed (transmission, resonances, far field, nonlinear, forces, optimization)?
- Python or Scheme preference?
- Is this a minimal reproducible baseline or a production-scale adaptation?
- Broadband pulse or single-frequency/CW objective?
- Is MPI required now or only after serial validation?
- What metric should be matched first (flux, Q, mode coefficient, pattern)?

### Canonical workflow
1. Map user intent to the smallest relevant tutorial/example pair.
2. Run that example unmodified first (serial).
3. Confirm expected output artifact (printed metric, HDF5, PNG, spectrum).
4. Modify one parameter block at a time (geometry, source, monitors, then runtime).
5. Re-run and compare against baseline metric after each change.
6. Escalate to docs/source only if adaptation behavior diverges from expected.

### Minimal working example
Fast baseline commands:
```bash
python python/examples/straight-waveguide.py
meep scheme/examples/straight-waveguide.ctl
```

Parallel follow-up command (after serial sanity):
```bash
mpirun -np 4 python -m mpi4py python/examples/bend-flux.py
```

### Pitfalls and fixes
- Starting from a complex tutorial (adjoint/metasurface) before a basic baseline slows debugging.
- Reflectance/transmittance scripts generally need normalization/reference runs; skipping them gives wrong spectra.
- Using a CW source for broadband tasks causes non-convergent Fourier interpretation; prefer Gaussian pulse for spectra.
- Runtime too short smears narrow spectral features (Fourier uncertainty).
- Near2far surfaces crossing PML or in inhomogeneous media invalidate far-field transforms.
- Mode-decomposition parity constraints must match geometry/source symmetry assumptions.

### Convergence and validation checks
- Sweep resolution and ensure key spectra/modes stabilize.
- Sweep PML thickness/padding for scattering/radiation problems.
- Increase runtime or decay strictness for narrowband resonances.
- Validate conservation-style checks where applicable (incident vs reflected/transmitted/loss).
- Compare adapted results to documented reference trends before adding extra complexity.

## Scope
- Handle routing to runnable examples and tutorial references.
- Prefer minimal reproducible examples over abstract explanation when user asks "how".

## Primary documentation references
- `python/examples/README.md`
- `doc/docs/Python_Tutorials/Basics.md`
- `doc/docs/Scheme_Tutorials/Basics.md`
- `doc/docs/Python_Tutorials/Mode_Decomposition.md`
- `doc/docs/Scheme_Tutorials/Mode_Decomposition.md`
- `doc/docs/Python_Tutorials/Near_to_Far_Field_Spectra.md`
- `doc/docs/Scheme_Tutorials/Near_to_Far_Field_Spectra.md`
- `doc/docs/Python_Tutorials/Frequency_Domain_Solver.md`
- `doc/docs/Scheme_Tutorials/Frequency_Domain_Solver.md`

## Example map by problem type
- Basics / first field plots:
  - Python: `python/examples/straight-waveguide.py`, `python/examples/bent-waveguide.py`, `python/examples/bend-flux.py`
  - Scheme: `scheme/examples/straight-waveguide.ctl`, `scheme/examples/bent-waveguide.ctl`, `scheme/examples/bend-flux.ctl`
- Resonances / cavities / bands:
  - Python: `python/examples/ring.py`, `python/examples/holey-wvg-cavity.py`, `python/examples/holey-wvg-bands.py`
  - Scheme: `scheme/examples/ring.ctl`, `scheme/examples/holey-wvg-cavity.ctl`, `scheme/examples/holey-wvg-bands.ctl`
- Mode decomposition / eigenmodes / diffraction:
  - Python: `python/examples/mode-decomposition.py`, `python/examples/oblique-source.py`, `python/examples/binary_grating.py`
  - Scheme: `scheme/examples/mode-decomposition.ctl`, `scheme/examples/oblique-source.ctl`, `scheme/examples/binary_grating.ctl`
- Near-to-far radiation:
  - Python: `python/examples/antenna-radiation.py`, `python/examples/cavity-farfield.py`, `python/examples/binary_grating_n2f.py`
  - Scheme: `scheme/examples/antenna-radiation.ctl`, `scheme/examples/cavity-farfield.ctl`, `scheme/examples/binary_grating_n2f.ctl`
- Materials / nonlinear / gyrotropic / multilevel:
  - Python: `python/examples/material-dispersion.py`, `python/examples/3rd-harm-1d.py`, `python/examples/faraday-rotation.py`, `python/examples/multilevel-atom.py`
  - Scheme: `scheme/examples/material-dispersion.ctl`, `scheme/examples/3rd-harm-1d.ctl`, `scheme/examples/faraday-rotation.ctl`, `scheme/examples/multilevel-atom.ctl`
- Frequency-domain solver / optimization:
  - Python: `python/examples/solve-cw.py`, `python/examples/adjoint_optimization/*`
  - Scheme: `scheme/examples/solve-cw.ctl`

## Workflow
- Start from primary docs and example map above.
- If details are missing, inspect `skills/meep-examples-and-tutorials/references/doc_map.md`.
- Escalate to `skills/meep-examples-and-tutorials/references/source_map.md` only for unresolved behavior.
- Cite exact file paths in responses.

## Test references
- `tests`
- `python/tests`

## Optional deeper inspection
- `libpympb`
- `python`
- `scheme`
- `src`

## Source entry points for unresolved issues (behavior-level)
- `python/examples/README.md` - high-level tutorial inventory.
- `python/examples` - canonical runnable Python scripts.
- `scheme/examples` - canonical runnable Scheme control files.
- `python/simulation.py` - API semantics behind example calls.
- `python/solver.py` - solver/mode helper behavior.
- `scheme/meep.scm.in` - Scheme interface defaults and run-function support.
- Prefer targeted search first: `rg -n "<symbol_or_keyword>" python/examples scheme/examples python scheme src`.
