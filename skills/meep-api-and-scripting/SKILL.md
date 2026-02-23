---
name: meep-api-and-scripting
description: Use this skill for Python/Scheme interface usage, scripting patterns, cross-interface mapping, and API symbol lookup strategy in this modified Meep repository (includes Lorentz-Bath/FDTD-Bath support for condensed-phase polaritons). Start from docs and examples, then escalate to binding/source files only for unresolved behavior.
---

# meep: API and Scripting

## High-Signal Playbook
### Route conditions
- Use `meep-getting-started` for first-run onboarding and minimal first simulation setup.
- Use `meep-simulation-workflows` for run/monitor/convergence workflows rather than symbol lookup.
- Use `meep-inputs-and-modeling` for material and geometry modeling decisions.
- Use `meep-build-and-install` for missing module, build, or dependency issues.

### Triage questions
- Is the user writing Python or Scheme control files?
- Are they asking for an API symbol, usage pattern, or behavior explanation?
- Do they need mode decomposition/eigenmode/adjoint/frequency-domain features?
- Is MPB availability relevant (eigenmodes, mode coefficients)?
- Is this a translation request between Scheme and Python APIs?
- Do they need docs-level guidance or source-level implementation details?

### Canonical workflow
1. Choose interface: Python first for new work (`doc/docs/Python_User_Interface.md`), Scheme for legacy control files (`doc/docs/Scheme_User_Interface.md`).
2. Resolve symbol in docs first (class/function sections, tutorial usage).
3. Confirm usage against a minimal example script.
4. If symbol behavior is unclear, inspect binding files (`python/meep.i`, `scheme/meep.i`) and high-level wrappers.
5. Escalate to C++ implementation only when docs+bindings do not resolve ambiguity.
6. Return answer with exact file-path references and shortest working pattern.

### Minimal working example
Python API baseline:
```python
import meep as mp

sim = mp.Simulation(
    cell_size=mp.Vector3(8, 4, 0),
    geometry=[],
    sources=[mp.Source(mp.GaussianSource(0.2, fwidth=0.1),
                       component=mp.Ez,
                       center=mp.Vector3(-2, 0))],
    boundary_layers=[mp.PML(1.0)],
    resolution=10,
)
sim.run(until=50)
```

Scheme control-file baseline:
```scheme
(set! geometry-lattice (make lattice (size 8 4 no-size)))
(set! sources (list
  (make source
    (src (make gaussian-src (frequency 0.2) (fwidth 0.1)))
    (component Ez)
    (center -2 0))))
(set! pml-layers (list (make pml (thickness 1.0))))
(set! resolution 10)
(run-until 50)
```

### Pitfalls and fixes
- Scheme interface is deprecated; prefer Python unless a legacy `.ctl` workflow is required (`doc/docs/Scheme_User_Interface.md`).
- `run` step functions must be function objects, not immediate function-call results (`doc/docs/The_Run_Function_Is_Not_A_Loop.md`).
- Source `component` labels denote current-source components, not hard-fixed fields (`doc/docs/Scheme_User_Interface.md#source`).
- `custom-src` in Scheme needs explicit `end-time` when using `run-sources` style logic.
- Eigenmode source and mode decomposition rely on MPB and parity conventions; check docs before interpreting coefficients.
- Scheme frequency-domain solver requires complex fields (`force-complex-fields? true`) (`doc/docs/Scheme_User_Interface.md`).

### Convergence and validation checks
- For equivalent Python/Scheme scripts, compare the same scalar metric (e.g., flux peak) under identical discretization.
- For nonzero `k_point` or cylindrical `m`, expect complex fields and higher memory use; validate run resources.
- For mode decomposition tasks, cross-check total flux from `get_fluxes` with mode-coefficient-derived power where applicable.
- For CW/frequency-domain tasks, tighten tolerances and verify field/error trends with resolution.

## Scope
- Handle language bindings, API lookup, and scripting-interface usage.
- Keep guidance docs-first and implementation-aware when needed.

## Primary documentation references
- `doc/docs/Python_User_Interface.md`
- `doc/docs/Scheme_User_Interface.md`
- `doc/docs/Mode_Decomposition.md`
- `doc/docs/Python_Tutorials/Eigenmode_Source.md`
- `doc/docs/Python_Tutorials/Frequency_Domain_Solver.md`
- `doc/docs/Python_Tutorials/Adjoint_Solver.md`
- `python/examples/README.md`

## Workflow
- Resolve via docs first.
- Use `skills/meep-api-and-scripting/references/doc_map.md` for complete topic inventory.
- Escalate to `skills/meep-api-and-scripting/references/source_map.md` only when docs are insufficient.
- Cite exact file paths in responses.

## Tutorials and examples
- `python/examples`
- `scheme/examples`
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
- `python/simulation.py` - high-level Python API surface and run/monitor plumbing.
- `python/source.py` - source classes and time-profile helpers.
- `python/solver.py` - solver wrappers and mode-related helpers.
- `python/meep.i` - SWIG interface definitions for Python bindings.
- `python/meep-python.hpp` - Python/C++ bridge declarations.
- `scheme/meep.scm.in` - Scheme-side input vars, defaults, and wrappers.
- `scheme/meep.i` - SWIG Scheme bindings.
- `scheme/meep.cpp` - Scheme/C++ binding implementation glue.
- `scheme/meep_op_renames.i` - compatibility symbol mapping.
- `src/meep.hpp` - core C++ API declarations.
- `src/meep_internals.hpp` - lower-level runtime/internal interfaces.
- Prefer targeted search first: `rg -n "<symbol_or_keyword>" python scheme src`.
