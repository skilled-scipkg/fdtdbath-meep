---
name: meep-advanced-topics
description: Use this skill for specialized low-frequency topics in this modified Meep repository (includes Lorentz-Bath/FDTD-Bath support for condensed-phase polaritons), consolidated from one-doc skills (special kz, symmetry, PML deep-dive, Yee lattice, units/nonlinearity, Scheme/Guile notes, theory math note, requirements, licensing, acknowledgements). Keep routing compact and docs-first.
---

# meep: Advanced Topics

## Scope
- Handle specialized topics that are important but too narrow as standalone one-doc skills.
- Route quickly to the exact document path and only escalate to source when documentation is insufficient.

## Topic routing
- 2d out-of-plane wavevector (`kz`) behavior:
  - `doc/docs/2d_Cell_Special_kz.md`
- Symmetry usage and reductions:
  - `doc/docs/Exploiting_Symmetry.md`
- PML theory, breakdown cases, and practical caveats:
  - `doc/docs/Perfectly_Matched_Layer.md`
- Yee lattice interpretation details:
  - `doc/docs/Yee_Lattice.md`
- Units/nonlinearity interpretation:
  - `doc/docs/Units_and_Nonlinearity.md`
- Scheme/Guile interface background:
  - `doc/docs/Guile_and_Scheme_Information.md`
- Frequency-domain eigensolver math note:
  - `doc/docs/Eigensolver_Math.md`
- Packaging/document requirements list:
  - `doc/requirements.txt`
- Licensing/copyright and project acknowledgements:
  - `doc/docs/License_and_Copyright.md`
  - `doc/docs/Acknowledgements.md`

## Workflow
- Start with the exact topic document above.
- If needed, inspect `skills/meep-advanced-topics/references/doc_map.md` for the consolidated inventory.
- Escalate to `skills/meep-advanced-topics/references/source_map.md` only when behavior details are not documented.
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

## Source entry points for unresolved issues
- `python/simulation.py` - runtime options (`k_point`, boundary/symmetry-facing knobs).
- `scheme/meep.scm.in` - Scheme input variables (Courant, complex fields, `k-point`, material classes).
- `src/fields.cpp` - Yee-grid field update behavior.
- `src/update_eh.cpp` - electric/magnetic update paths relevant to boundary/PML behavior.
- `src/update_pols.cpp` - polarization/material update internals.
- `src/step.cpp` - timestep control and orchestration.
- `src/meepgeom.cpp` - geometry periodicity/discretization behavior.
- `src/vec.cpp` - coordinate/vector utility internals.
- Prefer targeted source search first: `rg -n "<symbol_or_keyword>" python scheme src`.
