---
name: meep-index
description: Use this skill as the docs-first router for all Meep requests. Classify intent into a topic skill first (getting started, build/install, modeling, workflows, API/scripting, examples, advanced topics, or FDTD-Bath implementation tracing) and escalate to source only when topic docs are insufficient.
---

# meep Skills Index

## Route the request
- Classify the request into one topic skill below before answering deeply.
- Keep answers docs-first; escalate to source only after topic docs and examples are insufficient.
- Prefer workflow-level guidance over exhaustive per-function detail unless explicitly requested.

## Topic skill dispatch
- `meep-getting-started`: first-run onboarding, units, and minimal Python/Scheme simulations.
- `meep-build-and-install`: Conda/source build paths, dependencies, platform caveats, and validation.
- `meep-inputs-and-modeling`: geometry/material/source/boundary modeling decisions.
- `meep-simulation-workflows`: setup/run/monitor/termination/convergence workflows.
- `meep-api-and-scripting`: Python/Scheme interface usage, symbol lookup, and binding behavior.
- `meep-examples-and-tutorials`: selecting the best runnable Python/Scheme reference example.
- `meep-advanced-topics`: consolidated specialist docs (special `kz`, symmetry, PML, Yee lattice, units/nonlinearity, Scheme notes, theory note, requirements/license/acknowledgements).
- `meep-fdtd-bath-implementation`: implementation-focused mapping for this modified Lorentz-Bath/FDTD-Bath code path (manuscript -> Python `BathLorentzianSusceptibility` -> C++ `bath_lorentzian_susceptibility`).
- `paper_tutorial_lorentz_bath`: advanced publication-reproduction workflow for Lorentz vs Lorentz-Bath FDTD-Bath figures (geometry, 1D/2D spectra-dynamics, and MPI scaling).

## Documentation-first inputs
- `doc`

## Tutorials and examples roots
- `python/examples`
- `scheme/examples`
- `doc/docs/Python_Tutorials`
- `doc/docs/Scheme_Tutorials`

## Test roots for behavior checks
- `tests`
- `python/tests`

## Escalate only when needed
- Start in the routed topic skill file:
  - `skills/meep-getting-started/SKILL.md`
  - `skills/meep-build-and-install/SKILL.md`
  - `skills/meep-inputs-and-modeling/SKILL.md`
  - `skills/meep-simulation-workflows/SKILL.md`
  - `skills/meep-api-and-scripting/SKILL.md`
  - `skills/meep-examples-and-tutorials/SKILL.md`
  - `skills/meep-advanced-topics/SKILL.md`
  - `skills/meep-fdtd-bath-implementation/SKILL.md`
  - `skills/paper_tutorial_lorentz_bath/SKILL.md`
- If needed, inspect the matching topic doc map:
  - `skills/meep-getting-started/references/doc_map.md`
  - `skills/meep-build-and-install/references/doc_map.md`
  - `skills/meep-inputs-and-modeling/references/doc_map.md`
  - `skills/meep-simulation-workflows/references/doc_map.md`
  - `skills/meep-api-and-scripting/references/doc_map.md`
  - `skills/meep-examples-and-tutorials/references/doc_map.md`
  - `skills/meep-advanced-topics/references/doc_map.md`
  - `skills/meep-fdtd-bath-implementation/references/doc_map.md`
- If docs are still insufficient, inspect the matching topic source map:
  - `skills/meep-getting-started/references/source_map.md`
  - `skills/meep-build-and-install/references/source_map.md`
  - `skills/meep-inputs-and-modeling/references/source_map.md`
  - `skills/meep-simulation-workflows/references/source_map.md`
  - `skills/meep-api-and-scripting/references/source_map.md`
  - `skills/meep-examples-and-tutorials/references/source_map.md`
  - `skills/meep-advanced-topics/references/source_map.md`
  - `skills/meep-fdtd-bath-implementation/references/source_map.md`
- Use targeted symbol search while inspecting source: `rg -n "<symbol_or_keyword>" libpympb python scheme src`.

## Source directories for deeper inspection
- `libpympb`
- `python`
- `scheme`
- `src`
