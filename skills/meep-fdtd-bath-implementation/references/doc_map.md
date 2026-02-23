# meep documentation map: FDTD-Bath Implementation

Use this map from `skills/meep-fdtd-bath-implementation/SKILL.md`.

## Paper anchors (theory and implementation narrative)
- Paper title: `FDTD with Auxiliary Bath Fields for Condensed-Phase Polaritonics: Fundamentals and Implementation`
- Section anchor: `The Lorentz-Bath model`
  - Finite-difference update equations for bath oscillators and polarization.
- Section anchor: `Implementation Details`
  - Design choice: derive C++ `bath_lorentzian_susceptibility` and Python `BathLorentzianSusceptibility`.
- Code-listing anchor: `LB_custom`
  - Explicit arrays API: `bath_frequencies`, `bath_gammas`, `bath_couplings`.
- Code-listing anchor: `LB_lorentzian`
  - Simplified API: `bath_form`, `bath_width`, `bath_gamma`, `bath_dephasing`.

## Publication reproduction references
- `skills/paper_tutorial_lorentz_bath/SKILL.md`
- `skills/paper_tutorial_lorentz_bath/references/doc_map.md`
- `skills/paper_tutorial_lorentz_bath/playbooks/fig_001_model_and_geometry.md`
- `skills/paper_tutorial_lorentz_bath/playbooks/fig_002_1d_spectra_and_dynamics.md`
- `skills/paper_tutorial_lorentz_bath/playbooks/fig_003_2d_normal_incidence_spectra.md`
- `skills/paper_tutorial_lorentz_bath/playbooks/fig_004_mpi_scaling.md`

## Escalation
If documentation context is insufficient, use `skills/meep-fdtd-bath-implementation/references/source_map.md`.
