---
name: paper_tutorial_lorentz_bath
description: Reproduce the manuscript's Lorentz vs Lorentz-Bath FDTD-Bath results (Fig. 1-4 geometry, 1D/2D spectra and dynamics, and MPI scaling) from packaged local assets, and adapt the same parameterized procedures to related condensed-phase polaritonic systems.
---

Use this skill from repository root. Run simulations and postprocessing only inside `projects/YYYY-MM-DD-<scope>/` runtime folders.

## Core Simulation Strategy
1. Install the modified Meep stack before any figure run using the packaged scripts in `assets/fdtd_bath/installation_scripts/`.
2. Stage `assets/fdtd_bath/` and `assets/scripts/` into a dated runtime folder under `projects/`; never run workloads inside `skills/paper_tutorial_lorentz_bath/`.
3. Keep manuscript baseline constants unless a playbook explicitly changes them: `num_bath=100`, `d=0.99`, `k=0.01`, `bath_width=10*gamma`, `gamma=0.04`, `nu0=1.0 um^-1`, mirror `n=10`, mirror thickness `0.02 um`, slab/cavity thickness `1.0 um`.
4. For each figure, execute the fixed sequence: simulation batch -> plotting script -> validation script.
5. Figure 1 schematic image regeneration is constrained by a missing upstream asset (`fdtd_bath_demo.png`), so publication-grade reproducibility is enforced through scripted geometry/model-validation outputs.

## Minimal Execution Recipes
1. Create run directory and stage assets:
```bash
RUN_DATE="$(date +%F)"
RUN_SCOPE="lorentz-bath"
RUN_DIR="projects/${RUN_DATE}-${RUN_SCOPE}"
mkdir -p "$RUN_DIR"
cp -R skills/paper_tutorial_lorentz_bath/assets/fdtd_bath "$RUN_DIR/"
cp -R skills/paper_tutorial_lorentz_bath/assets/scripts "$RUN_DIR/"
```

2. Install the modified Meep stack (choose one environment script):
```bash
cd "$RUN_DIR/fdtd_bath"
bash installation_scripts/meep_install_CentOS9.sh
# or
bash installation_scripts/meep_install_hpc_anvil.sh
```

3. Run direct figure workflows:
```bash
cd "$RUN_DIR"
python scripts/extract_fig1_geometry.py --script-1d fdtd_bath/implementation_2025/1d_harmonic_broadlinewidth/run_sc_fabry_perot_1d.py --script-2d fdtd_bath/implementation_2025/2d_harmonic_broadlinewidth/run_sc_fabry_perot_2d_spectrum.py --output fig_001_geometry_constants.json

cd "$RUN_DIR/fdtd_bath/implementation_2025/1d_harmonic_broadlinewidth"
bash submit_outcav.sh && bash submit_incav_lorentz.sh && bash submit_incav_lb_lorentzian.sh
cd ../plotting && python plot_1d_demo.py
cd "$RUN_DIR" && python scripts/validate_figure2_outputs.py --root fdtd_bath/implementation_2025/1d_harmonic_broadlinewidth

cd "$RUN_DIR/fdtd_bath/implementation_2025/2d_harmonic_broadlinewidth"
bash submit_incav_lorentz_scan_sigma.sh && bash submit_incav_lb_lorentzian_scan_sigma.sh
cd ../plotting && python plot_2d_demo.py
cd "$RUN_DIR" && python scripts/validate_figure3_vs_figure2.py --one-d-root fdtd_bath/implementation_2025/1d_harmonic_broadlinewidth --two-d-root fdtd_bath/implementation_2025/2d_harmonic_broadlinewidth

bash scripts/prepare_figure4_trials.sh fdtd_bath/implementation_2025/benchmark_performance_2d
for trial in 1 2 3 4 5; do
  cd "$RUN_DIR/fdtd_bath/implementation_2025/benchmark_performance_2d/trial_${trial}_standard_dynamical_chunk"
  bash submit_all_benchmark.sh
  bash collect_simulation_timesteps.sh
done
cd "$RUN_DIR/fdtd_bath/implementation_2025/plotting" && python plot_mpi_benchmark.py
cd "$RUN_DIR" && python scripts/validate_figure4_timing.py --benchmark-root fdtd_bath/implementation_2025/benchmark_performance_2d
```

## Figure Routing
- `fig_001`: Scope is model-definition and geometry lock-in for all downstream runs, including mirror/slab constants and boundary assumptions; playbook: `skills/paper_tutorial_lorentz_bath/playbooks/fig_001_model_and_geometry.md`.
- `fig_002`: Scope is 1D free-space vs in-cavity Lorentz/Lorentz-Bath spectral behavior and cavity-energy persistence at weak/strong coupling; playbook: `skills/paper_tutorial_lorentz_bath/playbooks/fig_002_1d_spectra_and_dynamics.md`.
- `fig_003`: Scope is 2D normal-incidence reproduction plus quantitative 1D-vs-2D LP/UP consistency checks; playbook: `skills/paper_tutorial_lorentz_bath/playbooks/fig_003_2d_normal_incidence_spectra.md`.
- `fig_004`: Scope is 2D MPI strong-scaling and Lorentz-Bath overhead quantification over 1 to 240 CPUs with five-trial aggregation; playbook: `skills/paper_tutorial_lorentz_bath/playbooks/fig_004_mpi_scaling.md`.

## Beyond Manuscript Exploration
- Sweep `num_bath` (`50, 100, 200`) while preserving linewidth matching (`gamma0 + bath_dephasing = gamma`) to test bath-discretization robustness.
- Sweep `bath_width/gamma` (`6, 10, 14`) and track LP/UP FWHM shifts with the same windows used by `scripts/validate_figure2_outputs.py`.
- Extend 2D scans to finite angles via `-c` in `run_sc_fabry_perot_2d_spectrum.py`, while keeping a normal-incidence control at `0 deg`.
- For performance studies beyond the paper, keep physics fixed and vary only MPI decomposition/chunk balancing before re-running `validate_figure4_timing.py`.
