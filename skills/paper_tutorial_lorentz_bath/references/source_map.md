# Source Map: paper_tutorial_lorentz_bath

## Simulation drivers
- `skills/paper_tutorial_lorentz_bath/assets/fdtd_bath/implementation_2025/1d_harmonic_broadlinewidth/run_sc_fabry_perot_1d.py`
  - 1D free-space/cavity spectra and energy-dynamics outputs.
- `skills/paper_tutorial_lorentz_bath/assets/fdtd_bath/implementation_2025/2d_harmonic_broadlinewidth/run_sc_fabry_perot_2d_spectrum.py`
  - 2D normal-incidence spectra outputs.
- `skills/paper_tutorial_lorentz_bath/assets/fdtd_bath/implementation_2025/benchmark_performance_2d/trial_1_standard_dynamical_chunk/run_sc_fabry_perot_2d_spectrum.py`
  - High-resolution benchmark driver for per-step scaling.

## Submission wrappers
- 1D: `submit_outcav.sh`, `submit_incav_lorentz.sh`, `submit_incav_lb_lorentzian.sh`
- 2D: `submit_incav_lorentz_scan_sigma.sh`, `submit_incav_lb_lorentzian_scan_sigma.sh`
- Benchmark: `submit_all_benchmark.sh`, `run_benchmark_*`, `benchmark-*.sh`, `collect_simulation_timesteps.sh`

## Plotting and validation
- Plotting scripts:
  - `skills/paper_tutorial_lorentz_bath/assets/fdtd_bath/implementation_2025/plotting/plot_1d_demo.py`
  - `skills/paper_tutorial_lorentz_bath/assets/fdtd_bath/implementation_2025/plotting/plot_2d_demo.py`
  - `skills/paper_tutorial_lorentz_bath/assets/fdtd_bath/implementation_2025/plotting/plot_mpi_benchmark.py`
  - `skills/paper_tutorial_lorentz_bath/assets/fdtd_bath/implementation_2025/plotting/columnplots.py`
- Validation scripts:
  - `skills/paper_tutorial_lorentz_bath/assets/scripts/extract_fig1_geometry.py`
  - `skills/paper_tutorial_lorentz_bath/assets/scripts/validate_figure2_outputs.py`
  - `skills/paper_tutorial_lorentz_bath/assets/scripts/validate_figure3_vs_figure2.py`
  - `skills/paper_tutorial_lorentz_bath/assets/scripts/validate_figure4_timing.py`
