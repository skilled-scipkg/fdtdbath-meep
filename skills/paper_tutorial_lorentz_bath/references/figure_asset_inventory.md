# Figure Asset Inventory

## fig_001
- `implementation_2025/1d_harmonic_broadlinewidth/run_sc_fabry_perot_1d.py`
- `implementation_2025/2d_harmonic_broadlinewidth/run_sc_fabry_perot_2d_spectrum.py`
- Validation helper: `scripts/extract_fig1_geometry.py`
- Note: `fdtd_bath_demo.png` is referenced by the manuscript but not present in the packaged assets.

## fig_002
- `implementation_2025/1d_harmonic_broadlinewidth/run_sc_fabry_perot_1d.py`
- `implementation_2025/1d_harmonic_broadlinewidth/submit_outcav.sh`
- `implementation_2025/1d_harmonic_broadlinewidth/submit_incav_lorentz.sh`
- `implementation_2025/1d_harmonic_broadlinewidth/submit_incav_lb_lorentzian.sh`
- `implementation_2025/plotting/plot_1d_demo.py`
- `implementation_2025/plotting/columnplots.py`
- Validation helper: `scripts/validate_figure2_outputs.py`

## fig_003
- `implementation_2025/2d_harmonic_broadlinewidth/run_sc_fabry_perot_2d_spectrum.py`
- `implementation_2025/2d_harmonic_broadlinewidth/submit_incav_lorentz_scan_sigma.sh`
- `implementation_2025/2d_harmonic_broadlinewidth/submit_incav_lb_lorentzian_scan_sigma.sh`
- `implementation_2025/plotting/plot_2d_demo.py`
- `implementation_2025/plotting/columnplots.py`
- Validation helper: `scripts/validate_figure3_vs_figure2.py`

## fig_004
- `implementation_2025/benchmark_performance_2d/trial_1_standard_dynamical_chunk/run_sc_fabry_perot_2d_spectrum.py`
- `implementation_2025/benchmark_performance_2d/trial_1_standard_dynamical_chunk/benchmark-lorentz-multinode.sh`
- `implementation_2025/benchmark_performance_2d/trial_1_standard_dynamical_chunk/benchmark-lb-multinode.sh`
- `implementation_2025/benchmark_performance_2d/trial_1_standard_dynamical_chunk/run_benchmark_singlenode_lorentz.sh`
- `implementation_2025/benchmark_performance_2d/trial_1_standard_dynamical_chunk/run_benchmark_singlenode_lb.sh`
- `implementation_2025/benchmark_performance_2d/trial_1_standard_dynamical_chunk/run_benchmark_multinode_lorentz.sh`
- `implementation_2025/benchmark_performance_2d/trial_1_standard_dynamical_chunk/run_benchmark_multinode_lb.sh`
- `implementation_2025/benchmark_performance_2d/trial_1_standard_dynamical_chunk/submit_all_benchmark.sh`
- `implementation_2025/benchmark_performance_2d/trial_1_standard_dynamical_chunk/collect_simulation_timesteps.sh`
- `implementation_2025/benchmark_performance_2d/trial_{1..5}_standard_dynamical_chunk/info_lb.txt`
- `implementation_2025/benchmark_performance_2d/trial_{1..5}_standard_dynamical_chunk/info_lorentz.txt`
- `implementation_2025/plotting/plot_mpi_benchmark.py`
- `implementation_2025/plotting/columnplots.py`
- Helpers: `scripts/prepare_figure4_trials.sh`, `scripts/validate_figure4_timing.py`
