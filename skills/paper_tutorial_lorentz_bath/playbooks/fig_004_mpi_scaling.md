# fig_004: MPI scaling benchmark for Lorentz vs Lorentz-Bath

## Scientific aim
Reproduce per-step strong-scaling behavior up to 240 CPUs and verify the expected order-of-magnitude runtime overhead for Lorentz-Bath relative to Lorentz.

## Runtime recipe
1. Create or reuse an execution folder under `projects/YYYY-MM-DD-<scope>/`:
```bash
RUN_DATE="$(date +%F)"
RUN_SCOPE="lorentz-bath"
RUN_DIR="projects/${RUN_DATE}-${RUN_SCOPE}"
mkdir -p "$RUN_DIR"
if [[ ! -d "$RUN_DIR/fdtd_bath" ]]; then
  cp -R skills/paper_tutorial_lorentz_bath/assets/fdtd_bath "$RUN_DIR/"
  cp -R skills/paper_tutorial_lorentz_bath/assets/scripts "$RUN_DIR/"
fi
```

2. Bootstrap trials 2-5 from the packaged trial-1 template:
```bash
cd "$RUN_DIR"
bash scripts/prepare_figure4_trials.sh fdtd_bath/implementation_2025/benchmark_performance_2d
```

3. Run and collect all benchmark trials (SLURM workflow):
```bash
for trial in 1 2 3 4 5; do
  cd "$RUN_DIR/fdtd_bath/implementation_2025/benchmark_performance_2d/trial_${trial}_standard_dynamical_chunk"
  bash submit_all_benchmark.sh
  # after jobs complete:
  bash collect_simulation_timesteps.sh
done
```

4. Generate benchmark figure:
```bash
cd "$RUN_DIR/fdtd_bath/implementation_2025/plotting"
python plot_mpi_benchmark.py
```

5. Validate scaling and overhead claims:
```bash
cd "$RUN_DIR"
python scripts/validate_figure4_timing.py \
  --benchmark-root fdtd_bath/implementation_2025/benchmark_performance_2d
```

## Validation checklist
- Each trial contains both `info_lb.txt` and `info_lorentz.txt` with CPU grid `1,2,4,8,16,24,48,72,96,192,240`.
- Aggregated timing curves are mostly non-increasing with CPU and remain within a broad near-ideal envelope.
- Median Lorentz-Bath/Lorentz overhead lies within the expected manuscript band (`10x` to `30x`, centered near `20x`).
- Output figure `implementation_2025/plotting/mpi_benchmark.pdf` is regenerated from the current trial timing tables.
