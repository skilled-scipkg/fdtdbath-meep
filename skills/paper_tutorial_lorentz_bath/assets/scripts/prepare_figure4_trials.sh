#!/usr/bin/env bash
set -euo pipefail

BENCH_ROOT="${1:-implementation_2025/benchmark_performance_2d}"
TEMPLATE_DIR="$BENCH_ROOT/trial_1_standard_dynamical_chunk"

if [[ ! -d "$TEMPLATE_DIR" ]]; then
  echo "Template directory missing: $TEMPLATE_DIR" >&2
  exit 1
fi

for trial in 2 3 4 5; do
  target="$BENCH_ROOT/trial_${trial}_standard_dynamical_chunk"
  mkdir -p "$target"
  for file in \
    benchmark-lb-multinode.sh \
    benchmark-lorentz-multinode.sh \
    collect_simulation_timesteps.sh \
    run_benchmark_multinode_lb.sh \
    run_benchmark_multinode_lorentz.sh \
    run_benchmark_singlenode_lb.sh \
    run_benchmark_singlenode_lorentz.sh \
    run_sc_fabry_perot_2d_spectrum.py \
    submit_all_benchmark.sh; do
    cp "$TEMPLATE_DIR/$file" "$target/$file"
  done
  chmod +x "$target"/*.sh
  echo "Prepared $target"
done
