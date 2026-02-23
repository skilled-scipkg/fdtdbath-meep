#!/usr/bin/env python3
"""Validate figure-4 MPI scaling tables and model-overhead trends."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

EXPECTED_CPUS = np.array([1, 2, 4, 8, 16, 24, 48, 72, 96, 192, 240], dtype=float)


def read_info(path: Path) -> np.ndarray:
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != 2:
        raise AssertionError(f"{path}: expected 2 columns")
    if data.shape[0] != EXPECTED_CPUS.shape[0]:
        raise AssertionError(f"{path}: expected {EXPECTED_CPUS.shape[0]} rows")
    if not np.allclose(data[:, 0], EXPECTED_CPUS):
        raise AssertionError(f"{path}: CPU grid mismatch")
    return data


def nonincreasing_fraction(x: np.ndarray) -> float:
    valid = np.isfinite(x)
    x = x[valid]
    if x.size < 2:
        return 0.0
    nonincreasing = np.sum(np.diff(x) <= 0)
    return float(nonincreasing) / float(x.size - 1)


def finite_column_min(stacked: np.ndarray) -> np.ndarray:
    out = np.full(stacked.shape[1], np.nan, dtype=float)
    for i in range(stacked.shape[1]):
        col = stacked[:, i]
        finite = col[np.isfinite(col)]
        if finite.size:
            out[i] = float(np.min(finite))
    return out


def check_scaling(ncpu: np.ndarray, timestep: np.ndarray, label: str) -> None:
    frac = nonincreasing_fraction(timestep)
    if frac < 0.75:
        raise AssertionError(f"{label}: scaling is not mostly non-increasing (fraction={frac:.2f})")

    ref_idx = np.where(ncpu == 24)[0]
    if ref_idx.size == 0:
        raise AssertionError("Missing 24-CPU anchor point")
    iref = int(ref_idx[0])
    if not np.isfinite(timestep[iref]) or timestep[iref] <= 0:
        raise AssertionError(f"{label}: invalid 24-CPU timing")

    ideal = timestep[iref] * ncpu[iref] / ncpu
    valid = np.isfinite(timestep) & (timestep > 0)
    efficiency = timestep[valid] / ideal[valid]
    if np.min(efficiency) < 0.2 or np.max(efficiency) > 3.0:
        raise AssertionError(f"{label}: scaling deviates strongly from ideal envelope")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark-root",
        required=True,
        help="Path to implementation_2025/benchmark_performance_2d",
    )
    args = parser.parse_args()

    root = Path(args.benchmark_root)
    lb_trials = []
    lor_trials = []

    for trial in range(1, 6):
        trial_dir = root / f"trial_{trial}_standard_dynamical_chunk"
        lb = read_info(trial_dir / "info_lb.txt")
        lor = read_info(trial_dir / "info_lorentz.txt")
        lb_trials.append(lb[:, 1])
        lor_trials.append(lor[:, 1])

    lb_min = finite_column_min(np.vstack(lb_trials))
    lor_min = finite_column_min(np.vstack(lor_trials))

    check_scaling(EXPECTED_CPUS, lb_min, "Lorentz-Bath")
    check_scaling(EXPECTED_CPUS, lor_min, "Lorentz")

    valid = np.isfinite(lb_min) & np.isfinite(lor_min) & (lb_min > 0) & (lor_min > 0)
    if np.count_nonzero(valid) < 8:
        raise AssertionError("Insufficient finite timing points for robust overhead check")
    ratio = lb_min[valid] / lor_min[valid]
    median_ratio = float(np.median(ratio))
    if not (10.0 <= median_ratio <= 30.0):
        raise AssertionError(f"Median LB/Lorentz overhead {median_ratio:.2f} outside expected range [10, 30]")

    print("Figure-4 validation passed")
    print(f"Median LB/Lorentz overhead: {median_ratio:.3f}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
