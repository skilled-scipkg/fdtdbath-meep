#!/usr/bin/env python3
"""Validate figure-3 2D spectra against figure-2 1D references."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def load_array(path: Path) -> np.ndarray:
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 4:
        raise AssertionError(f"{path}: expected at least 4 columns")
    finite = np.all(np.isfinite(data[:, :4]), axis=1)
    filtered = data[finite]
    if filtered.shape[0] < 50:
        raise AssertionError(f"{path}: too few finite rows in first four columns")
    return filtered


def peak_freqs(freq: np.ndarray, trans: np.ndarray) -> tuple[float, float]:
    left = (freq >= 0.7) & (freq <= 1.0)
    right = (freq >= 1.0) & (freq <= 1.3)
    if not np.any(left) or not np.any(right):
        raise AssertionError("Missing LP/UP windows")
    lp = float(freq[left][np.argmax(trans[left])])
    up = float(freq[right][np.argmax(trans[right])])
    return lp, up


def peak_width(freq: np.ndarray, trans: np.ndarray, lo: float, hi: float) -> float:
    mask = (freq >= lo) & (freq <= hi)
    f = freq[mask]
    t = trans[mask]
    if f.size == 0:
        raise AssertionError(f"Missing data in window [{lo}, {hi}]")
    half = 0.5 * float(np.max(t))
    above = np.where(t >= half)[0]
    return float(f[above[-1]] - f[above[0]]) if above.size > 1 else 0.0


def validate_pair(one_d_file: Path, two_d_file: Path, tol: float = 0.01) -> None:
    d1 = load_array(one_d_file)
    d2 = load_array(two_d_file)

    f1, t1 = d1[:, 0], d1[:, 3]
    f2, t2 = d2[:, 0], d2[:, 3]

    lp1, up1 = peak_freqs(f1, t1)
    lp2, up2 = peak_freqs(f2, t2)

    if abs(lp1 - lp2) > tol:
        raise AssertionError(f"LP mismatch {one_d_file.name} vs {two_d_file.name}: {lp1:.4f} vs {lp2:.4f}")
    if abs(up1 - up2) > tol:
        raise AssertionError(f"UP mismatch {one_d_file.name} vs {two_d_file.name}: {up1:.4f} vs {up2:.4f}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--one-d-root", required=True, help="Path to implementation_2025/1d_harmonic_broadlinewidth")
    parser.add_argument("--two-d-root", required=True, help="Path to implementation_2025/2d_harmonic_broadlinewidth")
    args = parser.parse_args()

    one_d = Path(args.one_d_root)
    two_d = Path(args.two_d_root)

    pairs = [
        ("incav_lorentz_sigma_0.002/spectrum.txt", "incav_lorentz_sigma_0.002_angle_0/spectrum.txt"),
        ("incav_lorentz_sigma_0.020/spectrum.txt", "incav_lorentz_sigma_0.02_angle_0/spectrum.txt"),
        ("incav_lb_lorentzian_sigma_0.002/spectrum.txt", "incav_lb_lorentzian_sigma_0.002_angle_0/spectrum.txt"),
        ("incav_lb_lorentzian_sigma_0.020/spectrum.txt", "incav_lb_lorentzian_sigma_0.02_angle_0/spectrum.txt"),
    ]

    for one_d_rel, two_d_rel in pairs:
        validate_pair(one_d / one_d_rel, two_d / two_d_rel)

    l020 = load_array(two_d / "incav_lorentz_sigma_0.02_angle_0/spectrum.txt")
    lb020 = load_array(two_d / "incav_lb_lorentzian_sigma_0.02_angle_0/spectrum.txt")

    for lo, hi in [(0.7, 1.0), (1.0, 1.3)]:
        width_l = peak_width(l020[:, 0], l020[:, 3], lo, hi)
        width_lb = peak_width(lb020[:, 0], lb020[:, 3], lo, hi)
        peak_l = float(np.max(l020[(l020[:, 0] >= lo) & (l020[:, 0] <= hi), 3]))
        peak_lb = float(np.max(lb020[(lb020[:, 0] >= lo) & (lb020[:, 0] <= hi), 3]))
        if not (width_lb < width_l):
            raise AssertionError(
                f"Expected narrower LB 2D peak in [{lo}, {hi}], got LB={width_lb:.6f}, Lorentz={width_l:.6f}"
            )
        if not (peak_lb > peak_l):
            raise AssertionError(
                f"Expected stronger LB 2D peak in [{lo}, {hi}], got LB={peak_lb:.6f}, Lorentz={peak_l:.6f}"
            )

    print("Figure-3 validation passed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
