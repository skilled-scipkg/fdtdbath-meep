#!/usr/bin/env python3
"""Validate figure-2 data integrity and key Lorentz vs Lorentz-Bath trends."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

SPECTRUM_FILES = [
    "outcav_lorentz/spectrum.txt",
    "outcav_lb_uniform/spectrum.txt",
    "outcav_lb_lorentzian/spectrum.txt",
    "incav_lorentz_sigma_0.002/spectrum.txt",
    "incav_lorentz_sigma_0.020/spectrum.txt",
    "incav_lb_lorentzian_sigma_0.002/spectrum.txt",
    "incav_lb_lorentzian_sigma_0.020/spectrum.txt",
]

ENERGY_FILES = [
    "incav_lorentz_sigma_0.002/energy_dynamics.txt",
    "incav_lorentz_sigma_0.020/energy_dynamics.txt",
    "incav_lb_lorentzian_sigma_0.002/energy_dynamics.txt",
    "incav_lb_lorentzian_sigma_0.020/energy_dynamics.txt",
]

FREE_SPACE_WINDOW = (0.8, 1.2)
SHOULDER_OFFSETS = (0.08, 0.20)


def load_array(path: Path, expected_cols: int) -> np.ndarray:
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != expected_cols:
        raise AssertionError(f"{path}: expected {expected_cols} columns, found {data.shape[1]}")
    if not np.all(np.isfinite(data)):
        raise AssertionError(f"{path}: found non-finite values")
    return data


def peak_metrics(freq: np.ndarray, trans: np.ndarray, lo: float, hi: float) -> tuple[float, float, float]:
    mask = (freq >= lo) & (freq <= hi)
    if not np.any(mask):
        raise AssertionError(f"No data points in window [{lo}, {hi}]")
    f = freq[mask]
    t = trans[mask]
    idx = int(np.argmax(t))
    peak_val = float(t[idx])
    peak_freq = float(f[idx])

    half = 0.5 * peak_val
    above = np.where(t >= half)[0]
    width = float(f[above[-1]] - f[above[0]]) if above.size > 1 else 0.0
    return peak_val, width, peak_freq


def absorption_peak_metrics(freq: np.ndarray, trans: np.ndarray, lo: float, hi: float) -> tuple[float, float, float]:
    mask = (freq >= lo) & (freq <= hi)
    if not np.any(mask):
        raise AssertionError(f"No data points in absorption window [{lo}, {hi}]")
    f = freq[mask]
    a = 1.0 - trans[mask]
    idx = int(np.argmax(a))
    peak_val = float(a[idx])
    peak_freq = float(f[idx])
    half = 0.5 * peak_val
    above = np.where(a >= half)[0]
    width = float(f[above[-1]] - f[above[0]]) if above.size > 1 else 0.0
    return peak_val, width, peak_freq


def persistence_time(t: np.ndarray, e: np.ndarray, threshold: float = 0.2) -> float:
    emax = float(np.max(e))
    if emax <= 0:
        raise AssertionError("Energy trace has non-positive maximum")
    norm = e / emax
    i_peak = int(np.argmax(norm))
    indices = np.where(norm[i_peak:] <= threshold)[0]
    if indices.size == 0:
        return float(t[-1])
    return float(t[i_peak + indices[0]])


def normalized_tail_mean(t: np.ndarray, e: np.ndarray, tmin: float = 20.0) -> float:
    emax = float(np.max(e))
    if emax <= 0:
        raise AssertionError("Energy trace has non-positive maximum")
    norm = e / emax
    mask = t >= tmin
    if not np.any(mask):
        raise AssertionError(f"No points found in tail window t >= {tmin}")
    return float(np.mean(norm[mask]))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, help="Path to implementation_2025/1d_harmonic_broadlinewidth")
    args = parser.parse_args()

    root = Path(args.root)
    for rel in SPECTRUM_FILES:
        path = root / rel
        if not path.exists():
            raise AssertionError(f"Missing required spectrum file: {path}")
        load_array(path, expected_cols=4)

    for rel in ENERGY_FILES:
        path = root / rel
        if not path.exists():
            raise AssertionError(f"Missing required energy file: {path}")
        load_array(path, expected_cols=2)

    out_lor = load_array(root / "outcav_lorentz/spectrum.txt", 4)
    out_lb_uniform = load_array(root / "outcav_lb_uniform/spectrum.txt", 4)
    out_lb_lorentzian = load_array(root / "outcav_lb_lorentzian/spectrum.txt", 4)

    lo, hi = FREE_SPACE_WINDOW
    f_ref = out_lor[:, 0]
    win = (f_ref >= lo) & (f_ref <= hi)
    if np.count_nonzero(win) < 200:
        raise AssertionError(f"Too few free-space points in [{lo}, {hi}]")
    f_win = f_ref[win]
    t_lor = out_lor[win, 3]
    t_lb_uniform = np.interp(f_win, out_lb_uniform[:, 0], out_lb_uniform[:, 3])
    t_lb_lorentzian = np.interp(f_win, out_lb_lorentzian[:, 0], out_lb_lorentzian[:, 3])

    rmse_uniform = float(np.sqrt(np.mean((t_lor - t_lb_uniform) ** 2)))
    if rmse_uniform > 0.03:
        raise AssertionError(
            f"Expected Lorentz and LB-uniform free-space overlap; RMSE={rmse_uniform:.4f} exceeds 0.03"
        )

    apeak_l, awidth_l, afreq_l = absorption_peak_metrics(f_win, t_lor, lo, hi)
    apeak_u, awidth_u, afreq_u = absorption_peak_metrics(f_win, t_lb_uniform, lo, hi)
    if abs(afreq_l - afreq_u) > 0.01:
        raise AssertionError(
            f"Lorentz and LB-uniform absorption centers differ too much: {afreq_l:.4f} vs {afreq_u:.4f}"
        )
    width_rel_err = abs(awidth_l - awidth_u) / max(1e-12, awidth_l)
    if width_rel_err > 0.35:
        raise AssertionError(
            f"Lorentz and LB-uniform linewidth mismatch too large: rel_err={width_rel_err:.3f}"
        )

    inner, outer = SHOULDER_OFFSETS
    shoulder = (np.abs(f_win - afreq_l) >= inner) & (np.abs(f_win - afreq_l) <= outer)
    if np.count_nonzero(shoulder) < 60:
        raise AssertionError("Insufficient shoulder points for Lorentz-Bath(L) tail check")
    shoulder_abs_l = float(np.mean(1.0 - t_lor[shoulder]))
    shoulder_abs_lbl = float(np.mean(1.0 - t_lb_lorentzian[shoulder]))
    if shoulder_abs_l <= 1e-6:
        raise AssertionError("Lorentz shoulder absorption is too small to evaluate tail suppression")
    if not (shoulder_abs_lbl < 0.70 * shoulder_abs_l):
        raise AssertionError(
            "Expected reduced free-space tails for LB-lorentzian; "
            f"shoulder absorptions LB-L={shoulder_abs_lbl:.6f}, Lorentz={shoulder_abs_l:.6f}"
        )

    l002 = load_array(root / "incav_lorentz_sigma_0.002/spectrum.txt", 4)
    lb002 = load_array(root / "incav_lb_lorentzian_sigma_0.002/spectrum.txt", 4)
    l020 = load_array(root / "incav_lorentz_sigma_0.020/spectrum.txt", 4)
    lb020 = load_array(root / "incav_lb_lorentzian_sigma_0.020/spectrum.txt", 4)

    for lo, hi in [(0.7, 1.0), (1.0, 1.3)]:
        p_l020, w_l020, _ = peak_metrics(l020[:, 0], l020[:, 3], lo, hi)
        p_lb020, w_lb020, _ = peak_metrics(lb020[:, 0], lb020[:, 3], lo, hi)
        if not (w_lb020 < w_l020):
            raise AssertionError(
                f"Expected narrower LB peak in window [{lo}, {hi}], got widths LB={w_lb020:.6f}, Lorentz={w_l020:.6f}"
            )
        if not (p_lb020 > p_l020):
            raise AssertionError(
                f"Expected stronger LB peak in window [{lo}, {hi}], got peaks LB={p_lb020:.6f}, Lorentz={p_l020:.6f}"
            )

    # Weak-coupling traces should stay close.
    p_l002, _, _ = peak_metrics(l002[:, 0], l002[:, 3], 0.7, 1.0)
    p_lb002, _, _ = peak_metrics(lb002[:, 0], lb002[:, 3], 0.7, 1.0)
    rel_diff = abs(p_lb002 - p_l002) / max(1e-12, p_l002)
    if rel_diff > 0.25:
        raise AssertionError(f"Sigma=0.002 Lorentz/LB peak mismatch too large ({rel_diff:.3f})")

    e_l002 = load_array(root / "incav_lorentz_sigma_0.002/energy_dynamics.txt", 2)
    e_lb002 = load_array(root / "incav_lb_lorentzian_sigma_0.002/energy_dynamics.txt", 2)
    e_l020 = load_array(root / "incav_lorentz_sigma_0.020/energy_dynamics.txt", 2)
    e_lb020 = load_array(root / "incav_lb_lorentzian_sigma_0.020/energy_dynamics.txt", 2)

    # Crossing-time metrics can be noisy under oscillatory decay; use normalized late-time
    # mean energy as the primary persistence metric.
    tail_l002 = normalized_tail_mean(e_l002[:, 0], e_l002[:, 1])
    tail_lb002 = normalized_tail_mean(e_lb002[:, 0], e_lb002[:, 1])
    tail_l020 = normalized_tail_mean(e_l020[:, 0], e_l020[:, 1])
    tail_lb020 = normalized_tail_mean(e_lb020[:, 0], e_lb020[:, 1])

    if not (tail_lb020 > 2.0 * tail_l020):
        raise AssertionError(
            f"Expected larger late-time LB energy at sigma=0.020, got LB={tail_lb020:.6f}, Lorentz={tail_l020:.6f}"
        )

    near_overlap_ratio = tail_lb002 / max(1e-12, tail_l002)
    if not (0.50 <= near_overlap_ratio <= 1.50):
        raise AssertionError(
            f"Expected near-overlap at sigma=0.002; tail-energy ratio={near_overlap_ratio:.3f}"
        )

    print("Figure-2 validation passed")
    print(
        "free-space checks: "
        f"uniform RMSE={rmse_uniform:.5f}, "
        f"absorption-center delta={abs(afreq_l - afreq_u):.5f}, "
        f"tail ratio LB-L/Lorentz={shoulder_abs_lbl / shoulder_abs_l:.3f}"
    )
    print(f"sigma=0.020 normalized tail mean (Lorentz, LB): {tail_l020:.6f}, {tail_lb020:.6f}")
    print(f"sigma=0.002 normalized tail ratio LB/Lorentz: {near_overlap_ratio:.3f}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
