#!/usr/bin/env python3
"""Extract and validate figure-1 geometry constants from staged run scripts."""

from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from pathlib import Path

EXPECTED = {
    "mirror_index": 10.0,
    "mirror_thickness_um": 0.02,
    "slab_thickness_um": 1.0,
    "incident_angle_default_deg": 0.0,
}


def parse_default_refractive_index(text: str) -> float:
    pattern = re.compile(r"--refractiveindex'.*?default=([0-9.]+)", re.DOTALL)
    match = pattern.search(text)
    if not match:
        raise ValueError("Could not parse default refractive index")
    return float(match.group(1))


def parse_layer_thicknesses(text: str) -> list[float]:
    match = re.search(r"layer_thicknesses\s*=\s*np\.array\(\[([^\]]+)\]\)", text)
    if not match:
        raise ValueError("Could not parse layer_thicknesses")
    expr = "[" + match.group(1) + "]"
    node = ast.parse(expr, mode="eval")
    compiled = compile(node, "<layer_thicknesses>", "eval")
    values = eval(compiled, {"__builtins__": {}}, {"t_wavelength": 1.0})
    return [float(v) for v in values]


def parse_layer_indexes(text: str, n1: float) -> list[float]:
    match = re.search(r"layer_indexes\s*=\s*np\.array\(\[([^\]]+)\]\)", text)
    if not match:
        raise ValueError("Could not parse layer_indexes")
    expr = "[" + match.group(1) + "]"
    node = ast.parse(expr, mode="eval")
    compiled = compile(node, "<layer_indexes>", "eval")
    values = eval(compiled, {"__builtins__": {}}, {"n1": n1})
    return [float(v) for v in values]


def parse_pml(text: str) -> float:
    match = re.search(r"pml_thickness\s*=\s*([0-9.]+)", text)
    if not match:
        raise ValueError("Could not parse pml_thickness")
    return float(match.group(1))


def parse_incident_angle_default(text: str) -> float:
    match = re.search(r"--incidentangle'.*?default=([0-9.]+)", text, re.DOTALL)
    if not match:
        raise ValueError("Could not parse default incident angle")
    return float(match.group(1))


def parse_model_signatures(text: str) -> dict[str, bool]:
    return {
        "lorentzian_present": "mp.LorentzianSusceptibility(" in text,
        "bath_lorentzian_present": "mp.BathLorentzianSusceptibility(" in text,
    }


def parse_boundary_signatures(text: str, script_kind: str) -> dict[str, bool]:
    out = {
        "x_directed_pml_present": "direction=mp.X" in text,
    }
    if script_kind == "2d":
        out["bloch_k_point_present"] = "k_point" in text
        out["mirror_y_symmetry_present"] = "mp.Mirror(mp.Y)" in text
    return out


def check_close(name: str, value: float, expected: float, tol: float = 1e-9) -> None:
    if abs(value - expected) > tol:
        raise AssertionError(f"{name}={value} does not match expected {expected}")


def parse_script(path: Path, script_kind: str) -> dict:
    text = path.read_text(encoding="utf-8")
    default_n1 = parse_default_refractive_index(text)
    layer_thicknesses = parse_layer_thicknesses(text)
    layer_indexes = parse_layer_indexes(text, n1=default_n1)

    payload = {
        "file": str(path),
        "default_refractive_index": default_n1,
        "layer_thicknesses_um": layer_thicknesses,
        "layer_indexes": layer_indexes,
        "mirror_thickness_um": layer_thicknesses[1],
        "slab_thickness_um": layer_thicknesses[2],
        "pml_thickness_um": parse_pml(text),
        "model_signatures": parse_model_signatures(text),
        "boundary_signatures": parse_boundary_signatures(text, script_kind=script_kind),
    }
    if script_kind == "2d":
        payload["incident_angle_default_deg"] = parse_incident_angle_default(text)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--script-1d", required=True)
    parser.add_argument("--script-2d", required=True)
    parser.add_argument("--output", default="figure1_geometry_constants.json")
    args = parser.parse_args()

    one_d = parse_script(Path(args.script_1d), script_kind="1d")
    two_d = parse_script(Path(args.script_2d), script_kind="2d")

    check_close("1d mirror index", one_d["default_refractive_index"], EXPECTED["mirror_index"])
    check_close("2d mirror index", two_d["default_refractive_index"], EXPECTED["mirror_index"])
    check_close("1d mirror thickness", one_d["mirror_thickness_um"], EXPECTED["mirror_thickness_um"])
    check_close("2d mirror thickness", two_d["mirror_thickness_um"], EXPECTED["mirror_thickness_um"])
    check_close("1d slab thickness", one_d["slab_thickness_um"], EXPECTED["slab_thickness_um"])
    check_close("2d slab thickness", two_d["slab_thickness_um"], EXPECTED["slab_thickness_um"])
    check_close(
        "2d default incident angle",
        two_d["incident_angle_default_deg"],
        EXPECTED["incident_angle_default_deg"],
    )

    if not one_d["model_signatures"]["lorentzian_present"] or not one_d["model_signatures"]["bath_lorentzian_present"]:
        raise AssertionError("1d driver must include both Lorentz and Lorentz-Bath model definitions")
    if not two_d["model_signatures"]["lorentzian_present"] or not two_d["model_signatures"]["bath_lorentzian_present"]:
        raise AssertionError("2d driver must include both Lorentz and Lorentz-Bath model definitions")

    if not one_d["boundary_signatures"]["x_directed_pml_present"]:
        raise AssertionError("1d driver must include x-directed PML setup")
    if not two_d["boundary_signatures"]["x_directed_pml_present"]:
        raise AssertionError("2d driver must include x-directed PML setup")
    if not two_d["boundary_signatures"]["bloch_k_point_present"]:
        raise AssertionError("2d driver must include Bloch k-point handling")

    payload = {
        "validated": True,
        "expected": EXPECTED,
        "script_1d": one_d,
        "script_2d": two_d,
        "notes": [
            "Figure-1 schematic image fdtd_bath_demo.png is not present in the provided manifests.",
            "Use this geometry report as the authoritative figure-1 parameter reference.",
        ],
    }

    out = Path(args.output)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
