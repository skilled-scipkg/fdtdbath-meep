# fig_001: Lorentz vs Lorentz-Bath model and cavity setup geometry

## Scientific aim
Pin down the exact 1D/2D Fabry-Perot geometry and boundary assumptions reused by Figures 2-4, and verify that the simulation drivers expose both conventional Lorentz and Lorentz-Bath model pathways.

## Runtime recipe
1. Create an execution folder under `projects/YYYY-MM-DD-<scope>/` and stage assets if needed:
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

2. Extract geometry and model-signature metadata:
```bash
cd "$RUN_DIR"
python scripts/extract_fig1_geometry.py \
  --script-1d fdtd_bath/implementation_2025/1d_harmonic_broadlinewidth/run_sc_fabry_perot_1d.py \
  --script-2d fdtd_bath/implementation_2025/2d_harmonic_broadlinewidth/run_sc_fabry_perot_2d_spectrum.py \
  --output fig_001_geometry_constants.json
```

## Validation checklist
- `fig_001_geometry_constants.json` exists and reports `validated: true`.
- Geometry constants match manuscript values: mirror index `10.0`, mirror thickness `0.02 um`, slab thickness `1.0 um` in both 1D and 2D drivers.
- Boundary-condition controls are present: 1D `pml_thickness=0.5`; 2D `pml_thickness=2.0`, `incidentangle` argument with default `0.0`, and x-directed PML.
- Model distinction is explicit in run drivers: both `mp.LorentzianSusceptibility` and `mp.BathLorentzianSusceptibility` are detected.
- Known asset gap is documented: `fdtd_bath_demo.png` is not packaged, so the reproducible figure-1 artifact is the validated geometry/model metadata file.
