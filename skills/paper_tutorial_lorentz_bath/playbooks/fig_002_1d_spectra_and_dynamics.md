# fig_002: 1D free-space and in-cavity spectra plus energy dynamics

## Scientific aim
Reproduce 1D evidence that Lorentz-Bath(L) preserves free-space linewidth matching while producing narrower/stronger polariton peaks and longer-lived cavity energy oscillations at stronger coupling (`sigma=0.020`).

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

2. Run all 1D simulation batches:
```bash
cd "$RUN_DIR/fdtd_bath/implementation_2025/1d_harmonic_broadlinewidth"
bash submit_outcav.sh
bash submit_incav_lorentz.sh
bash submit_incav_lb_lorentzian.sh
```

3. Generate the publication panel:
```bash
cd ../plotting
python plot_1d_demo.py
```

4. Run figure-specific acceptance checks:
```bash
cd "$RUN_DIR"
python scripts/validate_figure2_outputs.py \
  --root fdtd_bath/implementation_2025/1d_harmonic_broadlinewidth
```

## Validation checklist
- All required `spectrum.txt` outputs exist, each with 4 finite columns (`freq`, `loss`, `R`, `T`).
- All required `energy_dynamics.txt` outputs exist, each with 2 finite columns (`time`, cavity EM energy).
- Figure 2a condition holds: free-space Lorentz and Lorentz-Bath(U) overlap in peak position/linewidth, while Lorentz-Bath(L) shows reduced tails.
- Figure 2b condition holds: at `sigma=0.020`, Lorentz-Bath(L) LP/UP peaks are narrower and stronger than Lorentz; at `sigma=0.002` the two models remain near-overlapping.
- Figure 2c condition holds: at `sigma=0.020`, Lorentz-Bath(L) retains larger late-time cavity EM energy than Lorentz.
