# fig_003: 2D normal-incidence cavity transmission spectra

## Scientific aim
Confirm dimensional robustness by reproducing 2D normal-incidence spectra and quantitatively matching 1D LP/UP peak positions within `0.01 um^-1`.

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

2. Run 2D sigma scans (submission scripts use `mpirun -np 48`):
```bash
cd "$RUN_DIR/fdtd_bath/implementation_2025/2d_harmonic_broadlinewidth"
bash submit_incav_lorentz_scan_sigma.sh
bash submit_incav_lb_lorentzian_scan_sigma.sh
```

3. Generate the publication panel:
```bash
cd ../plotting
python plot_2d_demo.py
```

4. Validate 2D-vs-1D consistency:
```bash
cd "$RUN_DIR"
python scripts/validate_figure3_vs_figure2.py \
  --one-d-root fdtd_bath/implementation_2025/1d_harmonic_broadlinewidth \
  --two-d-root fdtd_bath/implementation_2025/2d_harmonic_broadlinewidth
```

## Validation checklist
- Required 2D `spectrum.txt` files exist with finite data in the first four columns.
- LP and UP peak frequencies agree with the corresponding 1D runs within `0.01 um^-1` for all model/sigma pairings.
- At `sigma=0.02`, Lorentz-Bath(L) remains narrower and stronger than Lorentz in both LP and UP windows.
- Output figure `implementation_2025/plotting/2d_spectrum_demo.pdf` is regenerated from current run outputs.
