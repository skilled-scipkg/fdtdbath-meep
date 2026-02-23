# Simulation data for the initial implementation of the FDTD-Bath approach:

Tao E. Li. "FDTD with Auxiliary Bath Fields for Condensed-Phase Polaritonics: Fundamentals and Implementation". [arXiv:2505.23963](https://arxiv.org/abs/2505.23963) (2025)

- [**plotting/**](./plotting/): run each python script for plotting the corresponding data figures in the manuscript.

- [**1d_harmonic_broadlinewidth**](./1d_harmonic_broadlinewidth/): 1D simulation data. Users can simply run **./submit_*.sh** locally to perform the 1D simulations. Each 1D simulation takes approximately 1 minute to finish. The same submission bash scripts can also be used for slurm job submission in HPC.

- [**2d_harmonic_broadlinewidth**](./2d_harmonic_broadlinewidth/): 2D simulation data. The **./submit_*.sh** bash scripts require MPI parallel calculations with 48 CPUs for slurm job submission. Users can also slightly modify the bash scripts to reduce the required number of CPUs to run the simulation locally.

- [**benchmark_performance_2d**](./benchmark_performance_2d/): 2D simulation data for examining the MPI performance across multiple nodes.
