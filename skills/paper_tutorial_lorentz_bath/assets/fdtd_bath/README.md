# FDTD with Auxiliary Bath Fields (FDTD-Bath) for Condensed-Phase Polaritonics

This Github repo stores the input and post-processing files for the publications involving the FDTD-Bath approach developed in the [TEL research group](https://www.taoeli.org/).

## Installation

The FDTD-Bath approach is implemented on top of a modified MEEP codebase from the [TEL research group](https://www.taoeli.org/). This functionality has not been merged into upstream MEEP.

Install the packaged build from Conda channels (`tel-research` + `conda-forge`) with MPICH-enabled MPI:

Create a new Conda environment:

```bash
conda create -n fdtdbath-run --override-channels \
  -c tel-research -c conda-forge \
  "pymeep-fdtdbath * mpi_mpich_*"
```

Or install into an existing Conda environment:

```bash
conda install --override-channels \
  -c tel-research -c conda-forge \
  "pymeep-fdtdbath * mpi_mpich_*"
```

The package name is `pymeep-fdtdbath`; Python import remains `import meep as mp`.

Legacy source-install helper scripts are retained in [./installation_scripts/](./installation_scripts/) for environments where this packaged installation cannot be used:

- [./installation_scripts/meep_install_CentOS9.sh](./installation_scripts/meep_install_CentOS9.sh): tested on a clean CentOS9 system with **sudo** privileges.
- [./installation_scripts/meep_install_hpc_anvil.sh](./installation_scripts/meep_install_hpc_anvil.sh): tested on the Anvil HPC system, available through the NSF ACCESS program (https://allocations.access-ci.org/get-your-first-project).

## Data for FDTD-Bath publications

- [implementation_2025](./implementation_2025/)

This folder contains the input and postprocessing files for the initial implementation of the FDTD-Bath approach:

Tao E. Li. "FDTD with Auxiliary Bath Fields for Condensed-Phase Polaritonics: Fundamentals and Implementation". [arXiv:2505.23963](https://arxiv.org/abs/2505.23963) (2025)
