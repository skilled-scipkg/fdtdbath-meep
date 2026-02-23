# FDTD with Auxiliary Bath Fields (FDTD-Bath) for Condensed-Phase Polaritonics

This Github repo stores the input and post-processing files for the publications involving the FDTD-Bath approach developed in the [TEL research group](https://www.taoeli.org/).

## Installation

The FDTD-Bath approach is implemented on top of the open-source FDTD package [MEEP](https://meep.readthedocs.io/en/master/). However, at this moment, this approach has not been merged into to the main version of MEEP. Users need to install [our modified MEEP code](https://github.com/TaoELi/meep) from source: https://github.com/TaoELi/meep. The MEEP offical website contains the detailed discussion of [installing MEEP from source](https://meep.readthedocs.io/en/latest/Build_From_Source/).

Because it is very tedious to install the MEEP package from source, the following bash scripts are provided in [./installation_scripts/](./installation_scripts/) for smoothly installing the MEEP package in a few different Linux environments. 

- [./installation_scripts/meep_install_CentOS9.sh](./installation_scripts/meep_install_CentOS9.sh): Tested on a clean CentOS9 system with **sudo** privileges. 

- [./installation_scripts/meep_install_hpc_anvil.sh](./installation_scripts/meep_install_hpc_anvil.sh): Tested on the Anvil HPC system, which is available to U.S. researchers through the NSF ACCESS program; see https://allocations.access-ci.org/get-your-first-project. 

## Data for FDTD-Bath publications

- [implementation_2025](./implementation_2025/)

This folder contains the input and postprocessing files for the initial implementation of the FDTD-Bath approach:

Tao E. Li. "FDTD with Auxiliary Bath Fields for Condensed-Phase Polaritonics: Fundamentals and Implementation". [arXiv:2505.23963](https://arxiv.org/abs/2505.23963) (2025)
