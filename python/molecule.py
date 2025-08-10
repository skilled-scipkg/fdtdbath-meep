'''
This Python module defines realistic molecules coupled with the MEEP FDTD engine.

The important issue is that we plan to propagate molecules under atomic units, while MEEP uses its own unit system.
'''

import warnings
import numpy as np
import meep as mp
import math
from scipy.linalg import expm
from math import exp

class DummyMolecule():
    """
    A dummy molecule class that serves as a parent class for realistic molecules.
    """
    
    def __init__(self, dt, dx, center=mp.Vector3(), size=mp.Vector3(1, 1, 1)):
        '''
        Initialize the dummy molecule with default properties.
        This class is intended to be subclassed for specific molecular implementations.
        '''
        self.dt = dt  # time step for the simulation used in MEEP
        self.dx = dx  # spatial resolution for the simulation used in MEEP
        self.center = center
        self.size = size
        # each molecule has a corresponding MEEP Source object and will be used in defining the MEEP simulation object
        self.sources = []

    def _propagate(self, int_ep):
        """
        Propagate the quantum molecular dynamics given the light-matter coupling int_ep.
        This method should be overridden by subclasses to implement specific propagation logic.
        """
        raise NotImplementedError("This method should be overridden by subclasses.")

    def _init_sources(self):
        """
        Initialize the sources for this molecule.
        This method should be overridden by subclasses to implement specific source initialization logic.
        """
        raise NotImplementedError("This method should be overridden by subclasses.")

    def _update_source_amplitude(self):
        """
        Update the source amplitude after propagating this molecule for one time step.
        This method should be overridden by subclasses to implement specific source update logic.
        """
        raise NotImplementedError("This method should be overridden by subclasses.")
    
    def _calculate_ep_integral(self, sim):
        """
        Calculate the integral of the electric field over the molecule's volume.
        This method should be overridden by subclasses to implement specific integral calculation logic.
        """
        raise NotImplementedError("This method should be overridden by subclasses.")

    def update_molecule_with_em(self, sources_non_molecule=[]):
        def __step_function__(sim):
            """
            A step function which is used in meep simulation.run() to account for the 
            interaction between this molecule and the Maxwell field at each time step.

            - sources_non_molecule: A list of sources that are not part of this molecule.
            Typically any non-molecule sources such as Gaussian sources must be included here.
            """
            # 1. calculate the light-matter coupling: int_ep
            int_ep = self._calculate_ep_integral(sim)

            # 2. propagate the quantum molecular dynamics given int_ep
            self._propagate(int_ep)

            # 3. update the source amplitude after propagating this molecule for one time step
            self._update_source_amplitude()
            sim.change_sources(sources_non_molecule + self.sources)

        return __step_function__
    
    def units_helper(self, time_units_in_fs=0.1):
        """
        Helper function to explain the unit system used in MEEP and its connection to atomic units.
        """
        print("MEEP uses its own unit system, which is based on the speed of light in vacuum (c), ",
                "the permittivity of free space (epsilon_0), and the permeability of free space (mu_0). ",
                "To couple MEEP with molecular dynamics, we set [c] = [epsilon_0] = [mu_0] = [hbar] = 1. ",
                "By further defining the time unit as %.2E fs, we can convert the MEEP unit (mu) to atomic units (au) as follows:" % time_units_in_fs,
                "- angular frequency [omega]: 1 mu = %.4E eV" %(41.357 * 0.1 / time_units_in_fs),
                "- length: 1 mu = %.3E nm" %(4.771 * time_units_in_fs / 0.1),
                "- dipole moment: 1 mu = %.4E Debye" %(1.8970 * time_units_in_fs / 0.1))
        
    def print_dynamics(self):
        """
        Print the molecular dynamics properties.
        This method can be overridden by subclasses to provide specific molecular dynamics information.
        """
        raise NotImplementedError("This method should be overridden by subclasses to print specific molecular dynamics properties.")


def update_multiple_molecules_with_em(sources_non_molecule=[], molecules=[]):
    if sources_non_molecule is None:
        sources_non_molecule = []
    if molecules is None:
        molecules = []
    def __step_function__(sim):
        """
            A step function which is used in meep simulation.run() to account for the 
            interaction between multiple molecules and the Maxwell field at each time step.

            - sources_non_molecule: A list of sources that are not part of this molecule.
            Typically any non-molecule sources such as Gaussian sources must be included here.
        """
        sources = list(sources_non_molecule)
        for molecule in molecules:
            int_ep = molecule._calculate_ep_integral(sim)
            molecule._propagate(int_ep)
            molecule._update_source_amplitude()
            sources.extend(molecule.sources)
        sim.change_sources(sources)

    return __step_function__


class TLSMolecule(DummyMolecule):
    """
    A simple two-level system (TLS) molecule class that inherits from DummyMolecule.
    This class represents a basic model of a molecule with two energy levels.
    """

    def __init__(self, resolution, center=mp.Vector3(), size=mp.Vector3(1, 1, 1), frequency=1.0, dipole_moment=0.1, sigma=1.0, orientation=mp.Ez, dimensions=2):
        """
        Initialize the TLS molecule with an energy gap and dipole moment.
        
        Parameters:
        - resolution: The resolution of the MEEP simulation, which determines the time step and spatial resolution.
        - center: The center of the molecule in the simulation space (default is mp.Vector3(0, 0, 0)).
        - size: The size of the molecule in the simulation space (default is mp.Vector3(1, 1, 1)).
        - frequency: The frequency of the molecule [omega = 2.0 * pi * frequency].
        - dipole_moment: The transient dipole moment of the molecule.
        - sigma: The width of the Gaussian profile for the TLS polarization density.
        - orientation: The polarization orientation of the molecule (default is mp.Ex).
        - dimensions: The dimensionality of the system (default is 3 for 3D, can be set to 2 for 2D).
        """
        dt = 0.5 / resolution  # time step for the simulation
        dx = 1.0 / resolution  # spatial resolution for the simulation
        super().__init__(dt=dt, dx=dx, center=center, size=size)
        self.omega = frequency #* 2.0 * np.pi  # convert frequency to angular frequency
        self.dipole_moment = dipole_moment
        self.sigma = sigma
        self.orientation = orientation
        self.dimensions = dimensions
        # for 2D simulations, the TLS orientation is fixed to mp.Ez
        if self.dimensions == 2 and self.orientation != mp.Ez:
            warnings.warn("In 2D simulations, the TLS orientation is fixed to mp.Ez. Setting orientation to mp.Ez.")
            self.orientation = mp.Ez
        if self.dimensions == 1:
            warnings.warn("1D simulations are not supported for TLSMolecule. Please use 2D or 3D simulations.")
            raise ValueError("1D simulations are not supported for TLSMolecule. Please use 2D or 3D simulations.")
        # the polarization prefactor in 3D
        self._polarization_prefactor_3d = 1.0 / (2.0 * np.pi)**(1.5) / self.sigma**5 * self.dipole_moment
        # the polarization prefactor in 2D
        self._polarization_prefactor_2d = 1.0 / (2.0 * np.pi)**(1.0) / self.sigma**2 * self.dipole_moment

        # set the Hamiltonian and density matrix for the TLS molecule
        self.Hs = np.matrix([[0, 0],[0, self.omega]], dtype=np.complex128)
        self.SIGMAX = np.matrix([[0,1],[1,0]], dtype=np.complex128)
        self.expHs = expm(-1j * self.dt * self.Hs / 2.0) 
        self.C = np.matrix([[1],[0]], dtype=np.complex128)
        self.rho = np.dot(self.C, self.C.conj().transpose() )
        self.rho_history = []
        # initialize the sources for the TLS molecule
        self._init_sources()
    
    def reset_tls_population(self, excited_population=0.1):
        """
        Reset the TLS population to a specified excited state population.
        
        Parameters:
        - excited_population: The population of the excited state (default is 0.1).
        """
        if excited_population < 0 or excited_population > 1:
            raise ValueError("Excited population must be between 0 and 1.")
        self.C = np.matrix([[(1 - excited_population)**0.5], [excited_population**0.5]], dtype=np.complex128)
        self.rho = np.dot(self.C, self.C.conj().transpose())

    def _init_sources(self):
        """ 
        Initialize the sources for the TLS molecule.
        This method sets up the polarization profile based on the molecule's properties.
        """
        # set the source corresponding to the TLS molecule
        local_size = 10.0
        if self.size[0] > local_size or self.size[1] > local_size or self.size[2] > local_size:
            raise ValueError("The size of the molecule is too large. Please use a smaller size.")
        x = np.arange(-local_size/2.0, local_size/2.0, self.dx)
        y = np.arange(-local_size/2.0, local_size/2.0, self.dx)
        z = np.arange(-local_size/2.0, local_size/2.0, self.dx)
        if self.dimensions == 2:
            [X, Y] = np.meshgrid(x, y)
            R2 = X**2 + Y**2
            Pz_raw = self._polarization_prefactor_2d * np.exp(-R2 / 2.0 / self.sigma**2)
            start_idx_x, end_idx_x = int((local_size - self.size[0])/2/self.dx) , int((local_size + self.size[0])/2/self.dx)
            start_idx_y, end_idx_y = int((local_size - self.size[1])/2/self.dx) , int((local_size + self.size[1])/2/self.dx)
            Pz0 = Pz_raw[start_idx_x:end_idx_x, start_idx_y:end_idx_y]
            Pz0 = np.reshape(Pz0, (math.ceil(self.size[0]/self.dx), math.ceil(self.size[1]/self.dx), 1))
            polarization_profile_2d = np.copy(Pz0.astype(np.complex128), order='C')

        elif self.dimensions == 3:
            [X, Y, Z] = np.meshgrid(x, y, z)
            R2 = X**2 + Y**2 + Z**2
            if self.orientation == mp.Ez:
                P3d_raw = self._polarization_prefactor_3d * Z**2 * np.exp(-R2 / 2.0 / self.sigma**2)
            elif self.orientation == mp.Ex:
                P3d_raw = self._polarization_prefactor_3d * X**2 * np.exp(-R2 / 2.0 / self.sigma**2)
            elif self.orientation == mp.Ey:
                P3d_raw = self._polarization_prefactor_3d * Y**2 * np.exp(-R2 / 2.0 / self.sigma**2)
            start_idx_x, end_idx_x = int((local_size - self.size[0])/2/self.dx) , int((local_size + self.size[0])/2/self.dx)
            start_idx_y, end_idx_y = int((local_size - self.size[1])/2/self.dx) , int((local_size + self.size[1])/2/self.dx)
            start_idx_z, end_idx_z = int((local_size - self.size[2])/2/self.dx) , int((local_size + self.size[2])/2/self.dx)
            P3d0 = P3d_raw[start_idx_x:end_idx_x, start_idx_y:end_idx_y, start_idx_z:end_idx_z]
            P3d0 = np.reshape(P3d0, (math.ceil(self.size[0]/self.dx), math.ceil(self.size[1]/self.dx), math.ceil(self.size[2]/self.dx)))
            polarization_profile_3d = np.copy(P3d0.astype(np.complex128), order='C')

            # 1. number of primary cells in the molecule box
            Nx = int(np.ceil(self.size[0] / self.dx))
            Ny = int(np.ceil(self.size[1] / self.dx))
            Nz = int(np.ceil(self.size[2] / self.dx))

            # 2. half-cell offsets for the Yee lattice
            off = {mp.Ex: (0.0, 0.5, 0.5),
                   mp.Ey: (0.5, 0.0, 0.5),
                   mp.Ez: (0.5, 0.5, 0.0)}[self.orientation]

            # coordinate arrays on the **Nx×Ny×Nz** grid
            x = (np.arange(Nx) - (Nx-1)/2 + off[0]) * self.dx
            y = (np.arange(Ny) - (Ny-1)/2 + off[1]) * self.dx
            z = (np.arange(Nz) - (Nz-1)/2 + off[2]) * self.dx
            X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
            R2 = X**2 + Y**2 + Z**2

            if   self.orientation == mp.Ex:
                P = self._polarization_prefactor_3d * X**2 * np.exp(-R2/(2*self.sigma**2))
            elif self.orientation == mp.Ey:
                P = self._polarization_prefactor_3d * Y**2 * np.exp(-R2/(2*self.sigma**2))
            else:  # Ez
                P = self._polarization_prefactor_3d * Z**2 * np.exp(-R2/(2*self.sigma**2))

            polarization_profile_3d = np.ascontiguousarray(P.astype(np.complex128))

        def amp_func_3d_x(R):
            return self._polarization_prefactor_3d * R.x * R.x * np.exp(-(R.x * R.x + R.y * R.y + R.z * R.z) / 2.0 / self.sigma**2)
        def amp_func_3d_y(R):
            return self._polarization_prefactor_3d * R.y * R.y * np.exp(-(R.x * R.x + R.y * R.y + R.z * R.z) / 2.0 / self.sigma**2)
        def amp_func_3d_z(R):
            return self._polarization_prefactor_3d * R.z * R.z * np.exp(-(R.x * R.x + R.y * R.y + R.z * R.z) / 2.0 / self.sigma**2)
        def amp_func_2d(R):
            return self._polarization_prefactor_2d * np.exp(-(R.x * R.x + R.y * R.y) / 2.0 / self.sigma**2)

        # note that using amp_data is x1.5 faster than using amp_func but for now we use amp_func for simplicity
        self.sources = [mp.Source(mp.CustomSource(src_func=lambda t: 1.0),
                            component=self.orientation,
                            center=self.center,
                            size=self.size,
                            amplitude=0.0,
                            amp_data=polarization_profile_2d if self.dimensions == 2 else polarization_profile_3d
                            # amp_func=amp_func_2d if self.dimensions == 2 else (amp_func_3d_x if self.orientation == mp.Ex else (amp_func_3d_y if self.orientation == mp.Ey else amp_func_3d_z))
                            ) ]

    def _propagate(self, int_ep):
        """
        TLS: Propagate the quantum molecular dynamics given the light-matter coupling int_ep.
        """
        # update the density matrix for one time step
        # the 2*pi factor in front of int_ep is to convert the energy to MEEP units, 
        # note that the omega in expHs has already included the 2*pi factor
        U = np.dot(self.expHs, np.dot(expm( (1j * self.dt * int_ep) * self.SIGMAX), self.expHs))
        self.rho = np.dot(np.dot(U, self.rho), U.conj().transpose())
        self.rho_history.append(self.rho.copy())

    def _update_source_amplitude(self):
        """
        TLS: Update the source amplitude after propagating this molecule for one time step.
        """
        amp = -2.0 * self.omega * np.imag(self.rho[0, 1])
        for s in self.sources:
            s.amplitude = amp

    def _calculate_ep_integral(self, sim):
        """
        TLS: Calculate the integral of the electric field over the molecule's volume.
        """
        int_ep = 0.
        if self.dimensions == 2:
            int_ep = sim.integrate_field_function([mp.Ez],
                            lambda R, ez : self._polarization_prefactor_2d * exp(-((R.x - self.center.x) * (R.x - self.center.x) + (R.y - self.center.y) * (R.y - self.center.y)) / (2.0 * self.sigma**2)) * (ez),
                            mp.Volume(size=self.size, center=self.center))
        elif self.dimensions == 3:
            if self.orientation == mp.Ez:
                int_ep = sim.integrate_field_function([mp.Ez],
                                lambda R, ez : self._polarization_prefactor_3d * (R.z-self.center.z)*(R.z-self.center.z) * exp(-((R.x-self.center.x)*(R.x-self.center.x) + (R.y-self.center.y)*(R.y-self.center.y) + (R.z-self.center.z)*(R.z-self.center.z)) / (2.0 * self.sigma**2)) * (ez),
                                mp.Volume(size=self.size, center=self.center))
            elif self.orientation == mp.Ex:
                int_ep = sim.integrate_field_function([mp.Ex],
                                lambda R, ex : self._polarization_prefactor_3d * (R.x-self.center.x)*(R.x-self.center.x) * exp(-((R.x-self.center.x)*(R.x-self.center.x) + (R.y-self.center.y)*(R.y-self.center.y) + (R.z-self.center.z)*(R.z-self.center.z)) / (2.0 * self.sigma**2)) * (ex),
                                mp.Volume(size=self.size, center=self.center))
            elif self.orientation == mp.Ey:
                int_ep = sim.integrate_field_function([mp.Ey],
                                lambda R, ey : self._polarization_prefactor_3d * (R.y-self.center.y)*(R.y-self.center.y) * exp(-((R.x-self.center.x)*(R.x-self.center.x) + (R.y-self.center.y)*(R.y-self.center.y) + (R.z-self.center.z)*(R.z-self.center.z)) / (2.0 * self.sigma**2)) * (ey),
                                mp.Volume(size=self.size, center=self.center))
        return int_ep
    
    def print_dynamics(self):
        for idx, rho in enumerate(self.rho_history):
            print(f"Time step {idx}: Excited-state population = {np.real(rho[1,1])}")


class SocketMolecule(DummyMolecule):
    """
    TBD
    """

    def __init__(self, resolution, center=mp.Vector3(), size=mp.Vector3(1, 1, 1), dimensions=2, sigma=1.0):
        """
        Initialize the Socket molecule
        """
        super().__init__(dt=0.5 / resolution, dx=1.0 / resolution, center=center, size=size)
        self.dimensions = dimensions
        self.sigma = sigma

        # the polarization prefactor in 3D excluding the transient dipole moment
        self._polarization_prefactor_3d = 1.0 / (2.0 * np.pi)**(1.5) / self.sigma**5 #* self.dipole_moment
        # the polarization prefactor in 2D excluding the transient dipole moment
        self._polarization_prefactor_2d = 1.0 / (2.0 * np.pi)**(1.0) / self.sigma**2 #* self.dipole_moment

        # initialize the sources for the RT-TDDFT molecule
        self._init_sources()

    def _init_sources(self):
        """ 
        Initialize the sources for the molecule.
        This method sets up the polarization profile based on the molecule's properties.

        For general purposes, each molecule contains three polarization directions (Ex, Ey, Ez) in 3D or one polarization direction (Ez) in 2D.
        """
        def amp_func_3d_x(R):
            return self._polarization_prefactor_3d * R.x**2 * np.exp(-(R.x**2 + R.y**2 + R.z**2) / 2.0 / self.sigma**2)
        def amp_func_3d_y(R):
            return self._polarization_prefactor_3d * R.y**2 * np.exp(-(R.x**2 + R.y**2 + R.z**2) / 2.0 / self.sigma**2)
        def amp_func_3d_z(R):
            return self._polarization_prefactor_3d * R.z**2 * np.exp(-(R.x**2 + R.y**2 + R.z**2) / 2.0 / self.sigma**2)
        def amp_func_2d(R):
            return self._polarization_prefactor_2d * np.exp(-(R.x**2 + R.y**2) / 2.0 / self.sigma**2)
        if self.dimensions == 2:
            self.sources = [mp.Source(mp.CustomSource(src_func=lambda t: 1.0),
                                                  component=mp.Ez,
                                                  center=self.center,
                                                  size=self.size,
                                                  amplitude=0.0,
                                                  amp_func=amp_func_2d)]
        elif self.dimensions == 3:
            self.sources = [mp.Source(mp.CustomSource(src_func=lambda t: 1.0),
                                                  component=mp.Ex,
                                                  center=self.center,
                                                  size=self.size,
                                                  amplitude=0.0,
                                                  amp_func=amp_func_3d_x),
                        mp.Source(mp.CustomSource(src_func=lambda t: 1.0),
                                                  component=mp.Ey,
                                                  center=self.center,
                                                  size=self.size,
                                                  amplitude=0.0,
                                                  amp_func=amp_func_3d_y),
                        mp.Source(mp.CustomSource(src_func=lambda t: 1.0),
                                                  component=mp.Ez,
                                                  center=self.center,
                                                  size=self.size,
                                                  amplitude=0.0,
                                                  amp_func=amp_func_3d_z)]
        else:
            raise ValueError("Socket Molecule only supports 2D and 3D simulations. Please use dimensions=2 or dimensions=3.")