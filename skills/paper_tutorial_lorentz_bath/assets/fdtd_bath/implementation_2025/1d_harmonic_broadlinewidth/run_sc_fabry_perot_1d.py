'''
Script to submit 1D strong coupling simulation with the Lorentz-Bath model
'''

import numpy as np
import meep as mp
import argparse 

# create the parser to retrive the parameters from outside input 
parser = argparse.ArgumentParser(description="Need job parameters to run strong coupling simulations in a 1D Fabry-Perot cavity")
parser.add_argument('-f', '--fprobe', type=float, default=1.0, help="Value for the incident Gaussian pulse frequency")
parser.add_argument('-w', '--wprobe', type=float, default=1.2, help="Value for the incident Gaussian pulse width")
parser.add_argument('-g', '--gamma', type=float,  default=0.04, help="Value for the Lorentz medium gamma, damping rate")
parser.add_argument('-s', '--sigma', type=float,  default=0.005, help="Value for the Lorentz medium sigma, absorption strength")
parser.add_argument('-r', '--refractiveindex', type=float,  default=10.0, help="Value for the refractive index of the mirrors")
parser.add_argument('-n', '--numbath', type=int,  default=0, help="Value for the number of the bath oscillators")
parser.add_argument('-b', '--bathwidth', type=float,  default=10.0, help="Value for the width of frequency range of the bath oscillators, in the units of gamma")
parser.add_argument('-d', '--dephasingratio', type=float,  default=1.0, help="Value for the ratio of dephasing mechanism versus the total decay of the bright mode")
parser.add_argument('-l', '--bathform', type=str,  default="lorentzian", help="Value for form of the bath oscillators, ")
parser.add_argument('-k', '--decayratio', type=float,  default=0.0, help="Value for the ratio between the bath decay rate versus the Lorentz oscillator decay rate")
parser.add_argument('-a', '--anharmonicity', type=float, default=0.0, help="Value for the anharmonicity parameter of the bath oscillators (not the bright mode)")
parser.add_argument('-i', '--intensity', type=float, default=1.0, help="Value for the intensity of the incident Gaussian pulse, amplitude = intensity^0.5")
parser.add_argument('-t', '--transform', type=bool, default=False, help="Value for determining whether to transform the Lorentz-Bath model to N+1 independent Lorentz oscillators")
parser.add_argument('-p', '--percentageoscillators', type=float, default=0.0, help="Value for determining removing the percentage of the side oscillators when transforming the Lorentz-Bath model to N+1 independent Lorentz oscillators")
parser.add_argument('-m', '--renormalizationcoeff', type=float, default=0.0, help="Value for reducing the relative strengths of the side oscillators when transforming the Lorentz-Bath model to N+1 independent Lorentz oscillators")
parser.add_argument('-e', '--endtime', type=float, default=200.0, help="Value for controlling the time of simulation")

# Parse the command-line arguments.
args = parser.parse_args()

resolution = 100 
frequency_probe = args.fprobe
frequency_lorentz = 1.0
frequency_width = args.wprobe
pulse_intensity = args.intensity
number_frequencies = 1200
pml_thickness = 0.5
power_decay_target = 1e-9

gamma_lorentz = args.gamma
sigma_lorentz = args.sigma

t_wavelength = 1.0
end_time = args.endtime

# reset the cavity as free space
n1 = args.refractiveindex

# set other parameters for the medium
num_bath = args.numbath
bath_width = args.bathwidth
dephasing_ratio = args.dephasingratio
bath_form = args.bathform
decay_ratio = args.decayratio
bath_anharmonicity = args.anharmonicity

# set parameters for the basis transformation
transform = args.transform
percentage_oscillators = args.percentageoscillators
renormalization_coeff = args.renormalizationcoeff

print("Receiving the input parameters from the terminal")
print("frequency_probe:", frequency_probe)
print("frequency_width:", frequency_width)
print("pulse intensity:", pulse_intensity)
print("gamma_lorentz:", gamma_lorentz)
print("sigma_lorentz:", sigma_lorentz)
print("refractive index of mirrors:", n1)
print("num_bath:", num_bath)
print("bath_width (as a function of gamma_lorentz):", bath_width)
print("bath_form:", bath_form)
print("bath decay_ratio (as a function of gamma_lorentz):", decay_ratio)
print("bath anharmonicity:", bath_anharmonicity)
print("dephasing from bright state to bath (as a function of gamma_lorentz):", dephasing_ratio)
if transform:
    print("transform lorentz-bath model to independent lorentz oscillators:", transform)
    print("percentage of N+1 independent oscillators to be removed (if a transform is taken):", percentage_oscillators)
    print("renormalization coefficient for the side oscillators (if a transform is taken):", renormalization_coeff)


layer_indexes = np.array([1.0, n1, 1, n1, 1.0])
layer_thicknesses = np.array([1.98, t_wavelength*0.02, t_wavelength, t_wavelength*0.02, 1.98])
layer_thicknesses[0] += pml_thickness
layer_thicknesses[-1] += pml_thickness
total_length = np.sum(layer_thicknesses)
layer_centers = np.cumsum(layer_thicknesses) - layer_thicknesses / 2
layer_centers = layer_centers - total_length / 2

source_location = mp.Vector3(layer_centers[0], 0, 0)
transmission_monitor_location = mp.Vector3(layer_centers[-1] - pml_thickness / 2, 0, 0)
reflection_monitor_location = mp.Vector3(layer_centers[0] + layer_thicknesses[0] / 4, 0, 0)

cell_y_size = 0 # effective 1D simulations
cell_size = mp.Vector3(total_length, cell_y_size, 0)

# middle layer cell_size
cell_size_confined = mp.Vector3(layer_thicknesses[2], cell_y_size, 0)
layer_center_confined = mp.Vector3(layer_centers[2], 0, 0)
# time spacing for storing the EM energy
dt_store = 0.1

pml_layers = [mp.PML(thickness=pml_thickness)]

geometry = [
        mp.Block(
            mp.Vector3(layer_thicknesses[i], mp.inf, mp.inf),
            center=mp.Vector3(layer_centers[i], 0, 0),
            material=mp.Medium(index=layer_indexes[i])
        )
        for i in range(layer_thicknesses.size)
]

sources = [
    mp.Source(
        mp.GaussianSource(frequency=frequency_probe, fwidth=frequency_width),
        component=mp.Ez,
        center=source_location,
        amplitude=pulse_intensity**0.5,
        size=mp.Vector3(0, cell_y_size, cell_y_size)
    )
]

sim = mp.Simulation(
    cell_size=cell_size,
    sources=sources,
    resolution=resolution,
    boundary_layers=pml_layers
)

incident_region = mp.FluxRegion(
    center=reflection_monitor_location, size=mp.Vector3(1.0, 1.0, 0),
    weight=1.0, direction=mp.X
)
incident_flux_monitor = sim.add_flux(
    frequency_probe, frequency_width, number_frequencies, incident_region
)

sim.run(until_after_sources=end_time)
# sim.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Ez,
#     transmission_monitor_location, power_decay_target))

flux_freq = mp.get_flux_freqs(incident_flux_monitor)
frequencies = np.array(flux_freq)
empty_flux_int = mp.get_fluxes(incident_flux_monitor)
empty_flux_data = sim.get_flux_data(incident_flux_monitor)
empty_intensities = np.array(empty_flux_int)

incident_flux_to_subtract = sim.get_flux_data(incident_flux_monitor)



geometry = [
    mp.Block(
        mp.Vector3(layer_thicknesses[i], mp.inf, mp.inf),
        center=mp.Vector3(layer_centers[i], 0, 0),
        material=mp.Medium(index=layer_indexes[i])
    )
    for i in range(layer_thicknesses.size)
]

if num_bath == 0:
    susceptibilities = [
    mp.LorentzianSusceptibility(
        frequency=frequency_lorentz,
        gamma=gamma_lorentz,
        sigma=sigma_lorentz,
        )
    ]
else:
    lb_model = mp.BathLorentzianSusceptibility(
            frequency=frequency_lorentz,
            gamma=gamma_lorentz*(1.0 - dephasing_ratio),
            sigma=sigma_lorentz,
            num_bath=num_bath,
            bath_form=bath_form,
            bath_width=bath_width * gamma_lorentz,
            bath_dephasing=dephasing_ratio * gamma_lorentz,
            bath_gamma=gamma_lorentz*decay_ratio,
            bath_anharmonicities=[bath_anharmonicity]*num_bath,
        )
    
    susceptibilities = [
        lb_model
    ]

    if transform:
        # convert the Lorentz-Bath model to N+1 independent Lorentz oscillators
        lorentzians = lb_model.obtain_independent_lorentzians()
        if percentage_oscillators > 0.0:
            # remove the percentage of the side oscillators
            num_lorentz = int(len(lorentzians) * (percentage_oscillators / 2))
            lorentzians = lorentzians[num_lorentz:-num_lorentz]
            print("### Removing %d side oscillators at both sides" %num_lorentz)
        if renormalization_coeff > 0.0:
            # define the second approach for removing the side oscillators
            delta  = 1.0 / renormalization_coeff * gamma_lorentz
            for l in lorentzians:
                freq = l["frequency"]
                renormalization_factor = delta**2 / ((freq - frequency_lorentz) ** 2 + delta**2)
                l["sigma_weight"] *= renormalization_factor
                print("### Renormalizing the Lorentzian oscillator at %f with a factor %.3E" %(freq, renormalization_factor))
        susceptibilities = [
                mp.LorentzianSusceptibility(
                frequency=l["frequency"],
                gamma=l["gamma"],
                sigma=l["sigma_weight"] * sigma_lorentz,
                )
                for l in lorentzians
            ]
    
lorentzian_material = mp.Medium(
    epsilon=1.0, E_susceptibilities=susceptibilities
)

idx_middle = 2
geometry[idx_middle] = mp.Block(
    mp.Vector3(layer_thicknesses[idx_middle], mp.inf, mp.inf),
    center=mp.Vector3(layer_centers[idx_middle], 0, 0),
    material=lorentzian_material
)

sim = mp.Simulation(
    cell_size=cell_size,
    sources=sources,
    resolution=resolution,
    geometry=geometry,
    boundary_layers=pml_layers
)

transmission_region = mp.FluxRegion(
    center=transmission_monitor_location, size=mp.Vector3(1.0, 1.0, 0),
    weight=1.0, direction=mp.X
)
transmission_flux_monitor = sim.add_flux(
    frequency_probe, frequency_width, number_frequencies, transmission_region
)
reflection_region = incident_region
reflection_flux_monitor = sim.add_flux(
    frequency_probe, frequency_width, number_frequencies, reflection_region
)
sim.load_minus_flux_data(reflection_flux_monitor, incident_flux_to_subtract)

# record the EM energy stored within the cavity
field_energy_snapshots = []

def capture_field_energy(sim):
    integration_box = mp.Volume(size=cell_size_confined, center=layer_center_confined)
    (Ex,Ey,Ez) = [sim.get_array(vol=integration_box, component=c, cmplx=True) for c in [mp.Ex, mp.Ey, mp.Ez]]
    eps = sim.get_array(vol=integration_box, component=mp.Dielectric)
    (Hx,Hy,Hz) = [sim.get_array(vol=integration_box, component=c, cmplx=True) for c in [mp.Hx, mp.Hy, mp.Hz]]
    perm = sim.get_array(vol=integration_box, component=mp.Permeability)
    # integrate the EM energy: https://meep.readthedocs.io/en/latest/Python_User_Interface/#output-functions
    (x, y, z, w) = sim.get_array_metadata(vol=integration_box)
    electric_energy_density = 0.5 * np.real(eps*(np.conj(Ex)*Ex + np.conj(Ey)*Ey + np.conj(Ez)*Ez))   
    magnetic_energy_density = 0.5 * np.real(1/perm*(np.conj(Hx)*Hx + np.conj(Hy)*Hy + np.conj(Hz)*Hz)) 
    energy_density = (electric_energy_density + magnetic_energy_density)          
    energy = np.sum(w*energy_density)  
    field_energy_snapshots.append(energy)                                 

sim.run(
    # mp.after_sources(mp.Harminv(mp.Ez, mp.Vector3(layer_centers[10], 0, 0), frequency, frequency_width)),
    mp.at_every(dt_store, capture_field_energy), 
    until_after_sources=end_time
)

# calculate the spectrum
incident_flux = np.array(mp.get_fluxes(incident_flux_monitor))
transmitted_flux = np.array(mp.get_fluxes(transmission_flux_monitor))
reflected_flux = np.array(mp.get_fluxes(reflection_flux_monitor))
freqs_output = mp.get_flux_freqs(reflection_flux_monitor)
reflectance = -reflected_flux / incident_flux
transmittance = transmitted_flux / incident_flux
loss = 1.0 - reflectance - transmittance

# calculate the polariton decay EM energy dynamics
t = np.linspace(0, (len(field_energy_snapshots)-1)*dt_store, len(field_energy_snapshots))
field_energy = np.array(field_energy_snapshots)

# save the data
data_freq = np.zeros((len(freqs_output), 4))
data_freq[:, 0] = freqs_output
data_freq[:, 1] = loss
data_freq[:, 2] = reflectance
data_freq[:, 3] = transmittance
np.savetxt("spectrum.txt", data_freq, header="freq, loss, reflectance, transmittance (in meep units)", fmt="%.6E")

data_energy = np.zeros((len(t), 2))
data_energy[:, 0] = t
data_energy[:, 1] = field_energy
np.savetxt("energy_dynamics.txt", data_energy, header="t, EM energy in cavity (in meep units)", fmt="%.6E")


'''
# visualization 
import columnplots as clp
axes = clp.initialize(4, 1, width=4, height=10)
clp.plotone(
    [freqs_output], [loss], axes[0],
    xlabel="frequency", ylabel="loss intensity",
    colors=["k", "r", "y", "g", "b"], alphaspacing=0.1,
    labels=["$N_{bath}$= %d" %num_bath],
)
clp.plotone(
    [freqs_output], [reflectance], axes[1],
    xlabel="frequency", ylabel="R intensity",
    colors=["k", "r", "y", "g", "b"], alphaspacing=0.1,
    labels=["$N_{bath}$= %d" %num_bath],
)
clp.plotone(
    [freqs_output], [transmittance], axes[2],
    xlabel="frequency", ylabel="T intensity",
    colors=["k", "r", "y", "g", "b"], alphaspacing=0.1,
    labels=["$N_{bath}$= %d" %num_bath],
)

clp.plotone([t], [field_energy], axes[3],
    xlabel="time", ylabel="EM energy",
    colors=["k", "r", "y", "g", "b"], alphaspacing=0.1,
    labels=["$N_{bath}$= %d" %num_bath],
)
clp.adjust()
'''
