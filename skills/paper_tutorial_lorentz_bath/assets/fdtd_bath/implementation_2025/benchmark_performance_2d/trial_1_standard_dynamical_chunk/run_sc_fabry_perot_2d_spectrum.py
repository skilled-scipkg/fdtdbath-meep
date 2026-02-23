'''
Script to submit 2D strong coupling simulation with the Lorentz-Bath model for evaluating the linear dispersion relation
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
parser.add_argument('-e', '--endtime', type=float, default=10.0, help="Value for controlling the time of simulation")
parser.add_argument('-c', '--incidentangle', type=float, default=0.0, help="Value for controlling the incident angle of the incoming Gaussian pulse")
parser.add_argument('-y', '--ylength', type=float, default=20.0, help="Value for the length of the y direction simulaiton")

# Parse the command-line arguments.
args = parser.parse_args()

resolution = 400 
frequency_probe = args.fprobe
frequency_lorentz = 1.0
frequency_width = args.wprobe
pulse_intensity = args.intensity
number_frequencies = 1200
pml_thickness = 2.0
power_decay_target = 1e-9

gamma_lorentz = args.gamma
sigma_lorentz = args.sigma

t_wavelength = 1.0
end_time = args.endtime

# reset the cavity mirror refractive index to be vacuum or materials
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

# set the rotation angle of the incoming Gaussian pulse
incident_angle = args.incidentangle
incident_angle_rad = np.radians(incident_angle)
k_point = mp.Vector3(1.0, 0, 0).rotate(mp.Vector3(0, 0, 1), incident_angle_rad)

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
print("In 2D simulations, the incoming Gaussian pulse is rotated by %f degrees with respective to the X direction" %incident_angle)


layer_indexes = np.array([1.0, n1, 1, n1, 1.0])
layer_thicknesses = np.array([0.48, t_wavelength*0.02, t_wavelength, t_wavelength*0.02, 0.48])
layer_thicknesses[0] += pml_thickness
layer_thicknesses[-1] += pml_thickness
total_length = np.sum(layer_thicknesses)
layer_centers = np.cumsum(layer_thicknesses) - layer_thicknesses / 2
layer_centers = layer_centers - total_length / 2

source_location = mp.Vector3(layer_centers[0] - layer_thicknesses[0] / 2 + pml_thickness*1.1, 0., 0)
transmission_monitor_location = mp.Vector3(layer_centers[-1] - pml_thickness / 2, 0., 0)
reflection_monitor_location = mp.Vector3(layer_centers[0] - layer_thicknesses[0] / 2 + pml_thickness*1.1, 0., 0)

cell_y_size = args.ylength # run 2D simulations
cell_size = mp.Vector3(total_length, cell_y_size, 0)

# middle layer cell_size
cell_size_confined = mp.Vector3(layer_thicknesses[2], cell_y_size, 0)
layer_center_confined = mp.Vector3(layer_centers[2], 0, 0)
# time spacing for storing the EM energy
dt_store = 0.1

pml_layers = [mp.PML(thickness=pml_thickness, direction=mp.X)]

geometry = [
        mp.Block(
            mp.Vector3(layer_thicknesses[i], mp.inf, mp.inf),
            center=mp.Vector3(layer_centers[i], 0, 0),
            material=mp.Medium(index=layer_indexes[i])
        )
        for i in range(layer_thicknesses.size)
]

'''
sources = [
    mp.Source(
        mp.GaussianSource(frequency=frequency_probe, fwidth=frequency_width),
        component=mp.Ez,
        center=source_location,
        amplitude=pulse_intensity**0.5,
        size=mp.Vector3(0, cell_y_size, cell_y_size)
    )
]
'''

if incident_angle_rad == 0:
    direction = mp.AUTOMATIC
    eig_parity = mp.EVEN_Y + mp.ODD_Z
    symmetries = [mp.Mirror(mp.Y)]
    eig_vol = None
else:
    direction = mp.NO_DIRECTION
    eig_parity = mp.ODD_Z
    symmetries = []
    eig_vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(0, 1 / resolution, 0))

sources = [
    mp.EigenModeSource(
        mp.GaussianSource(frequency=frequency_probe, fwidth=frequency_width),
        center=source_location,
        size=mp.Vector3(0, cell_y_size, 0),
        direction=direction,
        eig_kpoint=k_point,
        amplitude=pulse_intensity**0.5,
        eig_band=1,
        eig_parity=eig_parity,
        eig_vol=eig_vol,
    )
]

# 1. Simulation in vacuum

size_monitor = 0.1

# 2. Simulation with the cavity and material geometry

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

from meep.chunk_balancer import ChunkBalancer
from meep.timing_measurements import MeepTimingMeasurements

chunk_balancer = ChunkBalancer()

sim = mp.Simulation(
    cell_size=cell_size,
    sources=sources,
    resolution=resolution,
    geometry=geometry,
    boundary_layers=pml_layers,
    k_point=k_point,
)

# Start of the chunk balancing procedure. Sometimes this chunk balancing procedure is not needed, 
# and the whole chunk balancing block can be commented out; see the end comment below.

sim.run(until=0.04)


for nrun in range(3):
    # Compute and save chunk layout for next run
    timings = MeepTimingMeasurements.new_from_simulation(sim)
    chunk_layout = sim.chunk_layout
    chunk_volumes = sim.structure.get_chunk_volumes()
    chunk_owners = sim.structure.get_chunk_owners()
    next_chunk_layout = ChunkBalancer().compute_new_chunk_layout(
        timings,
        chunk_layout,
        chunk_volumes,
        chunk_owners,
        sensitivity=0.4)
    sim = mp.Simulation(
        cell_size=cell_size,
        sources=sources,
        resolution=resolution,
        geometry=geometry,
        boundary_layers=pml_layers,
        k_point=k_point,
        chunk_layout=next_chunk_layout
    )
    sim.run(until=0.04)

sim = mp.Simulation(
    cell_size=cell_size,
    sources=sources,
    resolution=resolution,
    geometry=geometry,
    boundary_layers=pml_layers,
    k_point=k_point,
    chunk_layout=next_chunk_layout
)
# End of the chunk balancing procedure, the above chunk balancing block can be commented out if not needed.

sim.run(
    until=end_time
)


