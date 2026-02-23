import numpy as np
import columnplots as clp
from scipy.optimize import curve_fit

# path_1d = "../1d_harmonic/"
path_1d = "../1d_harmonic_broadlinewidth/"

# first, load the data from files
def load_spectrum(file_path):
    """
    Load the spectrum data from a file.
    """
    # Assuming the file is in a specific format, e.g., two columns: wavelength and flux
    data = np.loadtxt(file_path)
    freq, loss, reflectance, transmittance = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    # Prepare a dictionary to hold the data
    data_dict = {
         'frequency': freq,
         'loss': loss,
         'reflectance': reflectance,
         'transmittance': transmittance
        }
    return data_dict

def load_energy_dynamics(file_path):
    """
    Load the spectrum data from a file.
    """
    # Assuming the file is in a specific format, e.g., two columns: wavelength and flux
    data = np.loadtxt(file_path)
    t, em_energy = data[:, 0], data[:, 1]
    # Prepare a dictionary to hold the data
    data_dict = {
         't': t,
         'em_energy': em_energy #/ np.max(em_energy)
        }
    return data_dict

def prepare_plotting_data():
    # load the data for the standard Lorentzian medium
    path_outcav_lorentz = path_1d + "/outcav_lorentz/spectrum.txt"
    data_outcav_lorentz = load_spectrum(path_outcav_lorentz)

    # load the data for the Lorentz+Bath model with uniform distribution
    path_outcav_lb_uniform = path_1d + "/outcav_lb_uniform/spectrum.txt"
    data_outcav_lb_uniform = load_spectrum(path_outcav_lb_uniform)

    # load the data for the Lorentz+Bath model with lorentzian distribution
    path_outcav_lb_lorentzian = path_1d + "/outcav_lb_lorentzian/spectrum.txt"
    data_outcav_lb_lorentzian = load_spectrum(path_outcav_lb_lorentzian)

    xs_outcav = [data_outcav_lorentz['frequency'], data_outcav_lb_uniform['frequency'], data_outcav_lb_lorentzian['frequency']]
    ys_outcav = [data_outcav_lorentz['transmittance'], data_outcav_lb_uniform['transmittance'], data_outcav_lb_lorentzian['transmittance']]

    return xs_outcav, ys_outcav


def prepare_data_scanning_gamma(sigma_lst=["0.001", "0.002", "0.003"]):
    xs_lorentz, ys_lorentz = [], []
    xs_lb_lorentzian, ys_lb_lorentzian = [], [] 

    for sigma in sigma_lst:
        # load the data for the standard Lorentzian medium
        path_incav_lorentz = path_1d + "/incav_lorentz_sigma_%s/spectrum.txt" %sigma
        data_incav_lorentz = load_spectrum(path_incav_lorentz)

        # load the data for the Lorentz+Bath model with lorentzian distribution
        path_incav_lb_lorentzian = path_1d + "/incav_lb_lorentzian_sigma_%s/spectrum.txt" %sigma
        data_incav_lb_lorentzian = load_spectrum(path_incav_lb_lorentzian)

        xs_lorentz.append(data_incav_lorentz['frequency'])
        ys_lorentz.append(data_incav_lorentz['transmittance'])
        xs_lb_lorentzian.append(data_incav_lb_lorentzian['frequency'])
        ys_lb_lorentzian.append(data_incav_lb_lorentzian['transmittance'])

    return xs_lorentz, ys_lorentz, xs_lb_lorentzian, ys_lb_lorentzian


def prepare_energy_dynamics_scanning_gamma(sigma_lst=["0.001", "0.002", "0.003"], polariton="lp"):
    xs_lorentz, ys_lorentz = [], []
    xs_lb_lorentzian, ys_lb_lorentzian = [], [] 

    for sigma in sigma_lst:
        # load the data for the standard Lorentzian medium
        path_incav_lorentz = path_1d + "/incav_lorentz_sigma_%s/energy_dynamics.txt" %(sigma)
        data_incav_lorentz = load_energy_dynamics(path_incav_lorentz)

        # load the data for the Lorentz+Bath model with lorentzian distribution
        path_incav_lb_lorentzian = path_1d + "/incav_lb_lorentzian_sigma_%s/energy_dynamics.txt" %(sigma)
        data_incav_lb_lorentzian = load_energy_dynamics(path_incav_lb_lorentzian)

        xs_lorentz.append(data_incav_lorentz['t'])
        ys_lorentz.append(data_incav_lorentz['em_energy'])
        xs_lb_lorentzian.append(data_incav_lb_lorentzian['t'])
        ys_lb_lorentzian.append(data_incav_lb_lorentzian['em_energy'])

    return xs_lorentz, ys_lorentz, xs_lb_lorentzian, ys_lb_lorentzian

# Assuming a decay-oscillation pattern
def model_vacuum_rabi_oscillation_decay(t, gamma, rabi):
    return np.exp(-(t-t[0])*gamma) * 0.5 #np.cos(2.0 * np.pi * (t-t[0])*rabi/2)**2

def fit_vacuum_rabi_oscillation_decay(t, y, rabi):
    from functools import partial
    # Initial guess for the parameters
    model_vacuum_rabi_oscillation_decay_fixed = partial(model_vacuum_rabi_oscillation_decay, rabi=rabi)
    params, pcov = curve_fit(model_vacuum_rabi_oscillation_decay_fixed, t, y / y[0])
    return params, pcov

# Part 1: 1D spectrum inside versus outside the cavity
# dataset for 1d spectrum outside the cavity
xs_outcav, ys_outcav = prepare_plotting_data()
labels = ["Lorentz model", "uniform bath", "Lorentzian bath"]
labels = ["Lorentz", "Lorentz-Bath(U)", "Lorentz-Bath(L)"]


colors = [clp.darkgray, clp.cyan, clp.red_economist]
linestyles = ["-", "--", "-"]
# Set up the figure and axes
axes = clp.initialize(3, 1, width=4.3, height=4.3*0.618*3, fontsize=12, fontname='Arial', 
                      labelthem=True, labelsize=12, labelthemPosition=[-0.13, 1.04],
                      LaTeX=False)

clp.plotone(xs_outcav, ys_outcav, axes[0],
            xlim=[0.8, 1.2], labels=labels,
            ylim=[0.4, 1], showlegend=True, 
            legendFaceColor=clp.lightgray_background,
            colors=colors, linestyles=["-", "--", "-"], 
            legendFontSize=9,
            xlabel=r"$\nu = \omega/2\pi$ [µm$^{-1}$]", ylabel="transmittance")

axes[0].fill_between(xs_outcav[0], 1.0, ys_outcav[0], color=clp.darkgray)

# axes[0].text(0.2, 0.85, "cavity off", fontsize=11, transform=axes[0].transAxes, ha="center", va="center", alpha=0.5)

# Part 2: 1D spectrum inside the cavity versus the Lorentzian medium sigma (absorption coefficient)
sigma_lst = ["0.002", "0.020"]
xs_lorentz, ys_lorentz, xs_lb_lorentzian, ys_lb_lorentzian = prepare_data_scanning_gamma(sigma_lst=sigma_lst)
# lift the Lorentzian data when the index is increased
ys_lorentz = [y + idx*1. for idx, y in enumerate(ys_lorentz)]
ys_lb_lorentzian = [y + idx*1. for idx, y in enumerate(ys_lb_lorentzian)]

colors = [clp.red_economist] * len(xs_lb_lorentzian) + [clp.darkgray] * len(xs_lorentz)

for i in range(len(xs_lb_lorentzian)):
    axes[1].fill_between(xs_lorentz[i], ys_lorentz[i]*0.0 + i, ys_lorentz[i], color=clp.darkgray)

# Set up the figure and axes
clp.plotone(xs_lb_lorentzian, ys_lb_lorentzian, axes[1],
            xlabel=r"$\nu = \omega/2\pi$ [µm$^{-1}$]", ylabel="transmittance",
            xlim=[0.45, 1.55], ylim=[0, 2.], lw=1.3,
            showlegend=False, linestyles=["-"] * len(xs_lb_lorentzian) + ["--"] * len(xs_lorentz),
            colors=colors)

axes[1].text(0.6, 0.0+0.1, r"$\sigma=2\times 10^{-3}$", fontsize=11, color="0.5")
axes[1].text(0.6, 1.+0.1, r"$\sigma=2\times 10^{-2}$", fontsize=11, color="0.5")
# axes[1].text(0.2, 0.85, "cavity on", fontsize=11, transform=axes[1].transAxes, ha="center", va="center", alpha=0.5)
# 
# add inset zoom in parts to better show the polaritons
for i in range(2):  
    axins = axes[1].inset_axes(
        [0.53, 0.14 + 0.5 * i, 0.3, 0.3], xlim=(0.9, 1.1), ylim=(0+i, 0.18+i), xticklabels=[], yticklabels=[])
    # axins.fill_between(xs_lorentz[i], ys_lorentz[i]*0.0 + i, ys_lorentz[i], color="0.5", alpha=0.3)
    clp.plotone([xs_lb_lorentzian[i]], [ys_lb_lorentzian[i]], axins, showlegend=False, 
            colors=[clp.red, "k"], linestyles=["-", "--"], lw=1.2)
    axes[1].indicate_inset_zoom(axins, edgecolor="0.5", alpha=0.5)
    axins.set_facecolor(clp.lightgray_background)
    axins.fill_between(xs_lorentz[i], ys_lorentz[i]*0.0 + i, ys_lorentz[i], color=clp.darkgray)

# reset y decimal formatter
import matplotlib.ticker as mticker
formatter = mticker.FormatStrFormatter('%.1f')
axes[1].yaxis.set_major_formatter(formatter)

# Part 3: Nonequilibrium polariton energy dynamics 
# sigma_lst = ["0.001", "0.004", "0.010"]
rescaling = 100.0
xs_lorentz, ys_lorentz, xs_lb_lorentzian, ys_lb_lorentzian = prepare_energy_dynamics_scanning_gamma(sigma_lst=sigma_lst, polariton="lp")
# lift the lorentzian bath data when the index is increased
ys_lb_lorentzian_lifted = [y + idx*0.01 for idx, y in enumerate(ys_lb_lorentzian)]
ys_lorentz_lifted = [y + idx*0.01 for idx, y in enumerate(ys_lorentz)]
ys_lb_lorentzian_lifted = [y * rescaling for y in ys_lb_lorentzian_lifted]
ys_lorentz_lifted = [y * rescaling for y in ys_lorentz_lifted]

clp.plotone(xs_lb_lorentzian + xs_lorentz, ys_lb_lorentzian_lifted + ys_lorentz_lifted, axes[2],
            xlabel="time [µm/$c$]", ylabel="EM energy [arb. units]",
            lw=1.3,
            showlegend=False, linestyles=["-"] * len(xs_lb_lorentzian) + ["--"] * len(xs_lorentz),
            colors=colors,
            xlim=[0, 80], ylim=[0, 0.02*rescaling],
            )

axes[2].text(145/200*80, (0.0+1e-3)*rescaling, r"$\sigma=2\times 10^{-3}$", fontsize=11, color="0.5")
axes[2].text(145/200*80, (0.01+1e-3)*rescaling, r"$\sigma=2\times 10^{-2}$", fontsize=11, color="0.5")

# set background color for all figures
for ax in axes:
    ax.set_facecolor(clp.lightgray_background)

clp.adjust(tight_layout=True, savefile="1d_spectrum_demo.pdf")
