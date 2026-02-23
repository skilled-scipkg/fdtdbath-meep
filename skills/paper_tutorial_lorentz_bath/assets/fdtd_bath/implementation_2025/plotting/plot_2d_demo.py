import numpy as np
import columnplots as clp
from scipy.optimize import curve_fit

path_1d = "../2d_harmonic_broadlinewidth/"

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

def prepare_data_scanning_gamma(sigma_lst=["0.001", "0.002", "0.003"]):
    xs_lorentz, ys_lorentz = [], []
    xs_lb_lorentzian, ys_lb_lorentzian = [], [] 

    for sigma in sigma_lst:
        # load the data for the standard Lorentzian medium
        path_incav_lorentz = path_1d + "/incav_lorentz_sigma_%s_angle_0/spectrum.txt" %sigma
        data_incav_lorentz = load_spectrum(path_incav_lorentz)

        # load the data for the Lorentz+Bath model with lorentzian distribution
        path_incav_lb_lorentzian = path_1d + "/incav_lb_lorentzian_sigma_%s_angle_0/spectrum.txt" %sigma
        data_incav_lb_lorentzian = load_spectrum(path_incav_lb_lorentzian)

        xs_lorentz.append(data_incav_lorentz['frequency'])
        ys_lorentz.append(data_incav_lorentz['transmittance'])
        xs_lb_lorentzian.append(data_incav_lb_lorentzian['frequency'])
        ys_lb_lorentzian.append(data_incav_lb_lorentzian['transmittance'])

    return xs_lorentz, ys_lorentz, xs_lb_lorentzian, ys_lb_lorentzian

labels = ["Lorentz model", "uniform bath", "Lorentzian bath"]
labels = ["Lorentz", "Lorentz-Bath(U)", "Lorentz-Bath(L)"]
colors = [clp.darkgray, clp.cyan, clp.red_economist]
linestyles = ["-", "--", "-"]
# Set up the figure and axes
ax = clp.initialize(1, 1, width=4.3, height=4.3*0.618, fontsize=12, fontname='Arial', LaTeX=False)

# Part 1: 2D spectrum inside the cavity versus the Lorentzian medium sigma (absorption coefficient)
sigma_lst = ["0.002", "0.02"]
xs_lorentz, ys_lorentz, xs_lb_lorentzian, ys_lb_lorentzian = prepare_data_scanning_gamma(sigma_lst=sigma_lst)
# lift the Lorentzian data when the index is increased
ys_lorentz = [y + idx*1. for idx, y in enumerate(ys_lorentz)]
ys_lb_lorentzian = [y + idx*1. for idx, y in enumerate(ys_lb_lorentzian)]

colors = [clp.red_economist] * len(xs_lb_lorentzian) + [clp.darkgray] * len(xs_lorentz)

for i in range(len(xs_lb_lorentzian)):
    ax.fill_between(xs_lorentz[i], ys_lorentz[i]*0.0 + i, ys_lorentz[i], color=clp.darkgray)

# Set up the figure and axes
clp.plotone(xs_lb_lorentzian, ys_lb_lorentzian, ax,
            xlabel=r"$\nu = \omega/2\pi$ [µm$^{-1}$]", ylabel="transmittance",
            xlim=[0.45, 1.55], ylim=[0, 2.], lw=1.3,
            showlegend=False, linestyles=["-"] * len(xs_lb_lorentzian) + ["--"] * len(xs_lorentz),
            colors=colors)

ax.text(0.6, 0.0+0.1, r"$\sigma_{\text{L}}=2\times 10^{-3}$", fontsize=11, color="0.5")
ax.text(0.6, 1.+0.1, r"$\sigma_{\text{L}}=2\times 10^{-2}$", fontsize=11, color="0.5")
# axes[1].text(0.2, 0.85, "cavity on", fontsize=11, transform=axes[1].transAxes, ha="center", va="center", alpha=0.5)
# 
# add inset zoom in parts to better show the polaritons
for i in range(2):  
    axins = ax.inset_axes(
        [0.53, 0.14 + 0.5 * i, 0.3, 0.3], xlim=(0.9, 1.1), ylim=(0+i, 0.18+i), xticklabels=[], yticklabels=[])
    # axins.fill_between(xs_lorentz[i], ys_lorentz[i]*0.0 + i, ys_lorentz[i], color="0.5", alpha=0.3)
    clp.plotone([xs_lb_lorentzian[i]], [ys_lb_lorentzian[i]], axins, showlegend=False, 
            colors=[clp.red, "k"], linestyles=["-", "--"], lw=1.2)
    ax.indicate_inset_zoom(axins, edgecolor="0.5", alpha=0.5)
    axins.set_facecolor(clp.lightgray_background)
    axins.fill_between(xs_lorentz[i], ys_lorentz[i]*0.0 + i, ys_lorentz[i], color=clp.darkgray)

# reset y decimal formatter
import matplotlib.ticker as mticker
formatter = mticker.FormatStrFormatter('%.1f')
ax.yaxis.set_major_formatter(formatter)

ax.set_facecolor(clp.lightgray_background)

clp.adjust(tight_layout=True, savefile="2d_spectrum_demo.pdf")
