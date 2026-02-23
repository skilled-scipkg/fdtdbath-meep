import numpy as np
import columnplots as clp

def collect_mpi_data(path):
    filename_lb = path + "/info_lb.txt"
    data_lb = np.loadtxt(filename_lb)
    filaname_lorentz = path + "/info_lorentz.txt"
    data_lorentz = np.loadtxt(filaname_lorentz)
    ncpu = data_lb[:, 0]
    timestep_lb = data_lb[:, 1]
    timestep_lorentz = data_lorentz[:, 1]
    return ncpu, timestep_lb, timestep_lorentz

def collect_mpi_data_after_trials(dynamical_chunk=True):
    additional_path = ""
    if dynamical_chunk:
        additional_path = "_dynamical_chunk"
    path_1 = "../benchmark_performance_2d/trial_1_standard%s" %additional_path
    ncpu, timestep_lb, timestep_lorentz = collect_mpi_data(path=path_1)
    for ntrial in range(2, 6):
        path = "../benchmark_performance_2d/trial_%d_standard%s" %(ntrial, additional_path)
        ncpu_local, timestep_lb_local, timestep_lorentz_local = collect_mpi_data(path=path)
        print("timestep_lb", timestep_lb)
        print("timestep_lb_local", timestep_lb_local)
        timestep_lb = np.nanmin([timestep_lb, timestep_lb_local], axis=0)
        timestep_lorentz = np.nanmin([timestep_lorentz, timestep_lorentz_local], axis=0)
        print("nanmin", timestep_lb)

    return ncpu, timestep_lb, timestep_lorentz


# ncpu, timestep_lb, timestep_lorentz = collect_mpi_data_after_trials(dynamical_chunk=False)
ncpu, timestep_lb, timestep_lorentz = collect_mpi_data_after_trials(dynamical_chunk=True)

print(timestep_lb / timestep_lorentz)

# the ideal linear scaling of the computational cost
timestep_lb_ideal = timestep_lb[5] * ncpu[5] / ncpu
timestep_lorentz_ideal = timestep_lorentz[5] * ncpu[5] / ncpu

xs = [ncpu]*2
ys = [timestep_lorentz, timestep_lb]

# Set up the figure and axes
ax = clp.initialize(1, 1, width=4.3, height=4.3*0.618*1.2, fontsize=12, fontname='Arial', LaTeX=False)

colors = [clp.darkgray, clp.red_economist]
labels = ["Lorentz model", "Lorentz-Bath(L)"]

# Set up the figure and axes
clp.plotone(xs, ys, ax, xlabel="CPU number", ylabel="time per step [s]", xlog=True, ylog=True,
            markers=["^", "o"], linestyles=["none", "none"], lw=1.3, xlim=[1, 300],
            colors=colors, legendFaceColor=clp.lightgray_background, labels=labels)

# add ideal scaling
xs_ideal = [ncpu, ncpu]
ys_ideal = [timestep_lorentz_ideal, timestep_lb_ideal]
clp.plotone(xs_ideal, ys_ideal, ax, lw=1.3,
            colors=[clp.darkgray, clp.red_economist], showlegend=False)

ax.xaxis.set_ticks_position('bottom')
# add secondary x axis
secax = ax.secondary_xaxis('top', functions=(lambda x: x/24, lambda x: x/24,))
secax.set_xlabel('node number')
secax.set_xticks([1, 2, 3, 4, 8, 10])

ax.set_facecolor(clp.lightgray_background)

clp.adjust(tight_layout=True, savefile="mpi_benchmark.pdf")
