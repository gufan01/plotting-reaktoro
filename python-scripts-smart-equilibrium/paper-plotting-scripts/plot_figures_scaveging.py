import numpy as np
import matplotlib.pyplot as plt
import os
from natsort import natsorted
import matplotlib as mpl
import json

# Auxiliary time related constants
second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
year = 365 * day

# Discretisation parameters
xl = 0.0          # the x-coordinate of the left boundary
xr = 100.0          # the x-coordinate of the right boundary
ncells = 100      # the number of cells in the discretization
nsteps = 5000
reltol = 1e-1     # relative tolerance
abstol = 1e-8     # absolute tolerance
tol    = 1e-1

D  = 1.0e-9       # the diffusion coefficient (in units of m2/s)
v  = 1.0/day      # the fluid pore velocity (in units of m/s)
dt = 30 * minute  # the time step (in units of s)
T = 60.0          # the temperature (in units of K)
P = 100           # the pressure (in units of Pa)
phi = 0.1         # the porosity

dirichlet = False # the parameter that defines whether Dirichlet BC must be used
smrt_solv = True  # the parameter that defines whether classic or smart
                  # EquilibriumSolver must be used

tag = "-dt-" + "{:d}".format(dt) + \
      "-ncells-" + str(ncells) + \
      "-nsteps-" + str(nsteps) + \
      "-reltol-" + "{:.{}e}".format(reltol, 1) + \
      "-abstol-" + "{:.{}e}".format(abstol, 1)

test_tag_smart = tag + "-smart"
test_tag_class = tag + "-reference"

#folder = 'cpp-reactivetransport-old-demo/results-pitzer-full-with-skipping-1e-13-both-solvers'
#folder = 'results-with-normalization'
#folder = 'results-deltan-with-normalization'
nsteps = 5000    # the number of steps in the reactive transport simulation
plot_at_selected_steps = [1, 10, 60, 120, 240, 480, 960, 1200, 2400, 3600, 4800, 5000]  # the time steps at which the results are plotted
folder_smart = 'results-scaveging-custering-primary-species-dt-4320-ncells-100-nsteps-5000-reltol-5.0e-03-pitzer-smart'
folder_class = 'results-scaveging-custering-primary-species-dt-4320-ncells-100-nsteps-5000-pitzer-reference'
folder_general = "comparison-plots-scaveging-custering-primary-species-dt-4320-ncells-100-nsteps-5000-reltol-5.0e-03-pitzer-smart"
#folder_smart = 'cpp-reactivetransport-old-demo/results-pitzer-dt-3594.24-ncells-100-nsteps-10000-eqreltol-1.0e-01-eqabstol-1.0e-08-smart'
#folder_class = 'cpp-reactivetransport-old-demo/results-pitzer-dt-3594.24-ncells-100-nsteps-10000-reference'
#folder_general = "plots-results-pitzer-dt-3594.24-ncells-100-nsteps-10000"
os.system('mkdir -p ' + folder_general)

fillz = len(str(123))

# Indices of the loaded data to plot
indx_pH = 0
indx_Hcation = 1
indx_HSanion = 2
indx_S2anion = 3
indx_SO4anion = 4
indx_HSO4anion = 5
indx_H2Saq = 6
indx_pyrrhotite = 7
indx_siderite = 8
indx_phase_vol_pyrrhotite = 9
indx_phase_vol_siderite = 10

# Plotting params
circ_area = 6 ** 2
# zoom = 0.5
zoom = 0.5
custom_font = { }
time_steps = np.linspace(0, nsteps, nsteps)

xcells = np.linspace(xl, xr, ncells)  # the x-coordinates of the plots

#font
#mpl.rcParams['font.family'] = 'sans-serif'
#mpl.rcParams['font.sans-serif'] = 'Fira Sans'
# mpl.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['font.serif'] ='TeXGyreSchola'
# mpl.rcParams['font.sans-serif'] ='Fira Code'
# mpl.rcParams['font.sans-serif'] = 'Century Gothic'
# mpl.rcParams['font.size'] = 12
# mpl.rcParams['legend.fontsize'] = 'medium'
#tick label spacing and tick width
# mpl.rcParams['xtick.major.pad'] = 4
# mpl.rcParams['ytick.major.pad'] = 5
# mpl.rcParams['xtick.major.width'] = 1
# mpl.rcParams['ytick.major.width'] = 1
#legend style
# mpl.rcParams['legend.frameon'] = True
# mpl.rcParams['legend.numpoints'] = 3
#mpl.rcParams['backend'] = 'PDF'
#mpl.rcParams['savefig.dpi'] = 200

def empty_marker(color):
    return {'facecolor': 'white', 'edgecolor': color, 's': circ_area, 'zorder': 2, 'linewidths': 1.5 }

def filled_marker(color):
    return {'color': color, 's': circ_area, 'zorder': 2, 'linewidths': 1.5 }

def line_empty_marker(color):
    return {'markerfacecolor': 'white', 'markeredgecolor':color, 'markersize': 4, 'markeredgewidth': 1.5 }

def line_filled_marker(color):
    return {'color': color, 'markersize': 4, 'markeredgewidth': 1.5 }

def line(color):
    return {'linestyle': '-', 'color': color, 'zorder': 1, 'linewidth': 2}

def titlestr(t):
    d = int(t / day)                 # The number of days
    h = int(int(t % day) / hour)     # The number of remaining hours
    m = int(int(t % hour) / minute)  # The number of remaining minutes
    return '{:>3d}d {:>2}h {:>2}m'.format(int(d), str(int(h)).zfill(2), str(int(m)).zfill(2))

def plot_figures_ph():

    for i in plot_at_selected_steps:
        ##print("On pH figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_ph = data_class[indx_pH]
        data_smart_ph = data_smart[indx_pH]
        #plt.axes(xlim=(-0.01,0.501))
        plt.axes(xlim=(-1, 75.1), ylim=(4, 9))
        plt.xlabel('Distance [m]')
        plt.ylabel('pH')
        plt.title(titlestr(t))
        plt.plot(xcells, data_class_ph, label='pH', **line('teal'))
        plt.plot(xcells[status[i-1]==0], data_smart_ph[status[i-1]==0], 'o', **line_empty_marker('teal'))
        plt.plot(xcells[status[i-1]==1], data_smart_ph[status[i-1]==1], 'o', **line_filled_marker('teal'))
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='lower right')
        plt.savefig(folder_general + '/pH-{}.pdf'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_pyrrhotite_siderite():

    for i in plot_at_selected_steps:
        #print("On pyrrhotite-siderite figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T

        data_class_pyrrhotite, data_class_siderite = data_class[indx_pyrrhotite], data_class[indx_siderite]
        data_smart_pyrrhotite, data_smart_siderite = data_smart[indx_pyrrhotite], data_smart[indx_siderite]
        plt.axes(xlim=(-1, 75.1), ylim=(-0.1, 5.5))
        plt.ylabel('Mineral Amount')
        plt.xlabel('Distance [m]')
        plt.title(titlestr(t))
        plt.plot(xcells, data_class_pyrrhotite, label='Pyrrhotite', **line('C0'))
        plt.plot(xcells, data_class_siderite, label='Siderite', **line('C1'))
        plt.plot(xcells[status[i-1]==0], data_smart_pyrrhotite[status[i-1]==0], 'o', **line_empty_marker('C0'))
        plt.plot(xcells[status[i-1]==1], data_smart_pyrrhotite[status[i-1]==1], 'o', **line_filled_marker('C0'))
        plt.plot(xcells[status[i-1]==0], data_smart_siderite[status[i-1]==0], 'o', **line_empty_marker('C1'))
        plt.plot(xcells[status[i-1]==1], data_smart_siderite[status[i-1]==1], 'o', **line_filled_marker('C1'))
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='center right')
        plt.savefig(folder_general + '/pyrrhotite-siderite-{}.pdf'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_aqueous_species():

    for i in plot_at_selected_steps:
        #print("On aqueous-species figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T

        data_class_Hcation  = data_class[indx_Hcation]
        data_class_HSanion  = data_class[indx_HSanion]
        data_class_S2anion = data_class[indx_S2anion]
        data_class_SO4anion     = data_class[indx_SO4anion]
        data_class_HSO4anion   = data_class[indx_HSO4anion]
        data_class_H2Saq = data_class[indx_H2Saq]

        data_smart_Hcation = data_smart[indx_Hcation]
        data_smart_HSanion = data_smart[indx_HSanion]
        data_smart_S2anion = data_smart[indx_S2anion]
        data_smart_SO4anion = data_smart[indx_SO4anion]
        data_smart_HSO4anion = data_smart[indx_HSO4anion]
        data_smart_H2Saq = data_smart[indx_H2Saq]

        plt.axes(xlim=(-1, 75.1), ylim=(1e-12, 1e1))
        plt.xlabel('Distance [m]')
        plt.ylabel('Concentration [molal]')
        plt.yscale('log')
        plt.title(titlestr(t))
        plt.plot(xcells, data_class_Hcation, label=r'$\mathrm{H^+}$', **line('darkviolet'))[0],
        plt.plot(xcells, data_class_HSanion, label=r'$\mathrm{HS^{-}}$', **line('C0'))[0],
        plt.plot(xcells, data_class_S2anion, label=r'$\mathrm{S^{2-}}$', **line('C1'))[0],
        plt.plot(xcells, data_class_SO4anion, label=r'$\mathrm{HCO_3^{-}}$',**line('C2'))[0],
        plt.plot(xcells, data_class_HSO4anion, label=r'$\mathrm{HSO_4^{-}}$',**line('red'))[0],
        plt.plot(xcells, data_class_H2Saq, label=r'$\mathrm{H2S(aq)}$', **line('gold'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_Hcation[status[i-1]==0], 'o', **line_empty_marker('darkviolet'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_Hcation[status[i-1]==1], 'o', **line_filled_marker('darkviolet'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_HSanion[status[i-1]==0], 'o', **line_empty_marker('C0'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_HSanion[status[i-1]==1], 'o', **line_filled_marker('C0'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_S2anion[status[i-1]==0], 'o', **line_empty_marker('C1'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_S2anion[status[i-1]==1], 'o', **line_filled_marker('C1'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_SO4anion[status[i-1]==0], 'o', **line_empty_marker('C2'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_SO4anion[status[i-1]==1], 'o', **line_filled_marker('C2'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_HSO4anion[status[i-1]==0], 'o', **line_empty_marker('red'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_HSO4anion[status[i-1]==1], 'o', **line_filled_marker('red'))[0],
        plt.plot(xcells[status[i - 1] == 0], data_smart_H2Saq[status[i - 1] == 0], 'o', **line_empty_marker('gold'))[0],
        plt.plot(xcells[status[i - 1] == 1], data_smart_H2Saq[status[i - 1] == 1], 'o', **line_filled_marker('gold'))[0],
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='upper right')
        plt.savefig(folder_general + '/aqueous-species-{}.pdf'.format(i))
        plt.tight_layout()
        plt.close()

def plot_animation_ph():

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-1, 75.1), ylim=(4.0, 9.0))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('pH')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label='pH', **line('teal'))[0],
        ax.plot([], [], 'o', **line_empty_marker('teal'))[0],
        ax.plot([], [], 'o', **line_filled_marker('teal'))[0],
    ]
    ax.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
    ax.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
    ax.legend(loc='lower right')

    def init():
        return tuple(objects)

    def animate(i):
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_ph = data_class[indx_pH]
        data_smart_ph = data_smart[indx_pH]
        objects[0].set_data(xcells, data_class_ph)
        objects[1].set_data(xcells[status[i]==0], data_smart_ph[status[i]==0])
        objects[2].set_data(xcells[status[i]==1], data_smart_ph[status[i]==1])
        ax.set_title(titlestr(t))
        return tuple(objects)

    print("Generating the animation of pH behaviour ...")
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)
    anim.save(folder_general + '/pH.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])
    print("Finished!")

def plot_animation_pyrrhotite_siderite():

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-1, 75.1), ylim=(-0.1, 5.5))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Mineral Amount')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label='Pyrrhotite', **line('C0'))[0],
        ax.plot([], [], label='Siderite', **line('C1'))[0],
        ax.plot([], [], 'o', **line_empty_marker('C0'))[0],
        ax.plot([], [], 'o', **line_filled_marker('C0'))[0],
        ax.plot([], [], 'o', **line_empty_marker('C1'))[0],
        ax.plot([], [], 'o', **line_filled_marker('C1'))[0],
    ]
    ax.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
    ax.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
    ax.legend(loc='center right')


    def init():
        return tuple(objects)


    def animate(i):
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_pyrrhotite, data_class_siderite = data_class[indx_pyrrhotite], data_class[indx_pyrrhotite]
        data_smart_pyrrhotite, data_smart_siderite = data_smart[indx_siderite], data_smart[indx_siderite]
        objects[0].set_data(xcells, data_class_pyrrhotite)
        objects[1].set_data(xcells, data_class_siderite)
        objects[2].set_data(xcells[status[i]==0], data_smart_pyrrhotite[status[i]==0])
        objects[3].set_data(xcells[status[i]==1], data_smart_pyrrhotite[status[i]==1])
        objects[4].set_data(xcells[status[i]==0], data_smart_siderite[status[i]==0])
        objects[5].set_data(xcells[status[i]==1], data_smart_siderite[status[i]==1])
        ax.set_title(titlestr(t))
        return tuple(objects)

    print("Generating the animation of pyrrhotite-siderite behaviour ...")
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)
    anim.save(folder_general + '/pyrrhotite-siderite.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])
    print("Finished!")

def plot_animation_aqueous_species():

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-1, 75.1), ylim=(1e-12, 1e1))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Concentration [molal]')
    ax.set_yscale('log')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label=r'$\mathrm{H^+}$', **line('C0'))[0],
        ax.plot([], [], label=r'$\mathrm{HS^{-}}$', **line('C1'))[0],
        ax.plot([], [], label=r'$\mathrm{S^{2-}}$',**line('C2'))[0],
        ax.plot([], [], label=r'$\mathrm{HCO_3^{-}}$',**line('red'))[0],
        ax.plot([], [], label=r'$\mathrm{HSO_4^{-}}$', **line('darkviolet'))[0],
        ax.plot([], [], label=r'$\mathrm{H2S(aq)}$', **line('gold'))[0],
        ax.plot([], [], 'o', **line_empty_marker('darkviolet'))[0],
        ax.plot([], [], 'o', **line_filled_marker('darkviolet'))[0],
        ax.plot([], [], 'o', **line_empty_marker('C0'))[0],
        ax.plot([], [], 'o', **line_filled_marker('C0'))[0],
        ax.plot([], [], 'o', **line_empty_marker('C1'))[0],
        ax.plot([], [], 'o', **line_filled_marker('C1'))[0],
        ax.plot([], [], 'o', **line_empty_marker('C2'))[0],
        ax.plot([], [], 'o', **line_filled_marker('C2'))[0],
        ax.plot([], [], 'o', **line_empty_marker('red'))[0],
        ax.plot([], [], 'o', **line_filled_marker('red'))[0],
        ax.plot([], [], 'o', **line_empty_marker('gold'))[0],
        ax.plot([], [], 'o', **line_filled_marker('gold'))[0],
    ]
    ax.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
    ax.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
    ax.legend(loc='upper right')

    def init():
        return tuple(objects)

    def animate(i):
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T

        data_class_Hcation = data_class[indx_Hcation]
        data_class_HSanion = data_class[indx_HSanion]
        data_class_S2anion = data_class[indx_S2anion]
        data_class_SO4anion = data_class[indx_SO4anion]
        data_class_HSO4anion = data_class[indx_HSO4anion]
        data_class_H2Saq = data_class[indx_H2Saq]

        data_smart_Hcation = data_smart[indx_Hcation]
        data_smart_HSanion = data_smart[indx_HSanion]
        data_smart_S2anion = data_smart[indx_S2anion]
        data_smart_SO4anion = data_smart[indx_SO4anion]
        data_smart_HSO4anion = data_smart[indx_HSO4anion]
        data_smart_H2Saq = data_smart[indx_H2Saq]

        objects[0].set_data(xcells, data_class_HSanion)
        objects[1].set_data(xcells, data_class_S2anion)
        objects[2].set_data(xcells, data_class_SO4anion)
        objects[3].set_data(xcells, data_class_HSO4anion)
        objects[4].set_data(xcells, data_class_H2Saq)
        objects[5].set_data(xcells, data_class_Hcation)

        objects[6].set_data(xcells[status[i]==0], data_smart_HSanion[status[i]==0])
        objects[7].set_data(xcells[status[i]==1], data_smart_HSanion[status[i]==1])

        objects[8].set_data(xcells[status[i]==0], data_smart_S2anion[status[i]==0])
        objects[9].set_data(xcells[status[i]==1], data_smart_S2anion[status[i]==1])

        objects[10].set_data(xcells[status[i]==0], data_smart_SO4anion[status[i]==0])
        objects[11].set_data(xcells[status[i] == 0], data_smart_SO4anion[status[i] == 0])

        objects[12].set_data(xcells[status[i]==1], data_smart_HSO4anion[status[i]==1])
        objects[13].set_data(xcells[status[i]==0], data_smart_HSO4anion[status[i]==0])

        objects[14].set_data(xcells[status[i]==1], data_smart_H2Saq[status[i]==1])
        objects[15].set_data(xcells[status[i]==0], data_smart_H2Saq[status[i]==0])

        objects[16].set_data(xcells[status[i]==1], data_smart_Hcation[status[i]==1])
        objects[17].set_data(xcells[status[i] == 1], data_smart_Hcation[status[i] == 1])

        ax.set_title(titlestr(t))
        return tuple(objects)

    print("Generating the animation of aqueous species behaviour ...")
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)
    anim.save(folder_general + '/aqueous-species.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])
    print("Finished!")


def plot_computing_costs():

    step = 20
    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    with open(folder_smart + '/analysis-smart.json') as read_file:
        data_smart = json.load(read_file)
    with open(folder_class + '/analysis-conventional.json') as read_file:
        data_class = json.load(read_file)

    timing_class = data_class.get('computing_costs_per_time_step')
    timing_smart = data_smart.get('computing_costs_per_time_step')

    timings_transport = np.array(timing_class.get('transport')) * 1e6  # in microseconds
    timings_equilibrium_class = np.array(timing_class.get('equilibrium')) * 1e6  # in microseconds
    timings_equilibrium_smart = np.array(timing_smart.get('smart_equilibrium')) * 1e6  # in microseconds

    timings_equilibrium_smart_ideal = (np.array(timing_smart.get('smart_equilibrium'))
                                       - np.array(timing_smart.get('smart_equilibrium_nearest_neighbor_search'))) * 1e6  # in microseconds

    plt.xlabel('Time Step')
    plt.ylabel('Computing Cost [μs]')
    plt.xlim(left=0, right=nsteps)
    plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], timings_equilibrium_class[0:nsteps:step], label="Chemical Equilibrium (Conventional)", color='C0', linewidth=2)
    plt.plot(time_steps[0:nsteps:step], timings_equilibrium_smart[0:nsteps:step], label="Chemical Equilibrium (Smart)", color='C1', linewidth=2, alpha=1.0)
    plt.plot(time_steps[0:nsteps:step], timings_equilibrium_smart_ideal[0:nsteps:step], label="Chemical Equilibrium (Smart)", color='C3', linewidth=2, alpha=1.0)
    plt.plot(time_steps[0:nsteps:step], timings_transport[0:nsteps:step], label="Transport", color='C2', linewidth=2, alpha=1.0)
    leg = plt.legend(loc='lower right', bbox_to_anchor=(1, 0.13))
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/computing-costs-nolegend-with-smart-ideal.pdf')
    #plt.savefig(folder_general + '/computing-costs-nolegend-with-smart.pdf')
    #plt.savefig(folder_general + '/computing-costs-nolegend.pdf')

    plt.close()


def plot_on_demand_learning_countings():

    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    with open(folder_smart + '/analysis-smart.json') as read_file:
        data = json.load(read_file)

    status_learnings = data.get('cells_where_learning_was_required_at_step')

    # Collect the number of learnings on each step
    learnings = [len(x) for x in status_learnings]

    plt.xlabel('Time Step')
    plt.xlim(left=0, right=nsteps)
    plt.ylim(bottom=0, top=np.max(learnings)+1)
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(learnings, color='C0', linewidth=2)
    plt.tight_layout()
    plt.savefig(folder_general + '/on-demand-learning-countings.pdf')
    plt.close()


def plot_on_demand_learnings_total():

    step = 1

    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    with open(folder_smart + '/analysis-smart.json') as read_file:
        data = json.load(read_file)

    status_learnings = data.get('cells_where_learning_was_required_at_step')

    # Collect the number of learnings on each step and accumulate them w.r.t. time-step
    learnings = [len(x) for x in status_learnings]
    accum_learnings = [np.sum(learnings[0:i]) for i in range(0, nsteps)]


    plt.xlabel('Time Step')
    plt.xlim(left=0, right=nsteps)
    plt.ylim(bottom=0, top=np.max(accum_learnings)+5)
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], accum_learnings[0:nsteps:step], color='C0', linewidth=2)
    plt.tight_layout()
    plt.savefig(folder_general + '/on-demand-learning-total.pdf')
    plt.close()


def plot_computing_costs_vs_total_learnings():

    step = 20
    with open(folder_smart + '/analysis-smart.json') as read_file:
        data_smart = json.load(read_file)

    timing_smart = data_smart.get('computing_costs_per_time_step')

    timings_estimate = np.array(timing_smart.get('smart_equilibrium_estimate')) * 1e6       # in microseconds
    timings_search = np.array(timing_smart.get('smart_equilibrium_nearest_neighbor_search')) * 1e6  # in microseconds`
    #timings_taylor = (np.array(timing_smart.get('smart_equilibrium_estimate')) \
    #                  - np.array(timing_smart.get('smart_equilibrium_acceptance')) \
    #                  - np.array(timing_smart.get('smart_equilibrium_nearest_neighbor_search'))) * 1e6  # in microseconds`

    plt.xlabel('Time Step')
    plt.ylabel('Computing Cost [μs]')
    #plt.xscale('log')
    plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], timings_search[0:nsteps:step], label="Search Time", color='C0', linewidth=2)
    plt.plot(time_steps[0:nsteps:step], timings_estimate[0:nsteps:step], label="Estimate Time", color='C1', linewidth=2)

    leg = plt.legend(loc='lower right')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/search-traylor-vs-total-learnings.pdf')
    plt.close()


def plot_speedups():

    step = 80
    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    with open(folder_smart + '/analysis-smart.json') as read_file:
        data_smart = json.load(read_file)
    with open(folder_class + '/analysis-conventional.json') as read_file:
        data_class = json.load(read_file)

    timing_class = data_class.get('computing_costs_per_time_step')
    timing_smart = data_smart.get('computing_costs_per_time_step')

    speedup = np.array(timing_class.get('equilibrium')) / np.array(timing_smart.get('smart_equilibrium'))
    speedup_ideal = np.array(timing_class.get('equilibrium')) / \
                    (np.array(timing_smart.get('smart_equilibrium'))
                     - np.array(timing_smart.get('smart_equilibrium_nearest_neighbor_search'))
                     - np.array(timing_smart.get('smart_equilibrium_storage')))

    print("speedup         = ", np.sum(np.array(timing_class.get('equilibrium'))) / np.sum(np.array(timing_smart.get('smart_equilibrium'))))
    print("speedup (ideal) = ", np.sum(np.array(timing_class.get('equilibrium'))) / (np.sum(np.array(timing_smart.get('smart_equilibrium')) \
                                                                                            - np.array(timing_smart.get('smart_equilibrium_nearest_neighbor_search')) \
                                                                                            - np.array(timing_smart.get('smart_equilibrium_storage')))))


    plt.xlim(left=0, right=nsteps)
    plt.xlabel('Time Step')
    plt.ylabel('Speedup (-)')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], speedup[0:nsteps:step], label="Conventional vs. Smart ", color='C0', linewidth=2)
    plt.plot(time_steps[0:nsteps:step], speedup_ideal[0:nsteps:step], label="Conventional vs. Smart (Ideal search)",
             color='C2', linewidth=2, alpha=1.0)

    leg = plt.legend(loc='upper right')
    for line in leg.get_lines():
       line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/speedups.pdf')
    #plt.show()
    #plt.savefig(folder_general + '/speedups-nolegend-withupperbound.pdf')
    #plt.savefig(folder_general + '/speedups-nolegend.pdf')
    plt.close()

def count_trainings(status):
    counter = [st for st in status if st == 0]
    return (len(counter), len(status))


def calculate_error(tol):

    np.set_printoptions(precision=3)

    for i in plot_at_selected_steps:

        filearray_class = np.loadtxt(folder_class + '/' + files_class[i - 1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i - 1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T

        # Fetch reference and approximate data
        data_class_ph, data_smart_ph = data_class[indx_ph], data_smart[indx_ph]
        data_class_calcite, data_smart_calcite = data_class[indx_calcite], data_smart[indx_calcite]
        data_class_dolomite, data_smart_dolomite = data_class[indx_dolomite], data_smart[indx_dolomite]
        data_class_cacation, data_smart_cacation = data_class[indx_Cacation], data_smart[indx_Cacation]
        data_class_mgcation, data_smart_mgcation = data_class[indx_Mgcation], data_smart[indx_Mgcation]
        data_class_hco3anion, data_smart_hco3anion = data_class[indx_HCO3anion], data_smart[indx_HCO3anion]
        data_class_co2aq, data_smart_co2aq = data_class[indx_CO2aq], data_smart[indx_CO2aq]
        data_class_hcation, data_smart_hcation = data_class[indx_Hcation], data_smart[indx_Hcation]
        data_class_CO3anion, data_smart_CO3anion = data_class[indx_CO3anion], data_smart[indx_CO3anion]
        data_class_CaClcation, data_smart_CaClcation = data_class[indx_CaClcation], data_smart[indx_CaClcation]
        data_class_CaHCO3cation, data_smart_CaHCO3cation = data_class[indx_CaHCO3cation], data_smart[indx_CaHCO3cation]
        data_class_MgClcation, data_smart_MgClcation = data_class[indx_MgClcation], data_smart[indx_MgClcation]
        data_class_MgHCO3caanion, data_smart_MgHCO3caanion = data_class[indx_MgHCO3caanion], data_smart[indx_MgHCO3caanion]
        data_class_OHanion, data_smart_OHanion = data_class[indx_OHanion], data_smart[indx_OHanion]

        abs_error_ph = np.abs(data_class_ph - data_smart_ph)
        abs_error_calcite = np.abs(data_class_calcite - data_smart_calcite)
        abs_error_dolomite = np.abs(data_class_dolomite - data_smart_dolomite)
        abs_error_cacation = np.abs(data_class_cacation - data_smart_cacation)
        abs_error_mgcation = np.abs(data_class_mgcation - data_smart_mgcation)
        abs_error_hco3anion = np.abs(data_class_hco3anion - data_smart_hco3anion)
        abs_error_co2aq = np.abs(data_class_co2aq - data_smart_co2aq)
        abs_error_hcation = np.abs(data_class_hcation - data_smart_hcation)
        abs_error_CO3anion = np.abs(data_class_CO3anion - data_smart_CO3anion)
        abs_error_CaHCO3cation = np.abs(data_class_CaHCO3cation - data_smart_CaHCO3cation)
        abs_error_CaClcation = np.abs(data_class_CaClcation - data_smart_CaClcation)
        abs_error_MgClcation = np.abs(data_class_MgClcation - data_smart_MgClcation)
        abs_error_MgHCO3caanion = np.abs(data_class_MgHCO3caanion - data_smart_MgHCO3caanion)
        abs_error_OHanion = np.abs(data_class_OHanion - data_smart_OHanion)

        data_class_dolomite[data_class_dolomite<=tol] = 1.0
        data_class_calcite[data_class_calcite<=tol] = 1.0
        data_class_cacation[data_class_cacation<=tol] = 1.0
        data_class_ph[data_class_ph<=tol] = 1.0
        data_class_mgcation[data_class_mgcation<=tol] = 1.0
        data_class_hco3anion[data_class_hco3anion<=tol] = 1.0
        data_class_co2aq[data_class_co2aq<=tol] = 1.0
        data_class_hcation[data_class_hcation<=tol] = 1.0
        data_class_CO3anion[data_class_CO3anion<=tol] = 1.0
        data_class_CaHCO3cation[data_class_CaHCO3cation<=tol] = 1.0
        data_class_CaClcation[data_class_CaClcation<=tol] = 1.0
        data_class_MgClcation[data_class_MgClcation<=tol] = 1.0
        data_class_MgHCO3caanion[data_class_MgHCO3caanion<=tol] = 1.0
        data_class_OHanion[data_class_OHanion<=tol] = 1.0

        rel_error_ph = np.divide(abs_error_ph, data_class_ph)
        rel_error_calcite = np.divide(abs_error_calcite, data_class_calcite)
        rel_error_dolomite = np.divide(abs_error_dolomite, data_class_dolomite)
        rel_error_cacation = np.divide(abs_error_cacation, data_class_cacation)
        rel_error_mgcation = np.divide(abs_error_mgcation, data_class_mgcation)
        rel_error_hco3anion = np.divide(abs_error_hco3anion, data_class_hco3anion)
        rel_error_co2aq = np.divide(abs_error_co2aq, data_class_co2aq)
        rel_error_hcation = np.divide(abs_error_hcation, data_class_hcation)
        rel_error_CO3anion = np.divide(abs_error_CO3anion, data_class_CO3anion)
        rel_error_CaHCO3cation = np.divide(abs_error_CaHCO3cation, data_class_CaHCO3cation)
        rel_error_CaClcation = np.divide(abs_error_CaClcation, data_class_CaClcation)
        rel_error_MgClcation = np.divide(abs_error_MgClcation, data_class_MgClcation)
        rel_error_MgHCO3caanion = np.divide(abs_error_MgHCO3caanion, data_class_MgHCO3caanion)
        rel_error_OHanion = np.divide(abs_error_OHanion, data_class_OHanion)

        t = i * dt


        plt.axes(xlim=(-0.01, 0.501))
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('Relative Error [-]', fontproperties=prop)
        plt.title(titlestr(t), fontproperties=prop, fontsize=20)
        plt.plot(xcells, rel_error_calcite, label='Calcite', **line_error('C0'))
        plt.plot(xcells, rel_error_dolomite, label='Dolomite', **line_error('C1'))
        prop.set_size(12)
        plt.legend(loc='center right', prop=prop)
        plt.savefig(folder_general + f'/error-rel-calcite-dolomite-long-{i}.pdf')
        plt.tight_layout()
        plt.close()

        plt.axes(xlim=(-0.01, 0.501))
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('Relative Error [-]', fontproperties=prop)
        plt.yscale('log')
        plt.title(titlestr(t), fontproperties=prop, fontsize=20)
        plt.plot(xcells, rel_error_cacation, label=r'$Ca^{2+}$', **line_error('C0'))
        plt.plot(xcells, rel_error_mgcation, label=r'$Mg^{2+}$', **line_error('C1'))
        plt.plot(xcells, rel_error_hco3anion, label=r'$HCO_3^-$', **line_error('C2'))
        plt.plot(xcells, rel_error_co2aq, label=r'$CO_2(aq)$', **line_error('C3'))
        plt.plot(xcells, rel_error_hcation, label=r'$H^+$', **line_error('C4'))
        plt.plot(xcells, rel_error_CO3anion, label=r'$CO_3^{-2}$', **line_error('C5'))[0],
        plt.plot(xcells, rel_error_CaHCO3cation, label=r'$Ca(HCO_3)^+$', **line_error('C6'))[0],
        plt.plot(xcells, rel_error_CaClcation, label=r'$CaCl^+$',**line_error('C7'))[0],
        plt.plot(xcells, rel_error_MgClcation, label=r'$MgCl^+$',**line_error('C8'))[0],
        plt.plot(xcells, rel_error_MgHCO3caanion, label=r'$Mg(HCO_3))^+$',**line_error('C9'))[0],
        plt.plot(xcells, rel_error_OHanion, label=r'$OH^-$', **line_error('darkviolet'))[0],
        prop.set_size(12)
        plt.legend(loc='upper right', prop=prop)
        plt.savefig(folder_general + f'/error-rel-aqueous-species-long-{i}.pdf')
        plt.tight_layout()
        plt.close()

def mass_conservation():

    indices = [0, 19, 2399]
    for i in indices:

        filearray_res = np.loadtxt(folder_smart + f'/res-{i}.txt')
        filearray_state = np.loadtxt(folder_smart + f'/test-{i}.txt', skiprows=1)

        filearray_res = filearray_res.T
        filearray_state = filearray_state.T

        data_C = np.abs(filearray_state[indx_C])
        data_Ca = np.abs(filearray_state[indx_Ca])
        data_Cl = np.abs(filearray_state[indx_Cl])
        data_H = np.abs(filearray_state[indx_H])
        data_Mg = np.abs(filearray_state[indx_Mg])
        data_Na = np.abs(filearray_state[indx_Na])
        data_O = np.abs(filearray_state[indx_O])
        data_Si = np.abs(filearray_state[indx_Si])
        data_Z = np.abs(filearray_state[indx_Z])

        # sigma(Ca) = abs(Ca) if abs(Ca) < eps_b else 1.0
        eps_b = 1e-16
        data_C[data_C<=eps_b] = 1.0
        data_Ca[data_Ca<=eps_b] = 1.0
        data_Cl[data_Cl<=eps_b] = 1.0
        data_H[data_H<=eps_b] = 1.0
        data_Mg[data_Mg<=eps_b] = 1.0
        data_Na[data_Na<=eps_b] = 1.0
        data_O[data_O<=eps_b] = 1.0
        data_Si[data_Si<=eps_b] = 1.0
        data_Z[data_Z<=eps_b] = 1.0

        res_shift = 14
        data_res_C = np.abs(filearray_res[indx_C-res_shift])
        data_res_Ca = np.abs(filearray_res[indx_Ca-res_shift])
        data_res_Cl = np.abs(filearray_res[indx_Cl-res_shift])
        data_res_H = np.abs(filearray_res[indx_H-res_shift])
        data_res_Mg = np.abs(filearray_res[indx_Mg-res_shift])
        data_res_Na = np.abs(filearray_res[indx_Na-res_shift])
        data_res_O = np.abs(filearray_res[indx_O-res_shift])
        data_res_Si = np.abs(filearray_res[indx_Si-res_shift])

        # error[Ca] = abs(r[Ca])/sigma(Ca)
        error_C = np.divide(data_res_C, data_C)
        error_Ca = np.divide(data_res_Ca, data_Ca)
        error_Cl = np.divide(data_res_Cl, data_Cl)
        error_H = np.divide(data_res_H, data_H)
        error_Mg = np.divide(data_res_Mg, data_Mg)
        error_Na = np.divide(data_res_Na, data_Na)
        error_O = np.divide(data_res_O, data_O)
        error_Si = np.divide(data_res_Si, data_Si)

        #plt.axes(xlim=(-0.01, 0.501), ylim=(1e-16, 5e-7))
        plt.axes(xlim=(-0.01, 0.501))
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('Mass Balance Error [%]', fontproperties=prop)
        plt.yscale('log')
        plt.title(titlestr(dt * (i+1)), fontproperties=prop,fontsize=20)
        plt.plot(xcells, 100 * error_C, label='C', **line_error('C0'))
        plt.plot(xcells, 100 * error_Ca, label='Ca', **line_error('C1'))
        plt.plot(xcells, 100 * error_Cl, label='Cl', **line_error('C2'))
        plt.plot(xcells, 100 * error_H, label='H', **line_error('C3'))
        plt.plot(xcells, 100 * error_Mg, label='Mg', **line_error('C4'))
        plt.plot(xcells, 100 * error_Na, label='Na', **line_error('C5'))
        plt.plot(xcells, 100 * error_O, label='O', **line_error('C6'))
        plt.plot(xcells, 100 * error_Si, label='Si', **line_error('C7'))

        plt.legend(loc='lower right', prop=prop)
        plt.savefig(folder_general + f'/mass-balance-residual-{i}.pdf')
        plt.tight_layout()
        #plt.show()
        plt.close()

if __name__ == '__main__':

    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    with open(folder_smart + '/analysis-smart.json') as read_file:
        status_data = json.load(read_file)

    status_learnings = status_data.get('cells_where_learning_was_required_at_step')
    status = np.ones([nsteps, ncells])
    for i in range(0, nsteps):
        for j in range(0, len(status_learnings[i])):
            status[i][j] = 0


    # Count the percentage of the trainings needed
    training_counter = count_trainings(np.array(status).flatten())
    title = "%2.2f percent is training (%d out of %d cells)" % (
        100 * training_counter[0] / training_counter[1], training_counter[0], training_counter[1])
    text = ["Number of chemical equilibrium trainings: %d" % training_counter[0],
            "Number of smart chemical equilibrium estimations: %d" % (training_counter[1] - training_counter[0]),
            ("Percentage of smart chemical kineics predictions: %2.2f" % (
                    100 - 100 * training_counter[0] / training_counter[1])) + "$\%$"]
    print(title)

    print("Collecting files...")
    # Collect files with results corresponding to smart or reference (classical) solver
    files_smart = [file for file in natsorted( os.listdir(folder_smart) ) if ("analysis" not in file)]
    files_class = [file for file in natsorted( os.listdir(folder_class) ) if ("analysis" not in file)]

    plot_on_demand_learning_countings()
    plot_on_demand_learnings_total()
    plot_speedups()
    #plot_computing_costs_vs_total_learnings()
    plot_computing_costs()
    plot_figures_ph()
    plot_figures_aqueous_species()
    plot_figures_pyrrhotite_siderite()
    #plot_figures_pyrrhotite_siderite_moles()

    from matplotlib import animation

    animation_starts_at_frame = 0      # the first frame index to be considered
    animation_ends_at_frame = 10 * 500  # the last frame index to be considered
    animation_num_frames_to_jump = 10     # the number of frames to jump between current and next
    # Check for the correct end frame number
    assert animation_ends_at_frame <= nsteps, "WARNING: The number of the end frame must be smaller then number of steps! "

    # Provide the number of frames per second and the time (in milliseconds) to wait between each frame:

    animation_fps = 30 # the number of frames per second
    animation_interval_wait = 200    # the time (in milliseconds) to wait between each frame
    # Auxiliary animation options
    animation_frame_range = range(animation_starts_at_frame, animation_ends_at_frame, animation_num_frames_to_jump)

    #plot_animation_ph()
    plot_animation_pyrrhotite_siderite()
    plot_animation_aqueous_species()