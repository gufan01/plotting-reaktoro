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

from progress.bar import IncrementalBar

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
folder_smart = 'results-scaveging-custering-primary-species-dt-4320-ncells-100-nsteps-5000-reltol-1.0e-02-dk-smart'
folder_class = 'results-scaveging-custering-primary-species-dt-4320-ncells-100-nsteps-5000-dk-reference'
folder_general = 'plots-scaveging-custering-primary-species-dt-4320-ncells-100-nsteps-5000-reltol-1.0e-02-dk-smart'
#folder_smart = 'cpp-reactivetransport-old-demo/results-pitzer-dt-3594.24-ncells-100-nsteps-10000-eqreltol-1.0e-01-eqabstol-1.0e-08-smart'
#folder_class = 'cpp-reactivetransport-old-demo/results-pitzer-dt-3594.24-ncells-100-nsteps-10000-reference'
#folder_general = "plots-results-pitzer-dt-3594.24-ncells-100-nsteps-10000"
os.system('mkdir -p ' + folder_general)

fillz = len(str(123))

import os
from matplotlib import font_manager as fm, rcParams
#fpath = os.path.join(rcParams["datapath"], "/usr/share/fonts/truetype/ebgaramond/EBGaramond12-Regular.ttf")
#fpath = os.path.join(rcParams["datapath"], "/usr/share/fonts/truetype/ebgaramond/EBGaramond12-Bold.ttf")
#fpath = os.path.join(rcParams["datapath"], "/usr/share/fonts/truetype/ebgaramond/TeX-Gyre-Adventor.ttf")
fpath = os.path.join(rcParams["datapath"], "/home/skyas/Dropbox/texgyreadventor-regular.otf")
prop = fm.FontProperties(fname=fpath)
prop.set_size(14)

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'TeX Gyre Adventor'
mpl.rcParams['font.style'] = 'normal'
mpl.rcParams['font.size'] = 14

# Indices of the loaded data to plot
indx_pH = 0
indx_Hcation = 1
indx_HSanion = 2
indx_S2anion = 3
indx_CO3anion = 4
indx_HSO4anion = 5
indx_H2Saq = 6
indx_phase_amount_pyrrhotite = 7
indx_phase_amount_siderite = 8
indx_phase_vol_pyrrhotite = 9
indx_phase_vol_siderite = 10
indx_C = 11
indx_Ca = 12
indx_Cl = 13
indx_Fe = 14
indx_H = 15
indx_K = 16
indx_Mg = 17
indx_Na = 18
indx_O = 19
indx_S = 20
indx_Z = 21

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
#mpl.rcParams['backend'] = 'png'
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

def line_error(color):
    return {'linestyle': ':', 'marker': 'D', 'color': color, 'zorder': 1, 'linewidth': 0.5, 'markersize': 4}

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
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('pH', fontproperties=prop)
        plt.title(titlestr(t), fontproperties=prop,fontsize=20)
        plt.plot(xcells, data_class_ph, label='pH', **line('teal'))
        plt.plot(xcells[status[i-1]==0], data_smart_ph[status[i-1]==0], 'o', **line_empty_marker('teal'))
        plt.plot(xcells[status[i-1]==1], data_smart_ph[status[i-1]==1], 'o', **line_filled_marker('teal'))
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='lower right', prop=prop)
        plt.savefig(folder_general + '/pH-{}.png'.format(i))
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

        data_class_pyrrhotite, data_class_siderite = data_class[indx_phase_amount_pyrrhotite], data_class[indx_phase_amount_siderite]
        data_smart_pyrrhotite, data_smart_siderite = data_smart[indx_phase_amount_pyrrhotite], data_smart[indx_phase_amount_siderite]
        plt.axes(xlim=(-1, 75.1), ylim=(-0.1, 5.5))
        plt.ylabel('Mineral Amount', fontproperties=prop)
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.title(titlestr(t), fontproperties=prop,fontsize=20)
        plt.plot(xcells, data_class_pyrrhotite, label='Pyrrhotite', **line('C0'))
        plt.plot(xcells, data_class_siderite, label='Siderite', **line('C1'))
        plt.plot(xcells[status[i-1]==0], data_smart_pyrrhotite[status[i-1]==0], 'o', **line_empty_marker('C0'))
        plt.plot(xcells[status[i-1]==1], data_smart_pyrrhotite[status[i-1]==1], 'o', **line_filled_marker('C0'))
        plt.plot(xcells[status[i-1]==0], data_smart_siderite[status[i-1]==0], 'o', **line_empty_marker('C1'))
        plt.plot(xcells[status[i-1]==1], data_smart_siderite[status[i-1]==1], 'o', **line_filled_marker('C1'))
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='center right', prop=prop)
        plt.savefig(folder_general + '/pyrrhotite-siderite-{}.png'.format(i))
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
        data_class_CO3anion     = data_class[indx_CO3anion]
        data_class_HSO4anion   = data_class[indx_HSO4anion]
        data_class_H2Saq = data_class[indx_H2Saq]

        data_smart_Hcation = data_smart[indx_Hcation]
        data_smart_HSanion = data_smart[indx_HSanion]
        data_smart_S2anion = data_smart[indx_S2anion]
        data_smart_CO3anion = data_smart[indx_CO3anion]
        data_smart_HSO4anion = data_smart[indx_HSO4anion]
        data_smart_H2Saq = data_smart[indx_H2Saq]

        plt.axes(xlim=(-1, 75.1), ylim=(1e-12, 1e0))
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('Concentration [molal]', fontproperties=prop)
        plt.yscale('log')
        plt.title(titlestr(t), fontproperties=prop,fontsize=20)
        plt.plot(xcells, data_class_Hcation, label=r'$\mathrm{H^+}$', **line('darkviolet'))[0],
        plt.plot(xcells, data_class_HSanion, label=r'$\mathrm{HS^{-}}$', **line('C0'))[0],
        plt.plot(xcells, data_class_S2anion, label=r'$\mathrm{S^{2-}}$', **line('C1'))[0],
        plt.plot(xcells, data_class_CO3anion, label=r'$\mathrm{CO_3^{2-}}$',**line('C2'))[0],
        plt.plot(xcells, data_class_HSO4anion, label=r'$\mathrm{HSO_4^{-}}$',**line('red'))[0],
        plt.plot(xcells, data_class_H2Saq, label=r'$\mathrm{H_2S(aq)}$', **line('gold'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_Hcation[status[i-1]==0], 'o', **line_empty_marker('darkviolet'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_Hcation[status[i-1]==1], 'o', **line_filled_marker('darkviolet'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_HSanion[status[i-1]==0], 'o', **line_empty_marker('C0'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_HSanion[status[i-1]==1], 'o', **line_filled_marker('C0'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_S2anion[status[i-1]==0], 'o', **line_empty_marker('C1'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_S2anion[status[i-1]==1], 'o', **line_filled_marker('C1'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_CO3anion[status[i-1]==0], 'o', **line_empty_marker('C2'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_CO3anion[status[i-1]==1], 'o', **line_filled_marker('C2'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_HSO4anion[status[i-1]==0], 'o', **line_empty_marker('red'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_HSO4anion[status[i-1]==1], 'o', **line_filled_marker('red'))[0],
        plt.plot(xcells[status[i - 1] == 0], data_smart_H2Saq[status[i - 1] == 0], 'o', **line_empty_marker('gold'))[0],
        plt.plot(xcells[status[i - 1] == 1], data_smart_H2Saq[status[i - 1] == 1], 'o', **line_filled_marker('gold'))[0],
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='upper right', prop=prop)
        plt.savefig(folder_general + '/aqueous-species-{}.png'.format(i))
        plt.tight_layout()
        plt.close()

def plot_animation_ph():

    bar = IncrementalBar('On pH animation:', max = len(animation_frame_range))

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-1, 75.1), ylim=(4.0, 9.0))
    ax.set_xlabel('Distance [m]', fontproperties=prop)
    ax.set_ylabel('pH', fontproperties=prop)
    ax.set_title(titlestr(0.0), fontproperties=prop,fontsize=20)
    objects = [
        ax.plot([], [], label='pH', **line('teal'))[0],
        ax.plot([], [], 'o', **line_empty_marker('teal'))[0],
        ax.plot([], [], 'o', **line_filled_marker('teal'))[0],
    ]
    ax.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
    ax.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
    ax.legend(loc='lower right', prop=prop)

    def init():
        return tuple(objects)

    def animate(i):
        bar.next()
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
        ax.set_title(titlestr(t), fontproperties=prop,fontsize=20)
        return tuple(objects)

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)
    anim.save(folder_general + '/pH.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])
    bar.finish()

def plot_animation_pyrrhotite_siderite():

    bar = IncrementalBar('On pyrrhotite-siderite animation:', max = len(animation_frame_range))

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-1, 75.1), ylim=(-0.1, 5.5))
    ax.set_xlabel('Distance [m]', fontproperties=prop)
    ax.set_ylabel('Mineral Amount', fontproperties=prop)
    ax.set_title(titlestr(0.0), fontproperties=prop,fontsize=20)
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
    ax.legend(loc='center right', prop=prop)


    def init():
        return tuple(objects)


    def animate(i):
        bar.next()
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_pyrrhotite, data_class_siderite = data_class[indx_phase_amount_pyrrhotite], data_class[indx_phase_amount_pyrrhotite]
        data_smart_pyrrhotite, data_smart_siderite = data_smart[indx_phase_amount_siderite], data_smart[indx_phase_amount_siderite]
        objects[0].set_data(xcells, data_class_pyrrhotite)
        objects[1].set_data(xcells, data_class_siderite)
        objects[2].set_data(xcells[status[i]==0], data_smart_pyrrhotite[status[i]==0])
        objects[3].set_data(xcells[status[i]==1], data_smart_pyrrhotite[status[i]==1])
        objects[4].set_data(xcells[status[i]==0], data_smart_siderite[status[i]==0])
        objects[5].set_data(xcells[status[i]==1], data_smart_siderite[status[i]==1])
        ax.set_title(titlestr(t), fontproperties=prop,fontsize=20)
        return tuple(objects)

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)
    anim.save(folder_general + '/pyrrhotite-siderite.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])
    bar.finish()

def plot_animation_aqueous_species():

    bar = IncrementalBar('On aqueous species animation:', max = len(animation_frame_range))

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-1, 75.1), ylim=(1e-12, 1e0))
    ax.set_xlabel('Distance [m]', fontproperties=prop)
    ax.set_ylabel('Concentration [molal]', fontproperties=prop)
    ax.set_yscale('log')
    ax.set_title(titlestr(0.0), fontproperties=prop,fontsize=20)
    objects = [
        ax.plot([], [], label=r'$\mathrm{HS^{-}}$', **line('C1'))[0],
        ax.plot([], [], label=r'$\mathrm{S^{2-}}$',**line('C2'))[0],
        ax.plot([], [], label=r'$\mathrm{CO_3^{2-}}$',**line('red'))[0],
        ax.plot([], [], label=r'$\mathrm{HSO_4^{-}}$', **line('darkviolet'))[0],
        ax.plot([], [], label=r'$\mathrm{H_2S(aq)}$', **line('gold'))[0],
        ax.plot([], [], label=r'$\mathrm{H^+}$', **line('C0'))[0],
         ax.plot([], [], 'o', **line_empty_marker('C1'))[0],
        ax.plot([], [], 'o', **line_filled_marker('C1'))[0],
        ax.plot([], [], 'o', **line_empty_marker('C2'))[0],
        ax.plot([], [], 'o', **line_filled_marker('C2'))[0],
        ax.plot([], [], 'o', **line_empty_marker('red'))[0],
        ax.plot([], [], 'o', **line_filled_marker('red'))[0],
        ax.plot([], [], 'o', **line_empty_marker('darkviolet'))[0],
        ax.plot([], [], 'o', **line_filled_marker('darkviolet'))[0],
        ax.plot([], [], 'o', **line_empty_marker('gold'))[0],
        ax.plot([], [], 'o', **line_filled_marker('gold'))[0],
        ax.plot([], [], 'o', **line_empty_marker('C0'))[0],
        ax.plot([], [], 'o', **line_filled_marker('C0'))[0],
    ]
    ax.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
    ax.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
    ax.legend(loc='upper right', prop=prop)

    def init():
        return tuple(objects)

    def animate(i):
        bar.next()
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T

        data_class_Hcation = data_class[indx_Hcation]
        data_class_HSanion = data_class[indx_HSanion]
        data_class_S2anion = data_class[indx_S2anion]
        data_class_CO3anion = data_class[indx_CO3anion]
        data_class_HSO4anion = data_class[indx_HSO4anion]
        data_class_H2Saq = data_class[indx_H2Saq]

        data_smart_Hcation = data_smart[indx_Hcation]
        data_smart_HSanion = data_smart[indx_HSanion]
        data_smart_S2anion = data_smart[indx_S2anion]
        data_smart_CO3anion = data_smart[indx_CO3anion]
        data_smart_HSO4anion = data_smart[indx_HSO4anion]
        data_smart_H2Saq = data_smart[indx_H2Saq]

        objects[0].set_data(xcells, data_class_HSanion)
        objects[1].set_data(xcells, data_class_S2anion)
        objects[2].set_data(xcells, data_class_CO3anion)
        objects[3].set_data(xcells, data_class_HSO4anion)
        objects[4].set_data(xcells, data_class_H2Saq)
        objects[5].set_data(xcells, data_class_Hcation)

        objects[6].set_data(xcells[status[i]==0], data_smart_HSanion[status[i]==0])
        objects[7].set_data(xcells[status[i]==1], data_smart_HSanion[status[i]==1])

        objects[8].set_data(xcells[status[i]==0], data_smart_S2anion[status[i]==0])
        objects[9].set_data(xcells[status[i]==1], data_smart_S2anion[status[i]==1])

        objects[10].set_data(xcells[status[i]==0], data_smart_CO3anion[status[i]==0])
        objects[11].set_data(xcells[status[i]==1], data_smart_CO3anion[status[i]==1])

        objects[12].set_data(xcells[status[i]==0], data_smart_HSO4anion[status[i]==0])
        objects[13].set_data(xcells[status[i]==1], data_smart_HSO4anion[status[i]==1])

        objects[14].set_data(xcells[status[i]==0], data_smart_H2Saq[status[i]==0])
        objects[15].set_data(xcells[status[i]==1], data_smart_H2Saq[status[i]==1])

        objects[16].set_data(xcells[status[i]==0], data_smart_Hcation[status[i]==0])
        objects[17].set_data(xcells[status[i]==1], data_smart_Hcation[status[i]==1])

        ax.set_title(titlestr(t), fontproperties=prop,fontsize=20)
        return tuple(objects)

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)
    anim.save(folder_general + '/aqueous-species.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])
    bar.finish()


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
                                       - np.array(timing_smart.get('smart_equilibrium_search'))) * 1e6  # in microseconds

    plt.xlabel('Time Step', fontproperties=prop)
    plt.ylabel('Computing Cost [μs]', fontproperties=prop)
    plt.xlim(left=0, right=nsteps)
    plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], timings_equilibrium_class[0:nsteps:step], label="Chemical Equilibrium (Conventional)", color='C0', linewidth=2)
    plt.plot(time_steps[0:nsteps:step], timings_equilibrium_smart[0:nsteps:step], label="Chemical Equilibrium (Smart)", color='C1', linewidth=2, alpha=1.0)
    plt.plot(time_steps[0:nsteps:step], timings_equilibrium_smart_ideal[0:nsteps:step], label="Chemical Equilibrium (Smart)", color='C3', linewidth=2, alpha=1.0)
    plt.plot(time_steps[0:nsteps:step], timings_transport[0:nsteps:step], label="Transport", color='C2', linewidth=2, alpha=1.0)
    leg = plt.legend(loc='right', bbox_to_anchor=(0.5, 0.0, 0.5, 0.5), prop=prop)
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/computing-costs-nolegend-with-smart-ideal.png')
    #plt.savefig(folder_general + '/computing-costs-nolegend-with-smart.png')
    #plt.savefig(folder_general + '/computing-costs-nolegend.png')

    plt.close()


def plot_on_demand_learning_countings():

    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    with open(folder_smart + '/analysis-smart.json') as read_file:
        data = json.load(read_file)

    status_learnings = data.get('cells_where_learning_was_required_at_step')

    # Collect the number of learnings on each step
    learnings = [len(x) for x in status_learnings]

    plt.xlabel('Time Step', fontproperties=prop)
    plt.ylabel('On-demand learnings (at each step)', fontproperties=prop)
    plt.xlim(left=0, right=nsteps)
    plt.ylim(bottom=0, top=np.max(learnings)+1)
    yint = range(0, int(np.max(learnings)+1))
    plt.yticks(yint)
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(learnings, color='C0', linewidth=2)
    plt.tight_layout()
    plt.savefig(folder_general + '/on-demand-learning-countings.png')
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


    plt.xlabel('Time Step', fontproperties=prop)
    plt.ylabel('Accumulated on-demand learnings', fontproperties=prop)
    plt.xlim(left=0, right=nsteps)
    plt.ylim(bottom=0, top=np.max(accum_learnings)+5)
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], accum_learnings[0:nsteps:step], color='C0', linewidth=2)
    plt.tight_layout()
    plt.savefig(folder_general + '/on-demand-learning-total.png')
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

    plt.xlabel('Time Step', fontproperties=prop)
    plt.ylabel('Computing Cost [μs]', fontproperties=prop)
    #plt.xscale('log')
    plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], timings_search[0:nsteps:step], label="Search Time", color='C0', linewidth=2)
    plt.plot(time_steps[0:nsteps:step], timings_estimate[0:nsteps:step], label="Estimate Time", color='C1', linewidth=2)

    leg = plt.legend(loc='lower right', prop=prop)
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/search-traylor-vs-total-learnings.png')
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
                     - np.array(timing_smart.get('smart_equilibrium_search'))
                     - np.array(timing_smart.get('smart_equilibrium_storage')))

    print("speedup         = ", np.sum(np.array(timing_class.get('equilibrium'))) / np.sum(np.array(timing_smart.get('smart_equilibrium'))))
    print("speedup (ideal) = ", np.sum(np.array(timing_class.get('equilibrium'))) / (np.sum(np.array(timing_smart.get('smart_equilibrium')) \
                                                                                            - np.array(timing_smart.get('smart_equilibrium_search')) \
                                                                                            - np.array(timing_smart.get('smart_equilibrium_storage')))))


    plt.xlim(left=0, right=nsteps)
    plt.xlabel('Time Step', fontproperties=prop)
    plt.ylabel('Speedup (-)', fontproperties=prop)
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], speedup[0:nsteps:step], label="Conventional vs. Smart ", color='C0', linewidth=2)
    plt.plot(time_steps[0:nsteps:step], speedup_ideal[0:nsteps:step], label="Conventional vs. Smart (Ideal search)",
             color='C2', linewidth=2, alpha=1.0)

    leg = plt.legend(loc='upper right', prop=prop)
    for line in leg.get_lines():
       line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/speedups.png')
    #plt.show()
    #plt.savefig(folder_general + '/speedups-nolegend-withupperbound.png')
    #plt.savefig(folder_general + '/speedups-nolegend.png')
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
        data_class_ph = data_class[indx_pH]
        data_smart_ph = data_smart[indx_pH]

        data_class_pyrrhotite, data_class_siderite = data_class[indx_phase_amount_pyrrhotite], data_class[indx_phase_amount_siderite]
        data_smart_pyrrhotite, data_smart_siderite = data_smart[indx_phase_amount_pyrrhotite], data_smart[indx_phase_amount_siderite]

        data_class_Hcation  = data_class[indx_Hcation]
        data_class_CO3anion  = data_class[indx_CO3anion]
        data_class_HSanion  = data_class[indx_HSanion]
        data_class_S2anion = data_class[indx_S2anion]
        data_class_HSO4anion   = data_class[indx_HSO4anion]
        data_class_H2Saq = data_class[indx_H2Saq]

        data_smart_Hcation = data_smart[indx_Hcation]
        data_smart_CO3anion = data_smart[indx_CO3anion]
        data_smart_HSanion = data_smart[indx_HSanion]
        data_smart_S2anion = data_smart[indx_S2anion]
        data_smart_HSO4anion = data_smart[indx_HSO4anion]
        data_smart_H2Saq = data_smart[indx_H2Saq]


        abs_error_ph = np.abs(data_class_ph - data_smart_ph)
        abs_error_pyrrhotite = np.abs(data_class_pyrrhotite - data_smart_pyrrhotite)
        abs_error_siderite = np.abs(data_class_siderite - data_smart_siderite)
        abs_error_CO3anion = np.abs(data_class_CO3anion - data_smart_CO3anion)
        abs_error_HSanion = np.abs(data_class_HSanion - data_smart_HSanion)
        abs_error_S2anion = np.abs(data_class_S2anion - data_smart_S2anion)
        abs_error_HSO4anion = np.abs(data_class_HSO4anion - data_smart_HSO4anion)
        abs_error_H2Saq = np.abs(data_class_H2Saq - data_smart_H2Saq)
        abs_error_Hcation = np.abs(data_class_Hcation - data_smart_Hcation)

        data_class_siderite[data_class_siderite<=tol] = 1.0
        data_class_pyrrhotite[data_class_pyrrhotite<=tol] = 1.0
        data_class_ph[data_class_ph<=tol] = 1.0
        data_class_HSanion[data_class_HSanion<=tol] = 1.0
        data_class_S2anion[data_class_S2anion<=tol] = 1.0
        data_class_CO3anion[data_class_CO3anion<=tol] = 1.0
        data_class_HSO4anion[data_class_HSO4anion<=tol] = 1.0
        data_class_H2Saq[data_class_H2Saq<=tol] = 1.0
        data_class_Hcation[data_class_Hcation<=tol] = 1.0


        rel_error_ph = np.divide(abs_error_ph, data_class_ph)
        rel_error_pyrrhotite = np.divide(abs_error_pyrrhotite, data_class_pyrrhotite)
        rel_error_siderite = np.divide(abs_error_siderite, data_class_siderite)
        rel_error_CO3anion = np.divide(abs_error_CO3anion, data_class_CO3anion)
        rel_error_HSanion = np.divide(abs_error_HSanion, data_class_HSanion)
        rel_error_S2anion = np.divide(abs_error_S2anion, data_class_S2anion)
        rel_error_HSO4anion = np.divide(abs_error_HSO4anion, data_class_HSO4anion)
        rel_error_H2Saq = np.divide(abs_error_H2Saq, data_class_H2Saq)
        rel_error_Hcation = np.divide(abs_error_Hcation, data_class_Hcation)

        t = i * dt

        plt.axes(xlim=(-1, 75.1))
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('Relative Error [-]', fontproperties=prop)
        plt.title(titlestr(t), fontproperties=prop, fontsize=20)
        plt.plot(xcells, rel_error_pyrrhotite, label='Pyrrhotite', **line_error('C0'))
        plt.plot(xcells, rel_error_siderite, label='Siderite', **line_error('C1'))
        prop.set_size(12)
        plt.legend(loc='center right', prop=prop)
        plt.savefig(folder_general + f'/error-rel-pyrrhotite-siderite-{i}.png')
        plt.tight_layout()
        plt.close()

        plt.axes(xlim=(-1, 75.1))
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('Relative Error [-]', fontproperties=prop)
        plt.yscale('log')
        plt.title(titlestr(t), fontproperties=prop, fontsize=20)
        plt.plot(xcells, rel_error_Hcation, label=r'$\mathrm{H^+}$', **line_error('darkviolet'))[0],
        plt.plot(xcells, rel_error_HSanion, label=r'$\mathrm{HS^{-}}$', **line_error('C0'))[0],
        plt.plot(xcells, rel_error_S2anion, label=r'$\mathrm{S^{2-}}$', **line_error('C1'))[0],
        plt.plot(xcells, rel_error_CO3anion, label=r'$\mathrm{CO_3^{2-}}$',**line_error('C2'))[0],
        plt.plot(xcells, rel_error_HSO4anion, label=r'$\mathrm{HSO_4^{-}}$',**line_error('red'))[0],
        plt.plot(xcells, rel_error_H2Saq, label=r'$\mathrm{H_2S(aq)}$', **line_error('gold'))[0],
        prop.set_size(12)
        plt.legend(loc='upper right', prop=prop)
        plt.savefig(folder_general + f'/error-rel-aqueous-species-{i}.png')
        plt.tight_layout()
        plt.close()

def mass_conservation():

    #indices = [0, 19, 2399]
    indices = [60, 120, 960]
    for i in indices:

        filearray_res = np.loadtxt(folder_smart + f'/res-{i}.txt')
        filearray_state = np.loadtxt(folder_smart + f'/test-{i}.txt', skiprows=1)

        filearray_res = filearray_res.T
        filearray_state = filearray_state.T

        data_C = np.abs(filearray_state[indx_C])
        data_Ca = np.abs(filearray_state[indx_Ca])
        data_Cl = np.abs(filearray_state[indx_Cl])
        data_H = np.abs(filearray_state[indx_H])
        data_K = np.abs(filearray_state[indx_K])
        data_Mg = np.abs(filearray_state[indx_Mg])
        data_Na = np.abs(filearray_state[indx_Na])
        data_O = np.abs(filearray_state[indx_O])
        data_S = np.abs(filearray_state[indx_S])

        # sigma(Ca) = abs(Ca) if abs(Ca) < eps_b else 1.0
        eps_b = 1e-16
        data_C[data_C<=eps_b] = 1.0
        data_Ca[data_Ca<=eps_b] = 1.0
        data_Cl[data_Cl<=eps_b] = 1.0
        data_H[data_H<=eps_b] = 1.0
        data_K[data_K<=eps_b] = 1.0
        data_Mg[data_Mg<=eps_b] = 1.0
        data_Na[data_Na<=eps_b] = 1.0
        data_O[data_O<=eps_b] = 1.0
        data_S[data_S<=eps_b] = 1.0

        res_shift = indx_C
        data_res_C = np.abs(filearray_res[indx_C-res_shift])
        data_res_Ca = np.abs(filearray_res[indx_Ca-res_shift])
        data_res_Cl = np.abs(filearray_res[indx_Cl-res_shift])
        data_res_K = np.abs(filearray_res[indx_K-res_shift])
        data_res_H = np.abs(filearray_res[indx_H-res_shift])
        data_res_Mg = np.abs(filearray_res[indx_Mg-res_shift])
        data_res_Na = np.abs(filearray_res[indx_Na-res_shift])
        data_res_O = np.abs(filearray_res[indx_O-res_shift])
        data_res_S = np.abs(filearray_res[indx_S-res_shift])

        # error[Ca] = abs(r[Ca])/sigma(Ca)
        error_C = np.divide(data_res_C, data_C)
        error_Ca = np.divide(data_res_Ca, data_Ca)
        error_Cl = np.divide(data_res_Cl, data_Cl)
        error_K = np.divide(data_res_K, data_K)
        error_H = np.divide(data_res_H, data_H)
        error_Mg = np.divide(data_res_Mg, data_Mg)
        error_Na = np.divide(data_res_Na, data_Na)
        error_O = np.divide(data_res_O, data_O)
        error_S = np.divide(data_res_S, data_S)

        #plt.axes(xlim=(-0.01, 0.501), ylim=(1e-16, 5e-7))
        plt.axes(xlim=(-1, 75.1))
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
        plt.plot(xcells, 100 * error_S, label='S', **line_error('C7'))
        plt.plot(xcells, 100 * error_K, label='K', **line_error('C8'))

        plt.legend(loc='lower right', prop=prop)
        plt.savefig(folder_general + f'/mass-balance-residual-{i}.png')
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
    files_smart = [file for file in natsorted( os.listdir(folder_smart) ) if (("res" not in file) and ("analysis" not in file))]
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


    tolerance = 1e-2
    calculate_error(tolerance)
    #mass_conservation()

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

    plot_animation_ph()
    plot_animation_pyrrhotite_siderite()
    plot_animation_aqueous_species()
