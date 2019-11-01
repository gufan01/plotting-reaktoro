import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from joblib import Parallel, delayed
import os, sys
from natsort import natsorted
from matplotlib import animation
import json

# Animation options
animation_starts_at_frame    = 0      # the first frame index to be considered
#animation_ends_at_frame      = 10000  # the last frame index to be considered
#animation_num_frames_to_jump = 24     # the number of frames to jump between current and next
#result_tag = '-long-'
animation_ends_at_frame      = 300  # the last frame index to be considered
animation_num_frames_to_jump = 1     # the number of frames to jump between current and next
result_tag = '-short-'
animation_fps                = 30     # the number of frames per second
animation_interval_wait      = 200    # the time (in milliseconds) to wait between each frame

# Auxiliary animation options
animation_frame_range = range(animation_starts_at_frame, animation_ends_at_frame, animation_num_frames_to_jump)

# Auxiliary time related constants
second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
year = 365 * day

# Discretisation parameters
xl = 0.0          # the x-coordinate of the left boundary
xr = 1.0          # the x-coordinate of the right boundary
nsteps = 10000    # the number of steps in the reactive transport simulation
ncells = 100      # the number of cells in the discretization
eqreltol = 1e-1   # relative tolerance in equilibrium
eqabstol = 1e-12  # absolute tolerance in equilibrium
kinreltol = 1e-1  # relative tolerance in kinetics
kinabstol = 1e-5  # absolute tolerance in kinetics
cutoff = -1e-5    # absolute tolerance in kinetics

D  = 1.0e-9       # the diffusion coefficient (in units of m2/s)
v  = 1.0/day      # the fluid pore velocity (in units of m/s)
dt = 30 * minute  # the time step (in units of s)
T = 60.0          # the temperature (in units of K)
P = 100           # the pressure (in units of Pa)
phi = 0.1           # the porosity

tag_smart = "-dt-" + "{:d}".format(dt) + \
      "-ncells-" + str(ncells) + \
      "-nsteps-" + str(nsteps) + \
      "-eqrel-" + "{:.{}e}".format(eqreltol, 1) + \
      "-eqabs-" + "{:.{}e}".format(eqabstol, 1) + \
      "-kinrel-" + "{:.{}e}".format(kinreltol, 1) + \
      "-kinabs-" + "{:.{}e}".format(kinabstol, 1) + \
      "-cutoff-" + "{:.{}e}".format(cutoff, 1)

tag_class = "-dt-" + "{:d}".format(dt) + \
      "-ncells-" + str(ncells) + \
      "-nsteps-" + str(nsteps)

test_tag_smart = tag_smart + "-smart-kin-conv-eq"
test_tag_class = tag_class + "-conv-kin-conv-eq"

folder = 'rt-Jf'
folder_class = 'rt-exact-hessian'
#folder = 'rt-exact-hessian'
#folder_class = 'rt-exact-hessian'
folder_smart   = folder + test_tag_smart
folder_class   = folder_class + test_tag_class
folder_general = "results-smart-kinetics-videos" + result_tag + folder + tag_smart

os.system('mkdir -p ' + folder_general)

# Indices of the loaded data to plot
indx_ph        = 0
indx_Hcation   = 1
indx_Cacation  = 2
indx_Mgcation  = 3
indx_HCO3anion = 4
indx_CO2aq     = 5
indx_calcite   = 6
indx_dolomite  = 7

xcells = np.linspace(xl, xr, ncells)  # the x-coordinates of the plots

mpl.rcParams['font.sans-serif'] ='Century Gothic'

def line_empty_marker(color):
    return {'markerfacecolor': 'white', 'markeredgecolor':color, 'markersize': 6, 'markeredgewidth': 1.5 }

def line_filled_marker(color):
    return {'color': color, 'markersize': 6, 'markeredgewidth': 1.5 }

def line(color):
    return {'linestyle': '-', 'color': color, 'zorder': 1, 'linewidth': 2}

def titlestr(t):
    d = int(t / day)                 # The number of days
    h = int(int(t % day) / hour)     # The number of remaining hours
    m = int(int(t % hour) / minute)  # The number of remaining minutes
    return '{:>3d}d {:>2}h {:>2}m'.format(int(d), str(int(h)).zfill(2), str(int(m)).zfill(2))


def plot_animation_ph():

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-0.01, 0.501), ylim=(2.5, 12.0))
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
        print("On pH animation index: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_ph = data_class[indx_ph]
        data_smart_ph = data_smart[indx_ph]
        objects[0].set_data(xcells, data_class_ph)
        objects[1].set_data(xcells[status[i]==0], data_smart_ph[status[i]==0])
        objects[2].set_data(xcells[status[i]==1], data_smart_ph[status[i]==1])
        ax.set_title(titlestr(t))
        return tuple(objects)


    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)

    anim.save(folder_general + '/pH.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])


def plot_animation_calcite_dolomite():

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-0.01, 0.501), ylim=(-0.1, 2.1))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Mineral Volume [%$_{\mathsf{vol}}$]')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label='Calcite', **line('C0'))[0],
        ax.plot([], [], label='Dolomite', **line('C1'))[0],
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
        print("On calcite-dolomite animation index: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_calcite, data_class_dolomite = data_class[indx_calcite], data_class[indx_dolomite]
        data_smart_calcite, data_smart_dolomite = data_smart[indx_calcite], data_smart[indx_dolomite]
        objects[0].set_data(xcells, data_class_calcite * 100/(1 - phi))
        objects[1].set_data(xcells, data_class_dolomite * 100/(1 - phi))
        objects[2].set_data(xcells[status[i]==0], data_smart_calcite[status[i]==0] * 100/(1 - phi))
        objects[3].set_data(xcells[status[i]==1], data_smart_calcite[status[i]==1] * 100/(1 - phi))
        objects[4].set_data(xcells[status[i]==0], data_smart_dolomite[status[i]==0] * 100/(1 - phi))
        objects[5].set_data(xcells[status[i]==1], data_smart_dolomite[status[i]==1] * 100/(1 - phi))
        ax.set_title(titlestr(t))
        return tuple(objects)


    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)

    anim.save(folder_general + '/calcite-dolomite.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])


def plot_animation_aqueous_species():

    # Plot of mineral's volume the space coordinates
    fig = plt.figure()
    ax = plt.axes(xlim=(-0.01, 0.501), ylim=(0.5e-5, 2))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Concentration [molal]')
    ax.set_yscale('log')
    ax.set_title(titlestr(0.0))
    objects = [
        ax.plot([], [], label=r'$\mathrm{Ca^{2+}}$', **line('C0'))[0],
        ax.plot([], [], label=r'$\mathrm{Mg^{2+}}$', **line('C1'))[0],
        ax.plot([], [], label=r'$\mathrm{HCO_3^{-}}$',**line('C2'))[0],
        ax.plot([], [], label=r'$\mathrm{CO_2(aq)}$',**line('red'))[0],
        ax.plot([], [], label=r'$\mathrm{H^+}$', **line('darkviolet'))[0],
        ax.plot([], [], 'o', **line_empty_marker('C0'))[0],
        ax.plot([], [], 'o', **line_filled_marker('C0'))[0],
        ax.plot([], [], 'o', **line_empty_marker('C1'))[0],
        ax.plot([], [], 'o', **line_filled_marker('C1'))[0],
        ax.plot([], [], 'o', **line_empty_marker('C2'))[0],
        ax.plot([], [], 'o', **line_filled_marker('C2'))[0],
        ax.plot([], [], 'o', **line_empty_marker('red'))[0],
        ax.plot([], [], 'o', **line_filled_marker('red'))[0],
        ax.plot([], [], 'o', **line_empty_marker('darkviolet'))[0],
        ax.plot([], [], 'o', **line_filled_marker('darkviolet'))[0],
    ]
    ax.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
    ax.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
    ax.legend(loc='upper right')

    def init():
        return tuple(objects)

    def animate(i):
        print("On aqueous-species animation index: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T

        data_class_cacation  = data_class[indx_Cacation]
        data_class_mgcation  = data_class[indx_Mgcation]
        data_class_hco3anion = data_class[indx_HCO3anion]
        data_class_co2aq     = data_class[indx_CO2aq]
        data_class_hcation   = data_class[indx_Hcation]

        data_smart_cacation  = data_smart[indx_Cacation]
        data_smart_mgcation  = data_smart[indx_Mgcation]
        data_smart_hco3anion = data_smart[indx_HCO3anion]
        data_smart_co2aq     = data_smart[indx_CO2aq]
        data_smart_hcation   = data_smart[indx_Hcation]

        objects[0].set_data(xcells, data_class_cacation)
        objects[1].set_data(xcells, data_class_mgcation)
        objects[2].set_data(xcells, data_class_hco3anion)
        objects[3].set_data(xcells, data_class_co2aq)
        objects[4].set_data(xcells, data_class_hcation)
        objects[5].set_data(xcells[status[i]==0], data_smart_cacation[status[i]==0])
        objects[6].set_data(xcells[status[i]==1], data_smart_cacation[status[i]==1])
        objects[7].set_data(xcells[status[i]==0], data_smart_mgcation[status[i]==0])
        objects[8].set_data(xcells[status[i]==1], data_smart_mgcation[status[i]==1])
        objects[9].set_data(xcells[status[i]==0], data_smart_hco3anion[status[i]==0])
        objects[10].set_data(xcells[status[i]==1], data_smart_hco3anion[status[i]==1])
        objects[11].set_data(xcells[status[i]==0], data_smart_co2aq[status[i]==0])
        objects[12].set_data(xcells[status[i]==1], data_smart_co2aq[status[i]==1])
        objects[13].set_data(xcells[status[i]==0], data_smart_hcation[status[i]==0])
        objects[14].set_data(xcells[status[i]==1], data_smart_hcation[status[i]==1])
        ax.set_title(titlestr(t))
        return tuple(objects)


    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)

    anim.save(folder_general + '/aqueous-species.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])


if __name__ == '__main__':

    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    with open(folder_smart + '/analysis-smart-kin-conventional-eq.json') as read_file:
        status_data = json.load(read_file)

    status_learnings = status_data.get('smart_kinetics_cells_where_learning_was_required_at_step')
    status = np.ones([nsteps, ncells])
    for i in range(0, nsteps):
        for j in range(0, len(status_learnings[i])):
            status[i][j] = 0

    print("Collecting files...")
    # Collect files with results corresponding to smart or reference (classical) solver
    files_smart = [file for file in natsorted(os.listdir(folder_smart)) if ("analysis" not in file)]
    files_class = [file for file in natsorted(os.listdir(folder_class)) if ("analysis" not in file)]

    plot_animation_ph()
    plot_animation_calcite_dolomite()
    plot_animation_aqueous_species()

    print("Finished plotting animations!")
