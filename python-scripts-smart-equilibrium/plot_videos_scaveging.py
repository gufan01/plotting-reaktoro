import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, sys
from natsort import natsorted
from matplotlib import animation
import json

# Animation options
animation_starts_at_frame    = 0      # the first frame index to be considered
#animation_ends_at_frame      = 10000  # the last frame index to be considered
#animation_num_frames_to_jump = 24     # the number of frames to jump between current and next
animation_ends_at_frame      = 10 * 30  # the last frame index to be considered
animation_num_frames_to_jump = 1     # the number of frames to jump between current and next
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
nsteps = 10000     # the number of steps in the reactive transport simulation
ncells = 100      # the number of cells in the discretization
reltol = 1e-1     # relative tolerance
abstol = 1e-8    # absolute tolerance

D  = 1.0e-9       # the diffusion coefficient (in units of m2/s)
v  = 1.0/day     # the fluid pore velocity (in units of m/s)
dt = 30 * minute   # the time step (in units of s)
T = 60.0          # the temperature (in units of K)
P = 100           # the pressure (in units of Pa)
phi = 0.1         # the porosity

tag = "-dt-" + "{:d}".format(dt) + \
      "-ncells-" + str(ncells) + \
      "-nsteps-" + str(nsteps) + \
      "-reltol-" + "{:.{}e}".format(reltol, 1) + \
      "-abstol-" + "{:.{}e}".format(abstol, 1)

test_tag_smart = tag + "-smart"
test_tag_class = tag + "-reference"

# /results-addAqueousPhase

#folder = 'results-deltan-with-normalization'
#folder_smart   = folder + test_tag_smart
#folder_class   = folder + test_tag_class
#folder_smart = 'cpp-reactivetransport-old-demo/results-debey-huckel-dt-3594.24-ncells-100-nsteps-10000-eqreltol-1.0e-01-eqabstol-1.0e-08-smart'
#folder_class = 'cpp-reactivetransport-old-demo/results-debey-huckel-dt-3594.24-ncells-100-nsteps-10000-reference'
#folder_general = "plots-results-debey-huckel-dt-3594.24-ncells-100-nsteps-10000"
folder_smart = 'cpp-reactivetransport-old-demo/results-pitzer-dt-3594.24-ncells-100-nsteps-10000-eqreltol-1.0e-01-eqabstol-1.0e-08-smart'
folder_class = 'cpp-reactivetransport-old-demo/results-pitzer-dt-3594.24-ncells-100-nsteps-10000-reference'
folder_general = "plots-results-pitzer-dt-3594.24-ncells-100-nsteps-10000"
os.system('mkdir -p ' + folder_general)

# Indices of the loaded data to plot
pH = 0
Hcation = 1
HSanion = 2
S2anion = 3
SO4anion = 4
HSO4anion = 5
H2Saq = 6
pyrrhotite = 7
siderite = 8

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
    ax = plt.axes(xlim=(-0.01, 0.501))
    plt.axes(xlim=(-0.01, 0.501), ylim=(4, 9))
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
        data_class_ph = data_class[pH]
        data_smart_ph = data_smart[pH]
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
    ax = plt.axes(xlim=(-0.01, 0.501))
    plt.axes(xlim=(-0.01, 0.501), ylim=(-0.1, 6.1))
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Mineral Volume')
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
        print("On calcite-dolomite animation index: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_calcite, data_class_dolomite = data_class[pyrrhotite], data_class[pyrrhotite]
        data_smart_calcite, data_smart_dolomite = data_smart[siderite], data_smart[siderite]
        objects[0].set_data(xcells, data_class_calcite)
        objects[1].set_data(xcells, data_class_dolomite)
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
    #ax = plt.axes(xlim=(-0.01, 0.501), ylim=(0.5e-5, 2))
    ax = plt.axes(xlim=(-0.01, 0.501), ylim=(1e-12, 1e1))
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
        print("On aqueous-species animation index: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T

        data_class_Hcation = data_class[Hcation]
        data_class_HSanion = data_class[HSanion]
        data_class_S2anion = data_class[S2anion]
        data_class_SO4anion = data_class[SO4anion]
        data_class_HSO4anion = data_class[HSO4anion]
        data_class_H2Saq = data_class[H2Saq]

        data_smart_Hcation = data_smart[Hcation]
        data_smart_HSanion = data_smart[HSanion]
        data_smart_S2anion = data_smart[S2anion]
        data_smart_SO4anion = data_smart[SO4anion]
        data_smart_HSO4anion = data_smart[HSO4anion]
        data_smart_H2Saq = data_smart[H2Saq]

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


    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=animation_frame_range, interval=animation_interval_wait, blit=True)

    anim.save(folder_general + '/aqueous-species.mp4', fps=animation_fps, extra_args=['-vcodec', 'libx264'])

def count_trainings(status):
    counter = [st for st in status if st == 0]
    return (len(counter), len(status))

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
    files_smart = [file for file in natsorted(os.listdir(folder_smart)) if ("analysis" not in file)]
    files_class = [file for file in natsorted(os.listdir(folder_class)) if ("analysis" not in file)]

    plot_animation_ph()
    plot_animation_calcite_dolomite()
    plot_animation_aqueous_species()

    print("Finished plotting animations!")
