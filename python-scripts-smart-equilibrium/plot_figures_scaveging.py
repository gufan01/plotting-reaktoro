import numpy as np
import matplotlib.pyplot as plt
import os
from natsort import natsorted
import matplotlib as mpl
import json


# Options for the figure plotting
plot_at_selected_steps = [1, 10, 60, 120, 240, 480, 960, 1200, 2400, 3600, 4800, 6400, 7200]  # the time steps at which the results are plotted
#plot_at_selected_steps = [1, 10, 60, 120, 240, 480, 960, 1000]  # the time steps at which the results are plotted
#plot_at_selected_steps = [1, 10, 20, 40, 80, 100]  # the time steps at which the results are plotted

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
folder_smart = 'cpp-reactivetransport-old-demo/results-debey-huckel-dt-3594.24-ncells-100-nsteps-10000-eqreltol-1.0e-01-eqabstol-1.0e-08-smart'
folder_class = 'cpp-reactivetransport-old-demo/results-debey-huckel-dt-3594.24-ncells-100-nsteps-10000-reference'
folder_general = "plots-results-debey-huckel-dt-3594.24-ncells-100-nsteps-10000"
#folder_smart = 'cpp-reactivetransport-old-demo/results-pitzer-dt-3594.24-ncells-100-nsteps-10000-eqreltol-1.0e-01-eqabstol-1.0e-08-smart'
#folder_class = 'cpp-reactivetransport-old-demo/results-pitzer-dt-3594.24-ncells-100-nsteps-10000-reference'
#folder_general = "plots-results-pitzer-dt-3594.24-ncells-100-nsteps-10000"
os.system('mkdir -p ' + folder_general)

fillz = len(str(123))

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


def plot_figures_ph():

    for i in plot_at_selected_steps:
        print("On pH figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_ph = data_class[pH]
        data_smart_ph = data_smart[pH]
        #plt.axes(xlim=(-0.01,1.001))
        plt.axes(xlim=(-0.01,1.001), ylim=(4, 9))
        plt.xlabel('Distance [m]')
        plt.ylabel('pH')
        plt.title(titlestr(t))
        plt.plot(xcells, data_class_ph, label='pH', **line('teal'))
        plt.plot(xcells[status[i-1]==0], data_smart_ph[status[i-1]==0], 'o', **line_empty_marker('teal'))
        plt.plot(xcells[status[i-1]==1], data_smart_ph[status[i-1]==1], 'o', **line_filled_marker('teal'))
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='lower right')
        plt.savefig(folder_general + '/pH-{}.png'.format(i))
        plt.tight_layout()
        plt.close()


def plot_figures_pyrrhotite_siderite():

    for i in plot_at_selected_steps:
        print("On pyrrhotite-siderite figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T

        data_class_pyrrhotite, data_class_siderite = data_class[pyrrhotite], data_class[siderite]
        data_smart_pyrrhotite, data_smart_siderite = data_smart[pyrrhotite], data_smart[siderite]
        plt.axes(xlim=(-0.01,1.001), ylim=(-0.1, 6.1))
        #plt.axes(xlim=(-0.01,1.001), ylim=(-0.1, 2.1))
        plt.ylabel('Mineral Volume')

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
        plt.savefig(folder_general + '/pyrrhotite-siderite-{}.png'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_aqueous_species():

    for i in plot_at_selected_steps:
        print("On aqueous-species figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T

        data_class_Hcation  = data_class[Hcation]
        data_class_HSanion  = data_class[HSanion]
        data_class_S2anion = data_class[S2anion]
        data_class_SO4anion     = data_class[SO4anion]
        data_class_HSO4anion   = data_class[HSO4anion]
        data_class_H2Saq = data_class[H2Saq]

        data_smart_Hcation = data_smart[Hcation]
        data_smart_HSanion = data_smart[HSanion]
        data_smart_S2anion = data_smart[S2anion]
        data_smart_SO4anion = data_smart[SO4anion]
        data_smart_HSO4anion = data_smart[HSO4anion]
        data_smart_H2Saq = data_smart[H2Saq]

        #plt.axes(xlim=(-0.01,1.001))
        plt.axes(xlim=(-0.01,1.001), ylim=(1e-12, 1e1))
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
        plt.savefig(folder_general + '/aqueous-species-{}.png'.format(i))
        plt.tight_layout()
        plt.close()


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

    plt.xlabel('Time Step')
    plt.xlim(left=0, right=nsteps)
    plt.ylim(bottom=0, top=np.max(learnings)+1)
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


    plt.xlabel('Time Step')
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
                     - np.array(timing_smart.get('smart_equilibrium_nearest_neighbor_search')))

    print("speedup         = ", np.average(speedup))
    print("speedup (ideal) = ", np.average(speedup_ideal))

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
    plt.savefig(folder_general + '/speedups.png')
    #plt.show()
    #plt.savefig(folder_general + '/speedups-nolegend-withupperbound.png')
    #plt.savefig(folder_general + '/speedups-nolegend.png')
    plt.close()

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
    files_smart = [file for file in natsorted( os.listdir(folder_smart) ) if ("analysis" not in file)]
    files_class = [file for file in natsorted( os.listdir(folder_class) ) if ("analysis" not in file)]

    plot_on_demand_learning_countings()
    plot_on_demand_learnings_total()
    plot_speedups()
    plot_computing_costs_vs_total_learnings()
    plot_computing_costs()
    plot_figures_ph()
    plot_figures_aqueous_species()
    plot_figures_pyrrhotite_siderite()
    #plot_figures_pyrrhotite_siderite_moles()
