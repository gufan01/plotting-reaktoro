import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

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

# /results-addAqueousPhase

folder = 'results-Pitzer/'
folder_smart   = "../" + folder + "results" + test_tag_smart
folder_class   = "../" + folder + "results" + test_tag_class
folder_general = "results-pitzer" + tag

fillz = len(str(123))
# Output properties 
output_quantities = """
    pH
    speciesMolality(H+)
    speciesMolality(Ca++)
    speciesMolality(Mg++)
    speciesMolality(HCO3-)
    speciesMolality(CO2(aq))
    phaseVolume(Calcite)
    phaseVolume(Dolomite)
""".split()

# Indices of the loaded data to plot
indx_ph        = 0
indx_Hcation   = 1
indx_Cacation  = 2
indx_Mgcation  = 3
indx_HCO3anion = 4
indx_CO2aq     = 5
indx_calcite   = 6
indx_dolomite  = 7

# Plotting params
circ_area = 20.0
zoom = 0.5
custom_font = {}


#font
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] ='Fira Sans'
mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 'medium'
#tick label spacing and tick width
mpl.rcParams['xtick.major.pad'] = 4
mpl.rcParams['ytick.major.pad'] = 5
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
#legend style
mpl.rcParams['legend.frameon'] = True
mpl.rcParams['legend.numpoints'] = 3
#mpl.rcParams['backend'] = 'PDF'
#mpl.rcParams['savefig.dpi'] = 200

def myfig(width, hight):
    plt.clf()
    fig = plt.figure(figsize=(width, hight))
    ax = fig.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    return fig, ax

def empty_marker(color):
    return {'facecolor': 'white', 'edgecolor': color, 's': circ_area, 'zorder': 2}

def filled_marker(color):
    return {'color': color, 's': circ_area, 'zorder': 2}

def line_empty_marker(color):
    return {'markerfacecolor': 'white', 'markeredgecolor':color, 'markersize': 5.0}

def line_filled_marker(color):
    return {'color': color, 'markersize': 5.0}

def line(color):
    return {'linestyle': '-', 'color': color, 'zorder': 1}

def titlestr(t):
    t = t / minute   # Convert from seconds to minutes
    h = int(t) / 60  # The number of hours
    m = int(t) % 60  # The number of remaining minutes
    return 'Time: %2dh %2dm' % (h, m)

def make_results_folders():
    os.system('mkdir -p ' + folder_general + '/figures/ph')
    os.system('mkdir -p ' + folder_general + '/figures/aqueous-species')
    os.system('mkdir -p ' + folder_general + '/figures/calcite-dolomite')
    os.system('mkdir -p ' + folder_general + '/videos')
    os.system('mkdir -p ' + folder_general + '/results')

def get_step(file):
    return int(file.split('.')[0].split('-')[1])

def compare_chemistry(file, status, selected_steps):

    # Get the number of the step
    step = get_step(file)

    inf = 1e16
    dstep = 5

    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] ='Fira Sans'
    mpl.rcParams['font.size'] = 12
    mpl.rcParams['legend.fontsize'] = 'medium'
    #tick label spacing and tick width
    mpl.rcParams['xtick.major.pad'] = 4
    mpl.rcParams['ytick.major.pad'] = 5
    mpl.rcParams['xtick.major.width'] = 1
    mpl.rcParams['ytick.major.width'] = 1
    #legend style
    mpl.rcParams['legend.frameon'] = True
    mpl.rcParams['legend.numpoints'] = 3

    #if step == 0 or step % 5 == 0:
    if step >= 0:

        # Log the printing of the step
        # print('Plotting figure', step, '...')

        # The current time of the data loaded
        t = (step + 1) * dt

        # Load the data from the filearray skipping the 1st row
        filearray_smart = np.loadtxt(folder_smart + '/' + file, skiprows=1)
        data_smart = filearray_smart.T
        filearray_class = np.loadtxt(folder_class + '/' + file, skiprows=1)
        data_class = filearray_class.T

        # Number of digits in the total number of the steps
        ndigits = len(str(nsteps))

        # Cells coordinates
        cells = np.linspace(xl, xr, ncells)

        shift = 0.01 * cells[-1]


        # Plot change of pH wrt the space coordinates
        fig, ax = myfig(6, 3)
        if zoom == 0.5: ax.set_xlim(left=cells[0]-shift, right=cells[int(ncells/2)]+shift) # 50%
        elif zoom == 1:  ax.set_xlim(left=cells[0]-shift, right=cells[-1]+shift) # 100 %
        ax.set_ylim(bottom=2.5, top=12.0)
        ax.set_title(titlestr(t), **custom_font)
        ax.set_xlabel('Distance [m]', **custom_font)
        ax.set_ylabel('pH', **custom_font)
        ph = data_class[indx_ph]
        ax.plot(cells, ph, label='pH', **line('teal'))
        ph = data_smart[indx_ph]
        # status 0 - training
        # status 1 - prediction
        ax.scatter(cells[status[step]==0], ph[status[step]==0], label='Learning / Smart', **empty_marker('teal'))
        ax.scatter(cells[status[step]==1], ph[status[step]==1], label='Prediction / Smart', **filled_marker('teal'))
        #ax.plot(float(np.inf), float(np.inf), 'o', label='Smart Prediction', **line_filled_marker('black'))
        #ax.plot(float(np.inf), float(np.inf), 'o', **line_empty_marker('black'))
        ax.legend(loc='lower right', prop={'size': 11})
        plt.tight_layout()
        ax.grid(color='lightgray', linestyle=':', linewidth=1)
        #plt.show()
        fig.savefig(folder_general + '/figures/ph/%s.png' % (str(step+1).zfill(fillz)))
        if (step + 1) in selected_steps: fig.savefig(folder_general + '/figures/ph/%s.pdf' % (str(step+1).zfill(fillz)))

        # Plot of mineral's volume the space coordinates
        fig, ax = myfig(6, 3)
        if zoom == 0.5: ax.set_xlim(left=cells[0]-shift, right=cells[int(ncells/2)]+shift) # 50%
        elif zoom == 1:  ax.set_xlim(left=cells[0]-shift, right=cells[-1]+shift) # 100 %
        ax.set_ylim(bottom=-0.25, top=2.1+shift)
        ax.set_title(titlestr(t), **custom_font)
        ax.set_xlabel('Distance [m]', **custom_font)
        ax.set_ylabel('Mineral Volume [%$_{\mathsf{vol}}$]', **custom_font)

        data_calcite, data_dolomite = data_class[indx_calcite], data_class[indx_dolomite]
        ax.plot(cells, data_calcite * 100, label='Calcite', **line('indianred'))
        ax.plot(cells, data_dolomite * 100, label='Dolomite', **line('royalblue'))

        data_calcite, data_dolomite  = data_smart[indx_calcite], data_smart[indx_dolomite]
        ax.scatter(cells[status[step]==0], data_calcite[status[step]==0] * 100, **empty_marker('indianred'))
        ax.scatter(cells[status[step]==1], data_calcite[status[step]==1] * 100, **filled_marker('indianred'))
        ax.scatter(cells[status[step]==0], data_dolomite[status[step]==0] * 100, **empty_marker('royalblue'))
        ax.scatter(cells[status[step]==1], data_dolomite[status[step]==1] * 100, **filled_marker('royalblue'))
        ax.plot(float(inf), float(inf), 'o', label='Prediction / Smart', **line_filled_marker('black'))
        ax.plot(float(inf), float(inf), 'o', label='Learning / Smart', **line_empty_marker('black'))

        ax.legend(loc='center right', prop={'size': 11}, numpoints=1)
        ax.grid(color='lightgray', linestyle=':', linewidth=1)
        plt.tight_layout()

        fig.savefig(folder_general + '/figures/calcite-dolomite/%s.png' % (str(step+1).zfill(fillz)))
        if (step + 1) in selected_steps: fig.savefig(folder_general + '/figures/calcite-dolomite/%s.pdf' % (str(step+1).zfill(fillz)))

        # Plot of aqueous species's concentration the space coordinates
        fig, ax = myfig(6, 3.5)
        if zoom == 0.5: ax.set_xlim(left=cells[0]-shift, right=cells[int(ncells/2)]+shift) # 50%
        elif zoom == 1: ax.set_xlim(left=cells[0]-shift, right=cells[-1]+shift) # 100 %
        ax.set_ylim(bottom=0.5e-5, top=2)
        ax.set_yscale('log')
        ax.set_title(titlestr(t))
        ax.set_xlabel('Distance [m]', **custom_font)
        ax.set_ylabel('Concentration [molal]', **custom_font)

        data_cacation  = data_class[indx_Cacation]
        data_mgcation  = data_class[indx_Mgcation]
        data_hco3anion = data_class[indx_HCO3anion]
        data_co2aq     = data_class[indx_CO2aq]
        data_hcation   = data_class[indx_Hcation]
        ax.plot(cells, data_cacation, label=r'$\mathrm{Ca^{2+}}$', **line('steelblue'))
        ax.plot(cells, data_mgcation, label=r'$\mathrm{Mg^{2+}}$', **line('darkorange'))
        ax.plot(cells, data_hco3anion, label=r'$\mathrm{HCO_3^{-}}$',**line('forestgreen'))
        ax.plot(cells, data_co2aq, label=r'$\mathrm{CO_2(aq)}$',**line('red'))
        ax.plot(cells, data_hcation, label=r'$\mathrm{H^+}$', **line('darkviolet'))

        data_cacation  = data_smart[indx_Cacation]
        data_mgcation  = data_smart[indx_Mgcation]
        data_hco3anion = data_smart[indx_HCO3anion]
        data_co2aq     = data_smart[indx_CO2aq]
        data_hcation   = data_smart[indx_Hcation]

        ax.plot(float(inf), float(inf), 'o', label='Prediction / Smart', **line_filled_marker('black'))
        ax.plot(float(inf), float(inf), 'o', label='Learning / Smart', **line_empty_marker('black'))

        ax.scatter(cells[status[step]==0], data_cacation[status[step]==0], **empty_marker('steelblue'))
        ax.scatter(cells[status[step]==1], data_cacation[status[step]==1], **filled_marker('steelblue'))

        ax.scatter(cells[status[step]==0], data_mgcation[status[step]==0], **empty_marker('darkorange'))
        ax.scatter(cells[status[step]==1], data_mgcation[status[step]==1], **filled_marker('darkorange'))

        ax.scatter(cells[status[step]==0], data_hco3anion[status[step]==0], **empty_marker('forestgreen'))
        ax.scatter(cells[status[step]==1], data_hco3anion[status[step]==1], **filled_marker('forestgreen'))

        ax.scatter(cells[status[step]==0], data_co2aq[status[step]==0], **empty_marker('red'))
        ax.scatter(cells[status[step]==1], data_co2aq[status[step]==1], **filled_marker('red'))

        ax.scatter(cells[status[step]==0], data_hcation[status[step]==0], **empty_marker('darkviolet'))
        ax.scatter(cells[status[step]==1], data_hcation[status[step]==1], **filled_marker('darkviolet'))

        ax.legend(loc='upper right', prop={'size': 11}, numpoints=1)
        ax.grid(color='lightgray', linestyle=':', linewidth=1)
        plt.tight_layout()

        fig.savefig(folder_general + '/figures/aqueous-species/%s.png' % (str(step+1).zfill(fillz)))
        if (step + 1) in selected_steps: fig.savefig(folder_general + '/figures/aqueous-species/%s.pdf' % (str(step+1).zfill(fillz)))

        plt.close('all')

def plot_cpu_times(files, title):

    data = [[], []]

    conv  = 0
    smart = 1
    time_transport = 0
    time_equilibrium = 1
    time_estimate = 2
    time_search = 3
    time_mat_vect = 4
    time_acceptance = 5
    time_learn = 6
    time_store = 7
    time_gibbs_min = 8
    tree_size = 9

    time_conv_learn = 2

    indx = 1e16

    # Load the data from the filearray skipping the 1st row
    for file in files:
        file_name = file[0]
        is_smart  = file[1]

        if is_smart:
            indx = smart
        else:
            indx = conv
        folder = folder_smart if is_smart else folder_class
        data_ = np.loadtxt(folder + '/' + file_name)
        data[indx] = data_.T

    time = np.linspace(0, dt * nsteps / minute, nsteps)

    fig, ax = myfig(6, 3)
    ax.set_xlim(left=0, right=time[-1])
    #ax.set_ylim(bottom=0.0, top=0.8)
    ax.set_xlabel('Time (minute)', **custom_font)
    ax.set_ylabel('CPU time, step-wise (sec)', **custom_font)
    ax.plot(time, data[conv][time_equilibrium], label="Chemical equilibrium  (Conventional)",
             color='steelblue')
    ax.plot(time, data[smart][time_equilibrium], label="Chemical equilibrium  (Smart)",
             color='darkorange')
    ax.plot(time, data[smart][time_transport], label="Reactive transport",
             color='forestgreen')
    ax.grid(color='lightgray', linestyle=':', linewidth=1)
    ax.legend(loc='lower right')
    ax.set_yscale('log')
    plt.tight_layout()
    fig.savefig(folder_general + '/results/cpu-time-cw' + tag + '.pdf')

    accum_rt_cost = np.zeros(nsteps)
    accum_sm_cost = np.zeros(nsteps)
    accum_conv_cost = np.zeros(nsteps)
    for i in range(nsteps):
        accum_rt_cost[i] = np.sum(data[smart][time_transport][0:i])
        accum_sm_cost[i] = np.sum(data[smart][time_equilibrium][0:i])
        accum_conv_cost[i] = np.sum(data[conv][time_equilibrium][0:i])

    fig, ax = myfig(6, 3)
    ax.set_xlim(left=0, right=time[-1])
    ax.set_ylim(bottom=0.0, top=1.05 * np.max(accum_conv_cost))
    ax.set_xlabel('Time (minute)', **custom_font)
    ax.set_ylabel('CPU time, cumulative (sec)', **custom_font)
    ax.plot(time, accum_conv_cost, label="Chemical equilibrium (Conventional)", color='steelblue')
    ax.plot(time, accum_sm_cost, label="Chemical equilibrium  (Smart)", color='darkorange')
    ax.plot(time, accum_rt_cost, label="Reactive transport", color='forestgreen')
    ax.grid(color='lightgray', linestyle=':', linewidth=1)
    ax.legend(loc='upper left')
    plt.tight_layout()
    fig.savefig(folder_general + '/results/accumulated' + tag + '.pdf')

    speedup       = data[conv][time_equilibrium] / data[smart][time_equilibrium]
    speedup_ideal = data[conv][time_equilibrium] / (data[smart][time_equilibrium] - data[smart][time_search] - data[smart][time_store])

    fig, ax = myfig(6, 3)
    ax.set_xlim(left=0, right=time[-1])
    ax.set_ylim(bottom=0.0, top=1.05 * np.max(speedup_ideal))
    ax.set_xlabel('Time (minute)', **custom_font)
    ax.set_ylabel('Speedup (-)', **custom_font)
    ax.plot(time, speedup_ideal, color='yellowgreen', label="Speedup: Conventional vs. Smart (Ideal search)")
    ax.plot(time, speedup, color='forestgreen', label="Speedup: Conventional vs. Smart ")
    ax.grid(color='lightgray', linestyle=':', linewidth=1)
    ax.legend(loc='upper right')
    plt.tight_layout()
    fig.savefig(folder_general + '/results/speedup' + tag + '.pdf')

    '''
    fig, ax = myfig(6, 3)
    ax.set_xlim(left=0, right=time[-1])
    ax.set_ylim(bottom=0.0, top=0.5*np.max(data[conv][time_equilibrium]))
    ax.set_xlabel('Time (minute)', **custom_font)
    ax.set_ylabel('CPU time, step-wise (sec)', **custom_font)
    ax.plot(time, data[conv][time_equilibrium], label="Chemical equilibrium (Conventional)",
             color='steelblue')
    ax.plot(time, data[smart][time_equilibrium] - data[smart][time_search] - data[smart][time_store],
             label="Chemical equilibrium (Smart, Ideal Search)",
             color='orangered')
    ax.plot(time, data[smart][time_transport], label="Reactive transport",
             color='forestgreen')
    ax.grid(color='lightgray', linestyle=':', linewidth=1)
    ax.legend(loc='upper right')
    ax.set_title(title, **custom_font)
    plt.tight_layout()
    fig.savefig(folder_general + '/results/cpu-time-cw-with-ideal' + tag + '.pdf')
    '''

    fig, ax = myfig(6, 3)
    ax.set_xlim(left=0, right=time[-1])
    #ax.set_ylim(bottom=0, top=0.015)
    ax.set_xlabel('Time (minute)', **custom_font)
    ax.set_ylabel('CPU time, step-wise (sec)', **custom_font)
    ax.plot(time, data[smart][time_equilibrium] - data[smart][time_learn], label=r"Smart, $O(K)$",
             color='darkorange')
    ax.plot(time, data[smart][time_equilibrium] - data[smart][time_search] - data[smart][time_store] - data[smart][time_learn],
             label=r"Smart, Ideal Search",
             color='orangered')
    ax.grid(color='lightgray', linestyle=':', linewidth=1)
    ax.legend(loc='upper left')
    plt.tight_layout()
    fig.savefig(folder_general + '/results/smart-cw' + tag + '.pdf')

    fig, ax = myfig(6, 3)
    #ax.set_xlim(left=0, right=data[indx_eq_search_smart][2][-1])
    #ax.set_ylim(bottom=0.0, top=np.max(data[indx_eq_search_smart][0]))
    ax.set_xlabel('Number of stored states', **custom_font)
    ax.set_ylabel('CPU time (sec)', **custom_font)
    ax.plot(data[smart][tree_size], data[smart][time_estimate], label="Estimate time", color='darkcyan')
    ax.plot(data[smart][tree_size], data[smart][time_search], label="Search time", color='darkorchid')
    ax.plot(data[smart][tree_size], 1e-6 * data[smart][tree_size], label="Linear growth", linestyle=':', color='rosybrown')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid(color='lightgray', linestyle=':', linewidth=1)
    ax.legend(loc='upper left')
    plt.tight_layout()
    fig.savefig(folder_general + '/results/cpu-time-smart-est-search' + tag + '.pdf')

def plot_learning_statuses(statuses, text):

    time = np.linspace(0, dt * nsteps / minute, nsteps)
    learnings = np.zeros(nsteps)

    for step in range(nsteps):
        (learnings[step], total_statuses) = count_trainings(status[step])

    # Plot change of pH wrt the space coordinates
    fig, ax = myfig(6, 3)
    ax.set_xlim(left=0, right=time[-1])
    ax.set_ylim(bottom=0, top=7.1)
    ax.set_xlabel('Time (minute)', **custom_font)
    ax.plot(time, learnings, color='steelblue')
    ax.grid(color='lightgray', linestyle=':', linewidth=1)
    ax.set_title('Number of On Demand Trainings per Time Step')
    #plt.rcParams.update({'font.size': 11})
    total_text = text[0] + '\n' + text[1] + '\n\n' + text[2]
    plt.tight_layout()
    ax.text(0.5, 0.9, text[0], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.5, 0.8, text[1], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.5, 0.65, text[2], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    fig.savefig(folder_general + '/results/on-demand-trainings' + tag + '.pdf')

def count_trainings(status):
    counter = [st for st in status if st == 0]
    return (len(counter), len(status))

if __name__ == '__main__':

    profiling = []
    status    = []

    make_results_folders()

    # Collect files with results corresponding to smart or reference (classical) solver
    files_smart = [file for file in sorted(os.listdir(folder_smart)) if ("profiling" not in file and "statuses" not in file)]
    files_class = [file for file in sorted(os.listdir(folder_class)) if ("profiling" not in file)]

    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    status = np.loadtxt(folder_smart + '/statuses.txt')

    # Collect files with profiling data
    profiling   = profiling + [(file, 1) for file in sorted(os.listdir(folder_smart)) if "profiling" in file]
    profiling   = profiling + [(file, 0) for file in sorted(os.listdir(folder_class)) if "profiling" in file]

    # Count the percentage of the trainings needed
    training_counter = count_trainings(np.array(status).flatten())
    title = "%2.2f percent is training (%d out of %d cells)" % (100 * training_counter[0] / training_counter[1], training_counter[0], training_counter[1])
    text = ["Number of chemical equilibrium trainings: %d" % training_counter[0],
            "Number of smart chemical equilibrium estimations: %d" % (training_counter[1] - training_counter[0]),
            ("Percentage of smart chemical equlibrium predictions: %2.2f" % (100 - 100 * training_counter[0] / training_counter[1])) + "$\%$"]
    plot_learning_statuses(status, text)
    print(title)
    # Plot profiling results
    plot_cpu_times(profiling, title)

    # Plot comparison of the chemistry of the conventional and smart solvers
    #for file in files_smart: compare_chemistry(file, status);
    selected_steps = [10, 60, 120, 240, 480, 960, 1200, 2400]
    Parallel(n_jobs=16)(delayed(compare_chemistry)(file, status, selected_steps) for file in files_smart)

    # Save selected plots into results folder

    for step in selected_steps:
        copy_cmd = ('cp ' + folder_general + '/figures/{0}/%s.pdf ' + folder_general +
                    '/results/{0}-%s-%dmin.pdf') % (str(step).zfill(fillz), str(step).zfill(fillz), (step)*dt/minute)
        os.system(copy_cmd.format('calcite-dolomite'))
        os.system(copy_cmd.format('ph'))
        os.system(copy_cmd.format('aqueous-species'))

    # Make  videos
    # --------------------------------------------------------------------------
    # Define the command for making videos

    ffmpeg_cmd_caldol = 'ffmpeg -loglevel quiet -y -r 480 -i ' + folder_general + '/figures/calcite-dolomite/%03d.png ' \
                 '-codec:v mpeg4 -flags:v +qscale -global_quality:v 0 ' \
                 + folder_general + '/videos/calcite-dolomite' + tag + '.mp4'
    ffmpeg_cmd_ph = 'ffmpeg -loglevel quiet -y -r 120 -i ' + folder_general + '/figures/ph/%03d.png ' \
                 '-codec:v mpeg4 -flags:v +qscale -global_quality:v 0 ' \
                 + folder_general + '/videos/ph' + tag + '.mp4'
    ffmpeg_cmd_aq = 'ffmpeg -loglevel quiet -y -r 120 -i ' + folder_general + '/figures/aqueous-species/%03d.png ' \
                                                                              '-codec:v mpeg4 -flags:v +qscale -global_quality:v 0 ' \
                    + folder_general + '/videos/aqueous-species' + tag + '.mp4'
    '''
    ffmpeg_cmd = 'ffmpeg -loglevel quiet -y -r 240 -i ' + folder_general + '/figures/{0}/%d.png ' \
                 '-codec:v mpeg4 -flags:v +qscale -global_quality:v 0 ' \
                 + folder_general + '/videos/{0}' + tag + '.mp4'
    '''
    os.system(ffmpeg_cmd_caldol)
    os.system(ffmpeg_cmd_ph)
    os.system(ffmpeg_cmd_aq)
