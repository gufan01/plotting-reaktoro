import numpy as np
import matplotlib.pyplot as plt
import os
from natsort import natsorted
import matplotlib as mpl
import json

import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", module="matplotlib")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
#mpl.style.use('v2.0')

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

#folder_smart   = "results-clusters-partition-major-species-dt-1800-ncells-100-nsteps-10000-eqreltol-3.0e-03-smart"
#folder_class   = "results-clusters-partition-major-species-dt-1800-ncells-100-nsteps-10000-eqreltol-3.0e-03-reference"
#folder_general = "plot-results-clusters-partition-major-species-dt-1800-ncells-100-nsteps-1000-eqreltol-3.0e-03"
#plot_at_selected_steps = [1, 10, 20, 40, 80, 100, 120, 240, 480, 960, 1200, 2400, 3600, 4800, 7200]  # the time steps at which the results are plotted
#plot_at_selected_steps = [1,  10, 20, 40, 80, 100, 120, 240, 480, 960, 1000]  # the time steps at which the results are plotted
#nsteps = 10000    # the number of steps in the reactive transport simulation

# results for the paper
folder_smart   = "results-commit-fe311182-with-fixes-dt-1800-ncells-100-nsteps-10000-reltol-1.0e-03-pitzer-full-smart"
folder_class   = "results-commit-fe311182-with-fixes-dt-1800-ncells-100-nsteps-10000-pitzer-full-reference"
folder_general = "plots-commit-fe311182-with-fixes"
#folder_smart   = "results-custering-primary-species-paper-dt-1800-ncells-100-nsteps-10000-reltol-1.0e-03-pitzer-full-smart"
#folder_class   = "results-custering-primary-species-paper-dt-1800-ncells-100-nsteps-10000-pitzer-full-reference"
#folder_general = "plots-results-custering-primary-species-paper-dt-1800-ncells-100-nsteps-10000-reltol-1.0e-03-pitzer-full-smart"
#folder_smart   = "results-custering-primary-species-paper-skipping-1e-12-dt-1800-ncells-100-nsteps-10000-reltol-1.0e-03-hkf-full-smart"
#folder_class   = "results-custering-primary-species-paper-skipping-1e-12-dt-1800-ncells-100-nsteps-10000-hkf-full-reference"
#folder_general = "plots-results-custering-primary-species-paper-skipping-1e-12-dt-1800-ncells-100-nsteps-10000-reltol-1.0e-03-hkf-full-smart"
plot_at_selected_steps = [1, 10, 20, 40, 80, 100, 120, 240, 480, 960, 1200, 2400, 3600, 4800, 7200]  # the time steps at which the results are plotted
#plot_at_selected_steps = [500, 800, 864, 912, 1000, 1100, 1200, 1300]  # the time steps at which the results are plotted
#plot_at_selected_steps = [1,  10, 20, 40, 80, 100, 120, 240, 480, 960, 1000]  # the time steps at which the results are plotted
nsteps = 10000    # the number of steps in the reactive transport simulation

res1 = np.array([ 4.41836e-27, 3.30966e-24, 3.31731e-24, 8.52991e-27, 7.66111e-27, 6.99508e-27, 4.13614e-24, 6.62122e-24, 6.70269e-27, 3.31277e-24, 8.27225e-25, 6.59304e-27, 8.33643e-25, 8.30854e-25, 8.27644e-25, 3.31573e-24, 2.06865e-23, 7.45191e-24, 4.1363e-24, 4.13707e-24, 3.30916e-24, 8.27825e-25, 6.62112e-24, 6.62136e-24, 6.62091e-24, 6.6213e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24, 6.62131e-24])
res10 =  np.array([ 8.32447e-25, 4.14941e-24, 3.32773e-24, 4.16427e-24, 4.14393e-24, 1.66787e-27, 3.30976e-24, 4.13686e-24, 8.35746e-25, 3.8273e-27, 4.16176e-24, 8.27967e-25, 6.6322e-24, 4.13707e-24, 3.06095e-23, 1.40786e-23, 6.61751e-24, 6.61766e-24, 1.65474e-23, 4.13956e-24, 4.13612e-24, 6.63391e-24, 8.43557e-25, 3.72268e-23, 3.31256e-24, 3.31217e-24, 3.30973e-24, 1.4066e-23, 6.63449e-24, 1.65444e-23, 1.07573e-23, 6.61876e-24, 3.80043e-27, 6.62427e-24, 4.13953e-24, 6.62114e-24, 4.13972e-24, 4.13955e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.1394e-24, 4.13921e-24, 4.13921e-24])
res2400 =  np.array([ 3.19215e-28, 1.58241e-27, 4.90336e-27, 8.43116e-25, 3.31013e-23, 8.44418e-25, 1.32482e-23, 1.65571e-23, 8.54e-25, 4.13798e-24, 8.31515e-25, 8.48677e-25, 4.67183e-27, 4.24045e-24, 1.61592e-26, 8.44661e-25, 6.63381e-24, 4.20111e-24, 8.15923e-26, 8.9222e-25, 1.41302e-23, 3.32604e-24, 3.32508e-24, 4.20289e-24, 6.54351e-26, 3.37435e-24, 1.41272e-23, 3.37389e-24, 6.48947e-26, 8.44385e-25, 8.44583e-25, 8.43906e-25, 4.20075e-24, 4.15288e-24, 8.92883e-25, 4.15261e-24, 3.37636e-24, 4.20224e-24, 4.15242e-24, 8.92965e-25, 6.6829e-24, 8.1928e-26, 3.37383e-24, 8.55815e-25, 8.13415e-26, 4.24074e-24, 3.37356e-24, 4.2186e-24, 1.65613e-23, 3.38982e-24, 4.15208e-24, 8.13182e-26, 3.32783e-24, 1.33027e-23, 6.68227e-24, 9.10202e-25, 3.39001e-24, 3.31148e-23, 8.92176e-25, 1.4127e-23, 8.93854e-25, 4.20122e-24, 9.09739e-25, 4.15334e-24, 4.2018e-24, 4.1916e-24, 3.37373e-24, 3.32537e-24, 8.92208e-25, 1.33161e-23, 1.40786e-23, 4.21698e-24, 4.20134e-24, 8.92362e-25, 8.92094e-25, 1.41269e-23, 3.32621e-24, 4.21725e-24, 1.31522e-26, 1.41432e-23, 3.33524e-24, 5.555e-26, 8.92748e-25, 2.07453e-23, 3.39109e-24, 3.39052e-24, 4.20219e-24, 6.48725e-26, 3.45423e-24, 4.21732e-24, 3.37422e-24, 3.32531e-24, 4.20266e-24, 1.656e-23, 4.15474e-24, 3.39057e-24, 3.32521e-24, 9.09042e-25, 4.20155e-24, 4.20155e-24])
os.system('mkdir -p ' + folder_general)


import os
from matplotlib import font_manager as fm, rcParams
#fpath = os.path.join(rcParams["datapath"], "/usr/share/fonts/truetype/ebgaramond/EBGaramond12-Regular.ttf")
#fpath = os.path.join(rcParams["datapath"], "/usr/share/fonts/truetype/ebgaramond/EBGaramond12-Bold.ttf")
#fpath = os.path.join(rcParams["datapath"], "/usr/share/fonts/truetype/ebgaramond/TeX-Gyre-Adventor.ttf")
fpath = os.path.join(rcParams["datapath"], "/home/skyas/Dropbox/texgyreadventor-regular.otf")
prop = fm.FontProperties(fname=fpath)
prop.set_size(14)


from matplotlib import font_manager
prop2 = font_manager.findfont('TeX Gyre Adventor', rebuild_if_missing=True)

fillz = len(str(123))

# Indices of the loaded data to plot
indx_ph = 0
indx_Hcation = 1
indx_Cacation = 2
indx_Mgcation = 3
indx_HCO3anion = 4
indx_CO2aq = 5
indx_calcite = 6
indx_dolomite = 7

indx_CO3anion = 8
indx_CaClcation = 9
indx_CaHCO3cation = 10
indx_MgClcation = 11
indx_MgHCO3caanion = 12
indx_OHanion = 13

indx_C = 14
indx_Ca = 15
indx_Cl = 16
indx_H = 17
indx_Mg = 18
indx_Na = 19
indx_O = 20
indx_Si = 21
indx_Z = 22

# Plotting params
circ_area = 6 ** 2
# zoom = 0.5
zoom = 0.5
custom_font = { }
time_steps = np.linspace(0, nsteps, nsteps)

xcells = np.linspace(xl, xr, ncells)  # the x-coordinates of the plots

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'TeX Gyre Adventor'
mpl.rcParams['font.style'] = 'normal'
mpl.rcParams['font.size'] = 14
#mpl.rcParams['axes.linewidth'] = 1.5
#mpl.rcParams['legend.fontsize'] =
#mpl.rcParams['legend.fontsize'] = 'medium'
#tick label spacing and tick width
# mpl.rcParams['xtick.major.pad'] = 4
# mpl.rcParams['ytick.major.pad'] = 5
# mpl.rcParams['xtick.major.width'] = 1
# mpl.rcParams['ytick.major.width'] = 1
#legend style
# mpl.rcParams['legend.frameon'] = True
# mpl.rcParams['legend.numpoints'] = 3
#mpl.rcParams['backend'] = 'pdf'
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

def line_error(color):
    return {'linestyle': ':', 'marker': 'D', 'color': color, 'zorder': 1, 'linewidth': 0.5, 'markersize': 4}

def titlestr(t):
    d = int(t / day)                 # The number of days
    h = int(int(t % day) / hour)     # The number of remaining hours
    m = int(int(t % hour) / minute)  # The number of remaining minutes
    return '{:>3d}d {:>2}h {:>2}m'.format(int(d), str(int(h)).zfill(2), str(int(m)).zfill(2))

def plot_figures_ph():

    for i in plot_at_selected_steps:
        #print("On pH figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        data_class_ph = data_class[indx_ph]
        data_smart_ph = data_smart[indx_ph]
        plt.axes(xlim=(-0.01, 0.501), ylim=(2.5, 12.0))
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('pH', fontproperties=prop)
        plt.title(titlestr(t), fontproperties=prop,fontsize=20)
        plt.plot(xcells, data_class_ph, label='pH', **line('teal'))
        plt.plot(xcells[status[i-1]==0], data_smart_ph[status[i-1]==0], 'o', **line_empty_marker('teal'))
        plt.plot(xcells[status[i-1]==1], data_smart_ph[status[i-1]==1], 'o', **line_filled_marker('teal'))
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='lower right', prop=prop)
        plt.savefig(folder_general + '/pH-{}.pdf'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_calcite_dolomite():

    for i in plot_at_selected_steps:
        #print("On calcite-dolomite figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T

        data_class_calcite, data_class_dolomite = data_class[indx_calcite], data_class[indx_dolomite]
        data_smart_calcite, data_smart_dolomite = data_smart[indx_calcite], data_smart[indx_dolomite]
        plt.axes(xlim=(-0.01, 0.501), ylim=(-0.1, 2.1))
        plt.ylabel('Mineral Volume [%vol]', fontproperties=prop)
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.title(titlestr(t), fontproperties=prop,fontsize=20)
        plt.plot(xcells, data_class_calcite * 100/(1 - phi), label='Calcite', **line('C0'))
        plt.plot(xcells, data_class_dolomite * 100/(1 - phi), label='Dolomite', **line('C1'))
        plt.plot(xcells[status[i-1]==0], data_smart_calcite[status[i-1]==0] * 100/(1 - phi), 'o', **line_empty_marker('C0'))
        plt.plot(xcells[status[i-1]==1], data_smart_calcite[status[i-1]==1] * 100/(1 - phi), 'o', **line_filled_marker('C0'))
        plt.plot(xcells[status[i-1]==0], data_smart_dolomite[status[i-1]==0] * 100/(1 - phi), 'o', **line_empty_marker('C1'))
        plt.plot(xcells[status[i-1]==1], data_smart_dolomite[status[i-1]==1] * 100/(1 - phi), 'o', **line_filled_marker('C1'))
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='center right', prop=prop)
        plt.savefig(folder_general + '/calcite-dolomite-{}.pdf'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_calcite_dolomite_moles():

    for i in plot_at_selected_steps:
        #print("On calcite-dolomite figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T

        data_class_calcite, data_class_dolomite = data_class[indx_mol_calcite], data_class[indx_mol_dolomite]
        data_smart_calcite, data_smart_dolomite = data_smart[indx_mol_calcite], data_smart[indx_mol_dolomite]
        plt.axes(xlim=(-0.01, 0.501))
        plt.ylabel('Concentration [mol/m3]')
        plt.xlabel('Distance [m]')
        plt.title(titlestr(t))
        plt.plot(xcells, data_class_calcite, label='Calcite', **line('C0'))
        plt.plot(xcells, data_class_dolomite, label='Dolomite', **line('C1'))
        plt.plot(xcells[status[i-1]==0], data_smart_calcite[status[i-1]==0], 'o', **line_empty_marker('C0'))
        plt.plot(xcells[status[i-1]==1], data_smart_calcite[status[i-1]==1], 'o', **line_filled_marker('C0'))
        plt.plot(xcells[status[i-1]==0], data_smart_dolomite[status[i-1]==0], 'o', **line_empty_marker('C1'))
        plt.plot(xcells[status[i-1]==1], data_smart_dolomite[status[i-1]==1], 'o', **line_filled_marker('C1'))
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='center right')
        plt.savefig(folder_general + '/calcite-dolomite-moles-{}.pdf'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_log10_aqueous_species():

    for i in plot_at_selected_steps:
        #print("On aqueous-species figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T
        indx_ph = 0

        data_class_cacation  = data_class[indx_Cacation]
        data_class_mgcation  = data_class[indx_Mgcation]
        data_class_hco3anion = data_class[indx_HCO3anion]
        data_class_co2aq     = data_class[indx_CO2aq]
        data_class_hcation   = data_class[indx_Hcation]
        data_class_CO3anion   = data_class[indx_CO3anion]
        data_class_CaHCO3cation = data_class[indx_CaHCO3cation]
        data_class_CaClcation = data_class[indx_CaClcation]
        data_class_MgClcation = data_class[indx_MgClcation]
        data_class_MgHCO3caanion = data_class[indx_MgHCO3caanion]
        data_class_OHanion       = data_class[indx_OHanion]

        data_smart_cacation  = data_smart[indx_Cacation]
        data_smart_mgcation  = data_smart[indx_Mgcation]
        data_smart_hco3anion = data_smart[indx_HCO3anion]
        data_smart_co2aq     = data_smart[indx_CO2aq]
        data_smart_hcation   = data_smart[indx_Hcation]
        data_smart_CO3anion   = data_smart[indx_CO3anion]
        data_smart_CaHCO3cation  = data_smart[indx_CaHCO3cation]
        data_smart_CaClcation = data_smart[indx_CaClcation]
        data_smart_MgClcation = data_smart[indx_MgClcation]
        data_smart_MgHCO3caanion = data_smart[indx_MgHCO3caanion]
        data_smart_OHanion       = data_smart[indx_OHanion]

        plt.axes(xlim=(-0.01, 0.501), ylim=(-10, 1))
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('Concentration [molal]', fontproperties=prop)
        plt.title(titlestr(t), fontproperties=prop,fontsize=20)
        plt.plot(xcells, np.log10(data_class_cacation), label='Ca+2', **line('C0'))[0],
        plt.plot(xcells, np.log10(data_class_mgcation), label='Mg+2', **line('C1'))[0],
        plt.plot(xcells, np.log10(data_class_hco3anion), label='HCO3-',**line('C2'))[0],
        plt.plot(xcells, np.log10(data_class_co2aq), label='CO2(aq)',**line('C3'))[0],
        plt.plot(xcells, np.log10(data_class_hcation), label='H+', **line('C4'))[0],
        plt.plot(xcells, np.log10(data_class_CO3anion), label='CO3-2', **line('C5'))[0],
        plt.plot(xcells, np.log10(data_class_CaHCO3cation), label='Ca(HCO3)+', **line('C6'))[0],
        plt.plot(xcells, np.log10(data_class_CaClcation), label='CaCl+',**line('C7'))[0],
        plt.plot(xcells, np.log10(data_class_MgClcation), label='MgCl+',**line('C8'))[0],
        plt.plot(xcells, np.log10(data_class_MgHCO3caanion), label='Mg(HCO3))+',**line('C9'))[0],
        plt.plot(xcells, np.log10(data_class_OHanion), label='OH-', **line('darkviolet'))[0],

        plt.plot(xcells[status[i-1]==0], np.log10(data_smart_cacation[status[i-1]==0]), 'o', **line_empty_marker('C0'))[0],
        plt.plot(xcells[status[i-1]==1], np.log10(data_smart_cacation[status[i-1]==1]), 'o', **line_filled_marker('C0'))[0],
        plt.plot(xcells[status[i-1]==0], np.log10(data_smart_mgcation[status[i-1]==0]), 'o', **line_empty_marker('C1'))[0],
        plt.plot(xcells[status[i-1]==1], np.log10(data_smart_mgcation[status[i-1]==1]), 'o', **line_filled_marker('C1'))[0],
        plt.plot(xcells[status[i-1]==0], np.log10(data_smart_hco3anion[status[i-1]==0]), 'o', **line_empty_marker('C2'))[0],
        plt.plot(xcells[status[i-1]==1], np.log10(data_smart_hco3anion[status[i-1]==1]), 'o', **line_filled_marker('C2'))[0],
        plt.plot(xcells[status[i-1]==0], np.log10(data_smart_co2aq[status[i-1]==0]), 'o', **line_empty_marker('C3'))[0],
        plt.plot(xcells[status[i-1]==1], np.log10(data_smart_co2aq[status[i-1]==1]), 'o', **line_filled_marker('C3'))[0],
        plt.plot(xcells[status[i-1]==0], np.log10(data_smart_hcation[status[i-1]==0]), 'o', **line_empty_marker('C4'))[0],
        plt.plot(xcells[status[i-1]==1], np.log10(data_smart_hcation[status[i-1]==1]), 'o', **line_filled_marker('C4'))[0],
        plt.plot(xcells[status[i-1]==0], np.log10(data_smart_CO3anion[status[i-1]==0]), 'o', **line_empty_marker('C5'))[0],
        plt.plot(xcells[status[i-1]==1], np.log10(data_smart_CO3anion[status[i-1]==1]), 'o', **line_filled_marker('C5'))[0], \
        plt.plot(xcells[status[i-1]==0], np.log10(data_smart_CaHCO3cation[status[i-1]==0]), 'o', **line_empty_marker('C6'))[0],
        plt.plot(xcells[status[i-1]==1], np.log10(data_smart_CaHCO3cation[status[i-1]==1]), 'o', **line_filled_marker('C6'))[0],
        plt.plot(xcells[status[i-1]==0], np.log10(data_smart_CaClcation[status[i-1]==0]), 'o', **line_empty_marker('C7'))[0],
        plt.plot(xcells[status[i-1]==1], np.log10(data_smart_CaClcation[status[i-1]==1]), 'o', **line_filled_marker('C7'))[0],
        plt.plot(xcells[status[i-1]==0], np.log10(data_smart_MgClcation[status[i-1]==0]), 'o', **line_empty_marker('C8'))[0],
        plt.plot(xcells[status[i-1]==1], np.log10(data_smart_MgClcation[status[i-1]==1]), 'o', **line_filled_marker('C8'))[0],
        plt.plot(xcells[status[i-1]==0], np.log10(data_smart_MgHCO3caanion[status[i-1]==0]), 'o', **line_empty_marker('C9'))[0],
        plt.plot(xcells[status[i-1]==1], np.log10(data_smart_MgHCO3caanion[status[i-1]==1]), 'o', **line_filled_marker('C9'))[0],
        plt.plot(xcells[status[i-1]==0], np.log10(data_smart_OHanion[status[i-1]==0]), 'o', **line_empty_marker('darkviolet'))[0],
        plt.plot(xcells[status[i-1]==1], np.log10(data_smart_OHanion[status[i-1]==1]), 'o', **line_filled_marker('darkviolet'))[0],
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        plt.legend(loc='upper right', prop=prop)
        plt.savefig(folder_general + '/aqueous-species-log10-{}.pdf'.format(i))
        plt.tight_layout()
        plt.close()

def plot_figures_aqueous_species():

    #for i in plot_at_selected_steps:
    for i in [1, 20, 2400]:
        #print("On aqueous-species figure at time step: {}".format(i))
        t = i * dt
        filearray_class = np.loadtxt(folder_class + '/' + files_class[i-1], skiprows=1)
        filearray_smart = np.loadtxt(folder_smart + '/' + files_smart[i-1], skiprows=1)
        data_class = filearray_class.T
        data_smart = filearray_smart.T

        data_class_cacation  = data_class[indx_Cacation]
        data_class_mgcation  = data_class[indx_Mgcation]
        data_class_hco3anion = data_class[indx_HCO3anion]
        data_class_co2aq     = data_class[indx_CO2aq]
        data_class_hcation   = data_class[indx_Hcation]
        data_class_CO3anion   = data_class[indx_CO3anion]
        data_class_CaHCO3cation = data_class[indx_CaHCO3cation]
        data_class_CaClcation = data_class[indx_CaClcation]
        data_class_MgClcation = data_class[indx_MgClcation]
        data_class_MgHCO3caanion = data_class[indx_MgHCO3caanion]
        data_class_OHanion       = data_class[indx_OHanion]

        data_smart_cacation  = data_smart[indx_Cacation]
        data_smart_mgcation  = data_smart[indx_Mgcation]
        data_smart_hco3anion = data_smart[indx_HCO3anion]
        data_smart_co2aq     = data_smart[indx_CO2aq]
        data_smart_hcation   = data_smart[indx_Hcation]
        data_smart_CO3anion   = data_smart[indx_CO3anion]
        data_smart_CaHCO3cation  = data_smart[indx_CaHCO3cation]
        data_smart_CaClcation = data_smart[indx_CaClcation]
        data_smart_MgClcation = data_smart[indx_MgClcation]
        data_smart_MgHCO3caanion = data_smart[indx_MgHCO3caanion]
        data_smart_OHanion       = data_smart[indx_OHanion]

        plt.axes(xlim=(-0.01, 0.501), ylim=(1e-11, 5))
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('Concentration [molal]', fontproperties=prop)
        plt.yscale('log')
        plt.title(titlestr(t), fontproperties=prop,fontsize=20)
        plt.plot(xcells, (data_class_cacation), label=r'$Ca^{2+}$', **line('C0'))[0],
        plt.plot(xcells, (data_class_mgcation), label=r'$Mg^{2+}$', **line('C1'))[0],
        plt.plot(xcells, (data_class_hco3anion), label=r'$HCO_3^-$',**line('C2'))[0],
        plt.plot(xcells, (data_class_co2aq), label=r'$CO_2(aq)$',**line('C3'))[0],
        plt.plot(xcells, (data_class_hcation), label=r'$H^+$', **line('C4'))[0],
        plt.plot(xcells, (data_class_CO3anion), label=r'$CO_3^{2-}$', **line('C5'))[0],
        plt.plot(xcells, (data_class_CaHCO3cation), label=r'$Ca(HCO_3)^+$', **line('C6'))[0],
        plt.plot(xcells, (data_class_CaClcation), label=r'$CaCl^+$',**line('C7'))[0],
        plt.plot(xcells, (data_class_MgClcation), label=r'$MgCl^+$',**line('C8'))[0],
        plt.plot(xcells, (data_class_MgHCO3caanion), label=r'$Mg(HCO_3))^+$',**line('C9'))[0],
        plt.plot(xcells, (data_class_OHanion), label=r'$OH^-$', **line('darkviolet'))[0],

        plt.plot(xcells[status[i-1]==0], data_smart_cacation[status[i-1]==0], 'o', **line_empty_marker('C0'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_cacation[status[i-1]==1], 'o', **line_filled_marker('C0'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_mgcation[status[i-1]==0], 'o', **line_empty_marker('C1'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_mgcation[status[i-1]==1], 'o', **line_filled_marker('C1'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_hco3anion[status[i-1]==0], 'o', **line_empty_marker('C2'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_hco3anion[status[i-1]==1], 'o', **line_filled_marker('C2'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_co2aq[status[i-1]==0], 'o', **line_empty_marker('C3'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_co2aq[status[i-1]==1], 'o', **line_filled_marker('C3'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_hcation[status[i-1]==0], 'o', **line_empty_marker('C4'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_hcation[status[i-1]==1], 'o', **line_filled_marker('C4'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_CO3anion[status[i-1]==0], 'o', **line_empty_marker('C5'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_CO3anion[status[i-1]==1], 'o', **line_filled_marker('C5'))[0], \
        plt.plot(xcells[status[i-1]==0], data_smart_CaHCO3cation[status[i-1]==0], 'o', **line_empty_marker('C6'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_CaHCO3cation[status[i-1]==1], 'o', **line_filled_marker('C6'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_CaClcation[status[i-1]==0], 'o', **line_empty_marker('C7'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_CaClcation[status[i-1]==1], 'o', **line_filled_marker('C7'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_MgClcation[status[i-1]==0], 'o', **line_empty_marker('C8'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_MgClcation[status[i-1]==1], 'o', **line_filled_marker('C8'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_MgHCO3caanion[status[i-1]==0], 'o', **line_empty_marker('C9'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_MgHCO3caanion[status[i-1]==1], 'o', **line_filled_marker('C9'))[0],
        plt.plot(xcells[status[i-1]==0], data_smart_OHanion[status[i-1]==0], 'o', **line_empty_marker('darkviolet'))[0],
        plt.plot(xcells[status[i-1]==1], data_smart_OHanion[status[i-1]==1], 'o', **line_filled_marker('darkviolet'))[0],
        plt.plot([], [], 'o', label='Smart Prediction', **line_filled_marker('black'))
        plt.plot([], [], 'o', label='Learning', **line_empty_marker('black'))
        prop.set_size(12)
        plt.legend(loc='upper right', prop=prop)
        plt.savefig(folder_general + '/aqueous-species-{}.pdf'.format(i))
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

    plt.xlabel('Time Step', fontproperties=prop)
    plt.ylabel('Computing Cost [μs]', fontproperties=prop)
    plt.xlim(left=0, right=nsteps)
    plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], timings_equilibrium_class[0:nsteps:step], label="Chemical Equilibrium (Conventional)", color='C0', linewidth=1.5)
    plt.plot(time_steps[0:nsteps:step], timings_equilibrium_smart[0:nsteps:step], label="Chemical Equilibrium (Smart)", color='C1', linewidth=1.5, alpha=1.0)
    plt.plot(time_steps[0:nsteps:step], timings_equilibrium_smart_ideal[0:nsteps:step], label="Chemical Equilibrium (Smart, Ideal)", color='C3', linewidth=1.5, alpha=1.0)
    plt.plot(time_steps[0:nsteps:step], timings_transport[0:nsteps:step], label="Transport", color='C2', linewidth=1.5, alpha=1.0)
    leg = plt.legend(loc='right', bbox_to_anchor=(0.5, 0.4, 0.5, 0.5), prop=prop)
    #for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/computing-costs-nolegend-with-smart-ideal.pdf')
    #plt.savefig(folder_general + '/computing-costs-nolegend-with-smart.pdf')
    #plt.savefig(folder_general + '/computing-costs-nolegend.pdf')

    plt.close()

def plot_on_demand_learnings_total_initial_steps(nsteps_initial):

    step = 1

    # Load the status data, where 0 stands for conventional learning and 1 for smart prediction
    with open(folder_smart + '/analysis-smart.json') as read_file:
        data = json.load(read_file)

    status_learnings = data.get('cells_where_learning_was_required_at_step')

    # Collect the number of learnings on each step
    learnings = [len(x) for x in status_learnings]
    total_learnings = [np.sum(learnings[0:i]) for i in range(0, nsteps)]

    print(f"number of learning occuring on first {nsteps_initial} stesp: ", learnings[0:nsteps_initial])
    plt.xlabel('Time Step', fontproperties=prop)
    plt.ylabel('Accumulated on-demand learnings (over all steps)', fontproperties=prop)
    plt.ylim(bottom=0, top=total_learnings[nsteps_initial] + 1)
    plt.ticklabel_format(style='plain', axis='x')
    xint = range(0, int(time_steps[nsteps_initial]))
    plt.xticks(xint)
    plt.plot(xint, total_learnings[0:nsteps_initial:step], color='C0', linewidth=1.5, marker='D')
    plt.tight_layout()
    plt.savefig(folder_general + '/on-demand-learning-total-initial-steps-'+ str(nsteps_initial) + '.pdf')
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
    plt.ticklabel_format(style='plain', axis='y')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(learnings, color='C0', linewidth=1.5)
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


    plt.xlabel('Time Step', fontproperties=prop)
    plt.ylabel('Accumulated on-demand learnings', fontproperties=prop)
    plt.xlim(left=0, right=nsteps)
    plt.ylim(bottom=0, top=np.max(accum_learnings)+5)
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], accum_learnings[0:nsteps:step], color='C0', linewidth=1.5)
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

    plt.xlabel('Time Step', fontproperties=prop)
    plt.ylabel('Computing Cost [μs]', fontproperties=prop)
    #plt.xscale('log')
    plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], timings_search[0:nsteps:step], label="Search Time", color='C0', linewidth=1.5)
    plt.plot(time_steps[0:nsteps:step], timings_estimate[0:nsteps:step], label="Estimate Time", color='C1', linewidth=1.5)

    leg = plt.legend(loc='lower right', prop=prop)
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
                    (np.array(timing_smart.get('smart_equilibrium')) \
                     - np.array(timing_smart.get('smart_equilibrium_nearest_neighbor_search')) \
                     - np.array(timing_smart.get('smart_equilibrium_storage')))

    print("conv CPU         = ", np.sum(np.array(timing_class.get('equilibrium'))))
    print("odml CPU         = ", np.sum(np.array(timing_smart.get('smart_equilibrium'))))
    print("odml CPU (ideal) = ", np.sum(np.array(timing_smart.get('smart_equilibrium')) \
                                       - np.array(timing_smart.get('smart_equilibrium_nearest_neighbor_search')) \
                                       - np.array(timing_smart.get('smart_equilibrium_storage'))))

    print("speedup         = ", np.sum(np.array(timing_class.get('equilibrium'))) / np.sum(np.array(timing_smart.get('smart_equilibrium'))))
    print("speedup (ideal) = ", np.sum(np.array(timing_class.get('equilibrium'))) / (np.sum(np.array(timing_smart.get('smart_equilibrium')) \
                                                                                           - np.array(timing_smart.get('smart_equilibrium_nearest_neighbor_search')) \
                                                                                           - np.array(timing_smart.get('smart_equilibrium_storage')))))

    #print("speedup         = ", np.average(speedup))
    #print("speedup (ideal) = ", np.average(speedup_ideal))

    plt.xlim(left=0, right=nsteps)
    plt.xlabel('Time Step', fontproperties=prop)
    plt.ylabel('Speedup [-]', fontproperties=prop)
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], speedup_ideal[0:nsteps:step], label="Conventional vs. Smart (Ideal search)",
             color='C2', linewidth=1.5, alpha=1.0)
    plt.plot(time_steps[0:nsteps:step], speedup[0:nsteps:step], label="Conventional vs. Smart ", color='C0', linewidth=1.5)

    leg = plt.legend(loc='upper right', prop=prop)
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

    #for i in plot_at_selected_steps:
    for i in [1, 20, 2400]:

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
        plt.savefig(folder_general + f'/error-rel-calcite-dolomite-{i}.pdf')
        plt.tight_layout()
        plt.close()

        plt.axes(xlim=(-0.01, 0.501))
        if i == 1: plt.axes(xlim=(-0.01, 0.501), ylim=[-5e-21, 1.1e-19])
        elif i == 20: plt.axes(xlim=(-0.01, 0.501), ylim=[-7e-4, 3.5e-2])
        elif i == 2400: plt.axes(xlim=(-0.01, 0.501), ylim=[-4e-4, 1.6e-2])

        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('Relative Error [-]', fontproperties=prop)
        #plt.yscale('log')
        plt.title(titlestr(t), fontproperties=prop, fontsize=20)
        plt.plot(xcells, rel_error_cacation, label=r'$Ca^{2+}$', **line_error('C0'))
        plt.plot(xcells, rel_error_mgcation, label=r'$Mg^{2+}$', **line_error('C1'))
        plt.plot(xcells, rel_error_hco3anion, label=r'$HCO_3^-$', **line_error('C2'))
        plt.plot(xcells, rel_error_co2aq, label=r'$CO_2(aq)$', **line_error('C3'))
        plt.plot(xcells, rel_error_hcation, label=r'$H^+$', **line_error('C4'))
        plt.plot(xcells, rel_error_CO3anion, label=r'$CO_3^{2-}$', **line_error('C5'))[0],
        plt.plot(xcells, rel_error_CaHCO3cation, label=r'$Ca(HCO_3)^+$', **line_error('C6'))[0],
        plt.plot(xcells, rel_error_CaClcation, label=r'$CaCl^+$',**line_error('C7'))[0],
        plt.plot(xcells, rel_error_MgClcation, label=r'$MgCl^+$',**line_error('C8'))[0],
        plt.plot(xcells, rel_error_MgHCO3caanion, label=r'$Mg(HCO_3))^+$',**line_error('C9'))[0],
        plt.plot(xcells, rel_error_OHanion, label=r'$OH^-$', **line_error('darkviolet'))[0],
        prop.set_size(12)
        if i == 1 or i == 20: plt.legend(loc='upper right', prop=prop)
        elif i == 2400: plt.legend(loc='upper left', prop=prop)
        plt.savefig(folder_general + f'/error-rel-aqueous-species-{i}.pdf')
        plt.tight_layout()
        plt.close()

        '''
        t = i * dt
        min_val = min(min(abs_error_calcite[0:int(ncells/2)] * 100 / (1 - phi)), min(abs_error_dolomite[0:int(ncells/2)] * 100 / (1 - phi)))
        max_val = 1.1*max(max(abs_error_calcite[0:int(ncells/2)] * 100 / (1 - phi)), max(abs_error_dolomite[0:int(ncells/2)] * 100 / (1 - phi)))
        #plt.axes(xlim=(-0.01, 0.501), ylim=(-0.0000001,max_val))
        plt.axes(xlim=(-0.01, 0.501))
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('Error [%vol]', fontproperties=prop)
        plt.title(titlestr(t), fontproperties=prop,fontsize=20)
        plt.plot(xcells, abs_error_calcite * 100 / (1 - phi), label='Calcite', **line_error('C0'))
        plt.plot(xcells, abs_error_dolomite * 100 / (1 - phi), label='Dolomite', **line_error('C1'))
        plt.legend(loc='center right', prop=prop)
        plt.savefig(folder_general + f'/error-abs-calcite-dolomite-{i}.pdf')
        plt.tight_layout()
        plt.close()

        min_val = 0.9*min([min(abs_error_cacation[0:int(ncells/2)]), min(abs_error_mgcation[0:int(ncells/2)]), min(abs_error_hco3anion[0:int(ncells/2)]), min(abs_error_co2aq[0:int(ncells/2)]), min(abs_error_hcation[0:int(ncells/2)])])
        max_val = 1.1*max([max(abs_error_cacation[0:int(ncells/2)]), max(abs_error_mgcation[0:int(ncells/2)]), max(abs_error_hco3anion[0:int(ncells/2)]), max(abs_error_co2aq[0:int(ncells/2)]), max(abs_error_hcation[0:int(ncells/2)])])
        plt.axes(xlim=(-0.01, 0.501))
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('Error [molal]', fontproperties=prop)
        plt.yscale('log')
        plt.title(titlestr(t), fontproperties=prop,fontsize=20)

        plt.plot(xcells, abs_error_cacation, label=r'$Ca^{2+}$', **line_error('C0'))
        plt.plot(xcells, abs_error_mgcation, label=r'$Mg^{2+}$', **line_error('C1'))
        plt.plot(xcells, abs_error_hco3anion, label=r'$HCO_3^-$', **line_error('C2'))
        plt.plot(xcells, abs_error_co2aq, label=r'$CO_2(aq)$', **line_error('C3'))
        plt.plot(xcells, abs_error_hcation, label=r'$H^+$', **line_error('C4'))
        plt.plot(xcells, abs_error_CO3anion, label=r'$CO_3^{2-}$', **line_error('C5'))[0],
        plt.plot(xcells, abs_error_CaHCO3cation, label=r'$Ca(HCO_3)^+$', **line_error('C6'))[0],
        plt.plot(xcells, abs_error_CaClcation, label=r'$CaCl^+$',**line_error('C7'))[0],
        plt.plot(xcells, abs_error_MgClcation, label=r'$MgCl^+$',**line_error('C8'))[0],
        plt.plot(xcells, abs_error_MgHCO3caanion, label=r'$Mg(HCO_3))^+$',**line_error('C9'))[0],
        plt.plot(xcells, abs_error_OHanion, label=r'$OH^-$', **line_error('darkviolet'))[0],

        plt.legend(loc='center right', prop=prop)
        plt.savefig(folder_general + f'/error-abs-aqueous-species-{i}.pdf')
        plt.tight_layout()
        plt.close()

        '''
        #
        # abs_error_cacation[abs_error_cacation == 0] = 1e-12
        # abs_error_mgcation[abs_error_mgcation == 0] = 1e-12
        # abs_error_hco3anion[abs_error_hco3anion == 0] = 1e-12
        # abs_error_co2aq[abs_error_co2aq == 0] = 1e-12
        # abs_error_hcation[abs_error_hcation == 0] = 1e-12
        #
        # plt.axes(xlim=(-0.01, 0.501))
        # plt.xlabel('Distance [m]')
        # plt.ylabel('Error [log 10 molal]', fontproperties=prop)
        # plt.title(titlestr(t), fontproperties=prop,fontsize=20)
        # plt.plot(xcells, np.log10(abs_error_cacation), label=r'$\mathrm{Ca^{2+}}$', **line_error('C0'))
        # plt.plot(xcells, np.log10(abs_error_mgcation), label=r'$\mathrm{Mg^{2+}}$', **line_error('C1'))
        # plt.plot(xcells, np.log10(abs_error_hco3anion), label=r'$\mathrm{HCO_3^{-}}$', **line_error('C2'))
        # plt.plot(xcells, np.log10(abs_error_co2aq), label=r'$\mathrm{CO_2(aq)}$', **line_error('red'))
        # plt.plot(xcells, np.log10(abs_error_hcation), label=r'$\mathrm{H^+}$', **line_error('darkviolet'))
        # plt.legend(loc='upper right', prop=prop)
        # plt.savefig(folder_general + f'/error-log10-abs-aqueous-species-{i}.pdf')
        # plt.tight_layout()
        # plt.close()

        '''
        for j in range(ncells):
            rel_error_calcite[j] = abs_error_calcite[j] / data_class_calcite[j] if data_class_calcite[j] > tol else 0.0
            rel_error_dolomite[j] = abs_error_dolomite[j] / data_class_dolomite[j] if data_class_dolomite[j] > tol else 0.0
            rel_error_cacation[j] = abs_error_cacation[j] / data_class_cacation[j] if data_class_cacation[j] > tol else 0.0
            rel_error_mgcation[j] = abs_error_mgcation[j] / data_class_mgcation[j] if data_class_mgcation[j] > tol else 0.0
            rel_error_hco3anion[j] = abs_error_hco3anion[j] / data_class_hco3anion[j] if data_class_hco3anion[j] > tol else 0.0
            rel_error_co2aq[j] = abs_error_co2aq[j] / data_class_co2aq[j] if data_class_co2aq[j] > tol else 0.0
            rel_error_hcation[j] = abs_error_hcation[j] / data_class_hcation[j] if data_class_hcation[j] > tol else 0.0


        plt.axes(xlim=(-0.01, 0.501))
        plt.xlabel('Distance [m]')
        plt.ylabel('Relative Error [-]')
        plt.title(titlestr(t))
        plt.plot(xcells, rel_error_calcite, label='Calcite', **line('C0'))
        plt.plot(xcells, rel_error_dolomite, label='Dolomite', **line('C1'))
        plt.legend(loc='center right')
        plt.savefig(folder_general + f'/error-rel-calcite-dolomite-short-{i}.pdf')
        plt.tight_layout()
        plt.close()

        plt.axes(xlim=(-0.01, 0.501))
        plt.xlabel('Distance [m]')
        plt.ylabel('Relative Error [-]')
        plt.yscale('log')
        plt.title(titlestr(t))
        plt.plot(xcells, rel_error_cacation, label=r'$\mathrm{Ca^{2+}}$', **line('C0'))
        plt.plot(xcells, rel_error_mgcation, label=r'$\mathrm{Mg^{2+}}$', **line('C1'))
        plt.plot(xcells, rel_error_hco3anion, label=r'$\mathrm{HCO_3^{-}}$', **line('C2'))
        plt.plot(xcells, rel_error_co2aq, label=r'$\mathrm{CO_2(aq)}$', **line('red'))
        plt.plot(xcells, rel_error_hcation, label=r'$\mathrm{H^+}$', **line('darkviolet'))
        plt.legend(loc='upper right')
        plt.savefig(folder_general + f'/error-rel-aqueous-species-short-{i}.pdf')
        plt.tight_layout()
        plt.close()
        '''
        '''
        # ------------------------------------------------------------------------------------------------------ #
        # L1-norm
        # ------------------------------------------------------------------------------------------------------ #

        # Calculate difference between reference and approximate data with l1 norm
        abs_error_l1_ph[i] = la.norm(data_class_ph - data_smart_ph, 1)
        abs_error_l1_calcite[i] = la.norm(data_class_calcite - data_smart_calcite, 1)
        abs_error_l1_dolomite[i] = la.norm(data_class_dolomite - data_smart_dolomite, 1)
        abs_error_l1_hcation[i] = la.norm(data_class_hcation - data_smart_hcation, 1)
        abs_error_l1_cacation[i] = la.norm(data_class_cacation - data_smart_cacation, 1)
        abs_error_l1_mgcation[i] = la.norm(data_class_mgcation - data_smart_mgcation, 1)
        abs_error_l1_hco3anion[i] = la.norm(data_class_hco3anion - data_smart_hco3anion, 1)
        abs_error_l1_cacation[i] = la.norm(data_class_cacation - data_smart_cacation, 1)
        abs_error_l1_co2aq[i] = la.norm(data_class_co2aq - data_smart_co2aq, 1)

        # Calculate with l1 norm of the reference data
        norm_l1_ph[i] = la.norm(data_class_ph, 1)
        norm_l1_calcite[i] = la.norm(data_class_calcite, 1)
        norm_l1_dolomite[i] = la.norm(data_class_dolomite, 1)
        norm_l1_hcation[i] = la.norm(data_class_hcation, 1)
        norm_l1_cacation[i] = la.norm(data_class_cacation, 1)
        norm_l1_mgcation[i] = la.norm(data_class_mgcation, 1)
        norm_l1_hco3anion[i] = la.norm(data_class_hco3anion, 1)
        norm_l1_cacation[i] = la.norm(data_class_cacation, 1)
        norm_l1_co2aq[i] = la.norm(data_class_co2aq, 1)

        # Calculate the error by |n_ref - n_approx|_l1 / |n_ref|_l1
        # if norm_l1_calcite[i] < tol:
        #    print(norm_l1_calcite[i])
        #    input()
        # if i >= 8000:
        #    print(f"|n_ref - n_approx|_l1 = {abs_error_l1_calcite[i]:>15} \t |n_ref|_l1 = {norm_l1_calcite[i]}")
        #    input()

        error_l1_ph[i] = abs_error_l1_ph[i] / norm_l1_ph[i] if norm_l1_ph[i] > tol else 0.0
        error_l1_calcite[i] = abs_error_l1_calcite[i] / norm_l1_calcite[i] if norm_l1_calcite[i] > tol else 0.0
        error_l1_dolomite[i] = abs_error_l1_dolomite[i] / norm_l1_dolomite[i] if norm_l1_dolomite[i] > tol else 0.0
        error_l1_hcation[i] = abs_error_l1_hcation[i] / norm_l1_hcation[i] if norm_l1_hcation[i] > tol else 0.0
        error_l1_cacation[i] = abs_error_l1_cacation[i] / norm_l1_cacation[i] if norm_l1_cacation[i] > tol else 0.0
        error_l1_mgcation[i] = abs_error_l1_mgcation[i] / norm_l1_mgcation[i] if norm_l1_mgcation[i] > tol else 0.0
        error_l1_hco3anion[i] = abs_error_l1_hco3anion[i] / norm_l1_hco3anion[i] if norm_l1_hco3anion[i] > tol else 0.0
        error_l1_cacation[i] = abs_error_l1_cacation[i] / norm_l1_cacation[i] if norm_l1_cacation[i] > tol else 0.0
        error_l1_co2aq[i] = abs_error_l1_co2aq[i] / norm_l1_co2aq[i] if norm_l1_co2aq[i] > tol else 0.0

        # ------------------------------------------------------------------------------------------------------ #
        # L2-norm
        # ------------------------------------------------------------------------------------------------------ #

        abs_error_rms_ph[i] = np.sqrt(np.mean((data_class_ph - data_smart_ph) ** 2))
        abs_error_rms_calcite[i] = np.sqrt(np.mean((data_class_calcite - data_smart_calcite) ** 2))
        abs_error_rms_dolomite[i] = np.sqrt(np.mean((data_class_dolomite - data_smart_dolomite) ** 2))
        abs_error_rms_hcation[i] = np.sqrt(np.mean((data_class_hcation - data_smart_hcation) ** 2))
        abs_error_rms_cacation[i] = np.sqrt(np.mean((data_class_cacation - data_smart_cacation) ** 2))
        abs_error_rms_mgcation[i] = np.sqrt(np.mean((data_class_mgcation - data_smart_mgcation) ** 2))
        abs_error_rms_hco3anion[i] = np.sqrt(np.mean((data_class_hco3anion - data_smart_hco3anion) ** 2))
        abs_error_rms_cacation[i] = np.sqrt(np.mean((data_class_cacation - data_smart_cacation) ** 2))
        abs_error_rms_co2aq[i] = np.sqrt(np.mean((data_class_co2aq - data_smart_co2aq) ** 2))

        norm_rms_ph[i] = np.sqrt(np.mean(data_class_ph ** 2))
        norm_rms_calcite[i] = np.sqrt(np.mean(data_class_calcite ** 2))
        norm_rms_dolomite[i] = np.sqrt(np.mean(data_class_dolomite ** 2))
        norm_rms_hcation[i] = np.sqrt(np.mean(data_class_hcation ** 2))
        norm_rms_cacation[i] = np.sqrt(np.mean(data_class_cacation ** 2))
        norm_rms_mgcation[i] = np.sqrt(np.mean(data_class_mgcation ** 2))
        norm_rms_hco3anion[i] = np.sqrt(np.mean(data_class_hco3anion ** 2))
        norm_rms_cacation[i] = np.sqrt(np.mean(data_class_cacation ** 2))
        norm_rms_co2aq[i] = np.sqrt(np.mean(data_class_co2aq ** 2))

        error_rms_ph[i] = abs_error_rms_ph[i] / norm_rms_ph[i] if norm_rms_ph[i] > tol else 0.0
        error_rms_calcite[i] = abs_error_rms_calcite[i] / norm_rms_calcite[i] if norm_rms_calcite[i] > tol else 0.0
        error_rms_dolomite[i] = abs_error_rms_dolomite[i] / norm_rms_dolomite[i] if norm_rms_dolomite[i] > tol else 0.0
        error_rms_hcation[i] = abs_error_rms_hcation[i] / norm_rms_hcation[i] if norm_rms_hcation[i] > tol else 0.0
        error_rms_cacation[i] = abs_error_rms_cacation[i] / norm_rms_cacation[i] if norm_rms_cacation[i] > tol else 0.0
        error_rms_mgcation[i] = abs_error_rms_mgcation[i] / norm_rms_mgcation[i] if norm_rms_mgcation[i] > tol else 0.0
        error_rms_hco3anion[i] = abs_error_rms_hco3anion[i] / norm_rms_hco3anion[i] if norm_rms_hco3anion[
                                                                                           i] > tol else 0.0
        error_rms_cacation[i] = abs_error_rms_cacation[i] / norm_rms_cacation[i] if norm_rms_cacation[i] > tol else 0.0
        error_rms_co2aq[i] = abs_error_rms_co2aq[i] / norm_rms_co2aq[i] if norm_rms_co2aq[i] > tol else 0.0

        # ------------------------------------------------------------------------------------------------------ #
        # Inf-norm
        # ------------------------------------------------------------------------------------------------------ #

        # Calculate difference between reference and approximate data with l1 norm
        abs_error_inf_ph[i] = la.norm(data_class_ph - data_smart_ph, np.inf)
        abs_error_inf_calcite[i] = la.norm(data_class_calcite - data_smart_calcite, np.inf)
        abs_error_inf_dolomite[i] = la.norm(data_class_dolomite - data_smart_dolomite, np.inf)
        abs_error_inf_hcation[i] = la.norm(data_class_hcation - data_smart_hcation, np.inf)
        abs_error_inf_cacation[i] = la.norm(data_class_cacation - data_smart_cacation, np.inf)
        abs_error_inf_mgcation[i] = la.norm(data_class_mgcation - data_smart_mgcation, np.inf)
        abs_error_inf_hco3anion[i] = la.norm(data_class_hco3anion - data_smart_hco3anion, np.inf)
        abs_error_inf_cacation[i] = la.norm(data_class_cacation - data_smart_cacation, np.inf)
        abs_error_inf_co2aq[i] = la.norm(data_class_co2aq - data_smart_co2aq, np.inf)

        # Calculate with l1 norm of the reference data
        norm_inf_ph[i] = la.norm(data_class_ph, np.inf)
        norm_inf_calcite[i] = la.norm(data_class_calcite, np.inf)
        norm_inf_dolomite[i] = la.norm(data_class_dolomite, np.inf)
        norm_inf_hcation[i] = la.norm(data_class_hcation, np.inf)
        norm_inf_cacation[i] = la.norm(data_class_cacation, np.inf)
        norm_inf_mgcation[i] = la.norm(data_class_mgcation, np.inf)
        norm_inf_hco3anion[i] = la.norm(data_class_hco3anion, np.inf)
        norm_inf_cacation[i] = la.norm(data_class_cacation, np.inf)
        norm_inf_co2aq[i] = la.norm(data_class_co2aq, np.inf)

        # Calculate the error by |n_ref - n_approx|_inf / |n_ref|_inf
        #if norm_inf_calcite[i] < tol:
        #    print(norm_inf_calcite[i])
        #    input()
        # if i >= 8000:
        #    print(f"|n_ref - n_approx|_inf = {abs_error_inf_calcite[i]:>15} \t |n_ref|_inf = {norm_inf_calcite[i]}")
        #    input()

        error_inf_ph[i] = abs_error_inf_ph[i] / norm_inf_ph[i] if norm_inf_ph[i] > tol else 0.0
        error_inf_calcite[i] = abs_error_inf_calcite[i] / norm_inf_calcite[i] if norm_inf_calcite[i] > tol else 0.0
        error_inf_dolomite[i] = abs_error_inf_dolomite[i] / norm_inf_dolomite[i] if norm_inf_dolomite[i] > tol else 0.0
        error_inf_hcation[i] = abs_error_inf_hcation[i] / norm_inf_hcation[i] if norm_inf_hcation[i] > tol else 0.0
        error_inf_cacation[i] = abs_error_inf_cacation[i] / norm_inf_cacation[i] if norm_inf_cacation[i] > tol else 0.0
        error_inf_mgcation[i] = abs_error_inf_mgcation[i] / norm_inf_mgcation[i] if norm_inf_mgcation[i] > tol else 0.0
        error_inf_hco3anion[i] = abs_error_inf_hco3anion[i] / norm_inf_hco3anion[i] if norm_inf_hco3anion[
                                                                                           i] > tol else 0.0
        error_inf_cacation[i] = abs_error_inf_cacation[i] / norm_inf_cacation[i] if norm_inf_cacation[i] > tol else 0.0
        error_inf_co2aq[i] = abs_error_inf_co2aq[i] / norm_inf_co2aq[i] if norm_inf_co2aq[i] > tol else 0.0
        '''
    # -----------------------------------------------------------------------------------------------------------------#
    # Plot l1-error
    # -----------------------------------------------------------------------------------------------------------------#
    '''
    step = 40
    plt.xlabel('Time Step')
    #plt.axes(ylim=(pow(10, -6), pow(10, -1)))
    #plt.yscale('log')
    plt.ylabel(r'Error in $\ell_1$-norm')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], error_l1_calcite[0:nsteps:step], label=r'$\| e_{\rm CaCO_3}\|_{\ell_1}$',
             color='C0', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_l1_dolomite[0:nsteps:step], label=r'$\| e_{\rm CaMg(CO_3)_2}\|_{\ell_1}$',
             color='C1', linewidth=1)
    leg = plt.legend(loc='lower right')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/errors-calcite-dolomite-l1.pdf')
    plt.close()

    plt.xlabel('Time Step')
    #plt.axes(ylim=(pow(10, -5), pow(10, 0)))
    #plt.yscale('log')
    plt.ylabel(r'Error in $\ell_1$-norm')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], error_l1_cacation[0:nsteps:step], label=r'$\| e_{\rm Ca^{2+}}\|_{\ell_1}$',
             color='C0', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_l1_mgcation[0:nsteps:step], label=r'$\| e_{\rm Mg^{2+}}\|_{\ell_1}$',
             color='C1', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_l1_hco3anion[0:nsteps:step], label=r'$\| e_{\rm HCO_3^{-}}\|_{\ell_1}$',
             color='C2', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_l1_co2aq[0:nsteps:step], label=r'$\| e_{\rm CO_2(aq)}\|_{\ell_1}$',
             color='red', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_l1_hcation[0:nsteps:step], label=r'$\| e_{\rm H^+}\|_{\ell_1}$',
             color='darkviolet', linewidth=1)
    leg = plt.legend(loc='lower right')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/errors-aqueous-l1.pdf')
    plt.close()

    step = 40
    plt.xlabel('Time Step')
    plt.axes(ylim=(pow(10, -6), pow(10, -2)))
    plt.yscale('log')
    plt.ylabel(r'Error in $\ell_1$-norm, $\log$-scale')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], error_l1_calcite[0:nsteps:step], label=r'$\| e_{\rm CaCO_3}\|_{\ell_1}$',
             color='C0', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_l1_dolomite[0:nsteps:step], label=r'$\| e_{\rm CaMg(CO_3)_2}\|_{\ell_1}$',
             color='C1', linewidth=1)
    leg = plt.legend(loc='lower right')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/errors-calcite-dolomite-l1-log.pdf')
    plt.close()

    plt.xlabel('Time Step')
    plt.axes(ylim=(pow(10, -5), pow(10, 0)))
    plt.yscale('log')
    plt.ylabel(r'Error in $\ell_1$-norm, $\log$-scale')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], error_l1_cacation[0:nsteps:step], label=r'$\| e_{\rm Ca^{2+}}\|_{\ell_1}$',
             color='C0', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_l1_mgcation[0:nsteps:step], label=r'$\| e_{\rm Mg^{2+}}\|_{\ell_1}$',
             color='C1', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_l1_hco3anion[0:nsteps:step], label=r'$\| e_{\rm HCO_3^{-}}\|_{\ell_1}$',
             color='C2', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_l1_co2aq[0:nsteps:step], label=r'$\| e_{\rm CO_2(aq)}\|_{\ell_1}$',
             color='red', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_l1_hcation[0:nsteps:step], label=r'$\| e_{\rm H^+}\|_{\ell_1}$',
             color='darkviolet', linewidth=1)
    leg = plt.legend(loc='lower right')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/errors-aqueous-l1-log.pdf')
    plt.close()


    # -----------------------------------------------------------------------------------------------------------------#
    # Plot l1-error in %
    # -----------------------------------------------------------------------------------------------------------------#

    plt.xlabel('Time Step')
    plt.ylabel(r'Error in $\ell_1$-norm [%]')
    #plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], 100 * error_l1_calcite[0:nsteps:step], label=r'$\| e_{\rm CaCO_3}\|$',
             color='C0', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], 100 * error_l1_dolomite[0:nsteps:step], label=r'$\| e_{\rm CaMg(CO_3)_2}\|$',
             color='C1', linewidth=1)
    leg = plt.legend(loc='lower right')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/errors-calcite-dolomite-l1-percent.pdf')
    plt.close()

    plt.xlabel('Time Step')
    plt.ylabel(r'Error in $\ell_1$-norm [%]')
    #plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], 100 * error_l1_cacation[0:nsteps:step],
             label=r'$\| e_{\rm Ca^{2+}}\|$', color='C0', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], 100 * error_l1_mgcation[0:nsteps:step],
             label=r'$\| e_{\rm Mg^{2+}}\|$', color='C1', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], 100 * error_l1_hco3anion[0:nsteps:step],
             label=r'$\| e_{\rm HCO_3^{-}}\|$', color='C2', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], 100 * error_l1_co2aq[0:nsteps:step], label=r'$\| e_{\rm CO_2(aq)}\|$',
             color='red', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], 100 * error_l1_hcation[0:nsteps:step], label=r'$\| e_{\rm H^+}\|$',
             color='darkviolet', linewidth=1)
    leg = plt.legend(loc='upper left')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/errors-aqueous-l1-percent.pdf')
    plt.close()
    '''
    '''
    # -----------------------------------------------------------------------------------------------------------------#
    # Plot RMS-error
    # -----------------------------------------------------------------------------------------------------------------#

    plt.xlabel('Time Step')
    plt.ylabel(r'RMS Error')
    plt.axes(ylim=(pow(10, -6), pow(10, -1)))
    plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], error_rms_calcite[0:nsteps:step], label=r'$\| e_{\rm Calcite}\|_{\rm RMS}$',
             color='C0', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_rms_dolomite[0:nsteps:step], label=r'$\| e_{\rm Dolomite}\|_{\rm RMS}$',
             color='C1', linewidth=1)
    leg = plt.legend(loc='lower right', bbox_to_anchor=(1, 0.13))
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/errors-calcite-dolomite-rms.pdf')
    plt.close()

    plt.xlabel('Time Step')
    plt.ylabel('RMS Error')
    plt.axes(ylim=(pow(10, -5), pow(10, 0)))
    plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], error_rms_cacation[0:nsteps:step], label=r'$\| e_{\rm Ca^{2+}}\|_{\rm RMS}$',
             color='C0', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_rms_mgcation[0:nsteps:step], label=r'$\| e_{\rm Mg^{2+}}\|_{\rm RMS}$',
             color='C1', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_rms_hco3anion[0:nsteps:step], label=r'$\| e_{\rm HCO_3^{-}}\|_{\rm RMS}$',
             color='C2', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_rms_co2aq[0:nsteps:step], label=r'$\| e_{\rm CO_2(aq)}\|_{\rm RMS}$',
             color='red', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_rms_hcation[0:nsteps:step], label=r'$\| e_{\rm H^+}\|_{\rm RMS}$',
             color='darkviolet', linewidth=1)
    leg = plt.legend(loc='lower right')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/errors-aqueous-rms.pdf')
    plt.close()

    # -----------------------------------------------------------------------------------------------------------------#
    # Plot inf-error
    # -----------------------------------------------------------------------------------------------------------------#

    step = 40
    plt.xlabel('Time Step')
    plt.axes(ylim=(pow(10, -6), pow(10, 0)))
    plt.yscale('log')
    plt.ylabel(r'Error in $\infty$-norm')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], error_inf_calcite[0:nsteps:step],
             label=r'$\| e_{\rm CaCO_3}\|_{\infty}$',
             color='C0', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_inf_dolomite[0:nsteps:step],
             label=r'$\| e_{\rm CaMg(CO_3)_2}\|_{\infty}$',
             color='C1', linewidth=1)
    leg = plt.legend(loc='lower right')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig('errors-calcite-dolomite-inf.pdf')
    plt.close()

    plt.xlabel('Time Step')
    plt.axes(ylim=(pow(10, -5), pow(10, 0)))
    plt.yscale('log')
    plt.ylabel(r'Error in $\infty$-norm')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], error_inf_cacation[0:nsteps:step], label=r'$\| e_{\rm Ca^{2+}}\|_{\infty}$',
             color='C0', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_inf_mgcation[0:nsteps:step], label=r'$\| e_{\rm Mg^{2+}}\|_{\infty}$',
             color='C1', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_inf_hco3anion[0:nsteps:step], label=r'$\| e_{\rm HCO_3^{-}}\|_{\infty}$',
             color='C2', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_inf_co2aq[0:nsteps:step], label=r'$\| e_{\rm CO_2(aq)}\|_{\infty}$',
             color='red', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], error_inf_hcation[0:nsteps:step], label=r'$\| e_{\rm H^+}\|_{\infty}$',
             color='darkviolet', linewidth=1)
    leg = plt.legend(loc='lower right')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/errors-aqueous-inf.pdf')
    plt.close()

    # -----------------------------------------------------------------------------------------------------------------#
    # Plot inf-error in %
    # -----------------------------------------------------------------------------------------------------------------#

    plt.xlabel('Time Step')
    plt.ylabel(r'Error in $\infty$-norm [%]')
    plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], 100 * error_inf_calcite[0:nsteps:step], label=r'$\| e_{\rm CaCO_3}\|$',
             color='C0', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], 100 * error_inf_dolomite[0:nsteps:step], label=r'$\| e_{\rm CaMg(CO_3)_2}\|$',
             color='C1', linewidth=1)
    leg = plt.legend(loc='lower right')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig('errors-calcite-dolomite-inf-percent.pdf')
    plt.close()

    plt.xlabel('Time Step')
    plt.ylabel(r'Error in $\infty$-norm [%]')
    plt.yscale('log')
    plt.ticklabel_format(style='plain', axis='x')
    plt.plot(time_steps[0:nsteps:step], 100 * error_inf_cacation[0:nsteps:step],
             label=r'$\| e_{\rm Ca^{2+}}\|$', color='C0', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], 100 * error_inf_mgcation[0:nsteps:step],
             label=r'$\| e_{\rm Mg^{2+}}\|$', color='C1', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], 100 * error_inf_hco3anion[0:nsteps:step],
             label=r'$\| e_{\rm HCO_3^{-}}\|$', color='C2', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], 100 * error_inf_co2aq[0:nsteps:step], label=r'$\| e_{\rm CO_2(aq)}\|$',
             color='red', linewidth=1)
    plt.plot(time_steps[0:nsteps:step], 100 * error_inf_hcation[0:nsteps:step], label=r'$\| e_{\rm H^+}\|$',
             color='darkviolet', linewidth=1)
    leg = plt.legend(loc='upper left')
    for line in leg.get_lines(): line.set_linewidth(2.0)
    plt.tight_layout()
    plt.savefig(folder_general + '/errors-aqueous-inf-percent.pdf')
    plt.close()

    print('max |e in pH|_l1 = %2.2e' % np.max(error_l1_ph))
    print('max |e in calcite|_l1 = %2.2e' % np.max(error_l1_calcite))
    print('max |e in dolomite|_l1 = %2.2e' % np.max(error_l1_dolomite))
    print('max |e in hcation|_l1 = %2.2e' % np.max(error_l1_hcation))
    print('max |e in cacation|_l1 = %2.2e' % np.max(error_l1_cacation))
    print('max |e in mgcation|_l1 = %2.2e' % np.max(error_l1_mgcation))
    print('max |e in hco3anion|_l1 = %2.2e' % np.max(error_l1_hco3anion))
    print('max |e in cacation|_l1 = %2.2e' % np.max(error_l1_cacation))
    print('max |e in co2aq|_l1 = %2.2e\n' % np.max(error_l1_co2aq))

    print('max |e in pH|_l1 = %2.2f percent' % (100 * np.max(error_l1_ph)))
    print('max |e in calcite|_l1 = %2.2f percent' % (100 * np.max(error_l1_calcite)))
    print('max |e in dolomite|_l1 = %2.2f percent' % (100 * np.max(error_l1_dolomite)))
    print('max |e in hcation|_l1 = %2.2f percent' % (100 * np.max(error_l1_hcation)))
    print('max |e in cacation|_l1 = %2.2f percent' % (100 * np.max(error_l1_cacation)))
    print('max |e in mgcation|_l1 = %2.2f percent' % (100 * np.max(error_l1_mgcation)))
    print('max |e in hco3anion|_l1 = %2.2f percent' % (100 * np.max(error_l1_hco3anion)))
    print('max |e in cacation|_l1 = %2.2f percent' % (100 * np.max(error_l1_cacation)))
    print('max |e in co2aq|_l1 = %2.2f percent\n' % (100 * np.max(error_l1_co2aq)))

    print('max |e in pH|_rms = %2.2e' % np.max(error_rms_ph))
    print('max |e in calcite|_rms = %2.2e' % np.max(error_rms_calcite))
    print('max |e in dolomite|_rms = %2.2e' % np.max(error_rms_dolomite))
    print('max |e in hcation|_rms = %2.2e' % np.max(error_rms_hcation))
    print('max |e in cacation|_rms = %2.2e' % np.max(error_rms_cacation))
    print('max |e in mgcation|_rms = %2.2e' % np.max(error_rms_mgcation))
    print('max |e in hco3anion|_rms = %2.2e' % np.max(error_rms_hco3anion))
    print('max |e in cacation|_rms = %2.2e' % np.max(error_rms_cacation))
    print('max |e in co2aq|_rms = %2.2e\n' % np.max(error_rms_co2aq))

    print('max |e in pH|_rms = %2.2f percent' % (100 * np.max(error_rms_ph)))
    print('max |e in calcite|_rms = %2.2f percent' % (100 * np.max(error_rms_calcite)))
    print('max |e in dolomite|_rms = %2.2f percent' % (100 * np.max(error_rms_dolomite)))
    print('max |e in hcation|_rms = %2.2f percent' % (100 * np.max(error_rms_hcation)))
    print('max |e in cacation|_rms = %2.2f percent' % (100 * np.max(error_rms_cacation)))
    print('max |e in mgcation|_rms = %2.2f percent' % (100 * np.max(error_rms_mgcation)))
    print('max |e in hco3anion|_rms = %2.2f percent' % (100 * np.max(error_rms_hco3anion)))
    print('max |e in cacation|_rms = %2.2f percent' % (100 * np.max(error_rms_cacation)))
    print('max |e in co2aq|_rms = %2.2f percent' % (100 * np.max(error_rms_co2aq)))

    print('max |e in pH|_inf = %2.2e' % np.max(error_inf_ph))
    print('max |e in calcite|_inf = %2.2e' % np.max(error_inf_calcite))
    print('max |e in dolomite|_inf = %2.2e' % np.max(error_inf_dolomite))
    print('max |e in hcation|_inf = %2.2e' % np.max(error_inf_hcation))
    print('max |e in cacation|_inf = %2.2e' % np.max(error_inf_cacation))
    print('max |e in mgcation|_inf = %2.2e' % np.max(error_inf_mgcation))
    print('max |e in hco3anion|_inf = %2.2e' % np.max(error_inf_hco3anion))
    print('max |e in cacation|_inf = %2.2e' % np.max(error_inf_cacation))
    print('max |e in co2aq|_inf = %2.2e\n' % np.max(error_inf_co2aq))

    print('max |e in pH|_inf = %2.2f in percent' % (100 * np.max(error_inf_ph)))
    print('max |e in calcite|_inf = %2.2f in percent' % (100 * np.max(error_inf_calcite)))
    print('max |e in dolomite|_inf = %2.2f in percent' % (100 * np.max(error_inf_dolomite)))
    print('max |e in hcation|_inf = %2.2f in percent' % (100 * np.max(error_inf_hcation)))
    print('max |e in cacation|_inf = %2.2f in percent' % (100 * np.max(error_inf_cacation)))
    print('max |e in mgcation|_inf = %2.2f in percent' % (100 * np.max(error_inf_mgcation)))
    print('max |e in hco3anion|_inf = %2.2f in percent' % (100 * np.max(error_inf_hco3anion)))
    print('max |e in cacation|_inf = %2.2f in percent' % (100 * np.max(error_inf_cacation)))
    print('max |e in co2aq|_inf = %2.2f in percent' % (100 * np.max(error_inf_co2aq)))
    '''

def mass_conservation():

    indices = [1, 20, 2400]
    #indices = [20]
    for i in indices:

        filearray_res = np.loadtxt(folder_smart + f'/res-{i-1}.txt')
        filearray_state = np.loadtxt(folder_smart + f'/test-{i-1}.txt', skiprows=1)

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

        #plt.axes(xlim=(-0.01, 0.501), ylim=(1e-16, 5e-9))
        plt.axes(xlim=(-0.01, 0.501))
        plt.xlabel('Distance [m]', fontproperties=prop)
        plt.ylabel('Mass Conservation Relative Error [-]', fontproperties=prop)
        plt.yscale('log')
        plt.title(titlestr(dt * (i+1)), fontproperties=prop,fontsize=20)
        plt.plot(xcells, error_C, label='C', **line_error('C0'))
        plt.plot(xcells, error_Ca, label='Ca', **line_error('C1'))
        plt.plot(xcells, error_Cl, label='Cl', **line_error('C2'))
        plt.plot(xcells, error_H, label='H', **line_error('C3'))
        plt.plot(xcells, error_Mg, label='Mg', **line_error('C4'))
        plt.plot(xcells, error_Na, label='Na', **line_error('C5'))
        plt.plot(xcells, error_O, label='O', **line_error('C6'))
        plt.plot(xcells, error_Si, label='Si', **line_error('C7'))

        plt.legend(loc='upper right', prop=prop)
        #if i == 20: plt.legend(loc='upper left', prop=prop)
        #else: plt.legend(loc='upper right', prop=prop)
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
    files_smart = [file for file in natsorted( os.listdir(folder_smart) ) if (("res" not in file) and ("analysis" not in file))]
    files_class = [file for file in natsorted( os.listdir(folder_class) ) if ("analysis" not in file)]

    plot_on_demand_learning_countings()
    # plot_on_demand_learnings_total()
    # plot_speedups()
    # plot_computing_costs_vs_total_learnings()
    # plot_on_demand_learnings_total_initial_steps(10)
    # plot_computing_costs()
    # #plot_figures_ph()
    # plot_figures_aqueous_species()
    # plot_figures_calcite_dolomite()
    # #plot_figures_calcite_dolomite_moles()
    #
    # tolerance = 1e-2
    # calculate_error(tolerance)
    # mass_conservation()
