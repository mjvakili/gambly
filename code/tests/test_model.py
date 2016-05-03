import numpy as np
import matplotlib.pyplot as plt
from halotools.sim_manager import CachedHaloCatalog
from halotools.empirical_models import HodModelFactory
from halotools.empirical_models import TrivialPhaseSpace, Zheng07Cens
from halotools.empirical_models import Zheng07Sats , HeavisideAssembias
from halotools.empirical_models import PrebuiltHodModelFactory
from halotools.mock_observables import tpcf
from halotools.empirical_models.factories.mock_helpers import three_dim_pos_bundle
from halotools.mock_observables import FoFGroups
from halotools.mock_observables.pair_counters import npairs_3d
from halotools.mock_observables.catalog_analysis_helpers import return_xyz_formatted_array
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors
from Corrfunc.utils import read_catalog
from halotools.empirical_models import NFWPhaseSpace
from matplotlib import lines as mlines

halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')

class AssembiasZheng07Sats(Zheng07Sats, HeavisideAssembias):

    def __init__(self, **kwargs):

        Zheng07Sats.__init__(self, threshold = -21)

        HeavisideAssembias.__init__(self,
            method_name_to_decorate = 'mean_occupation',
            lower_assembias_bound = 0.,
            upper_assembias_bound = np.inf,
            **kwargs)

cens_occ_model =  Zheng07Cens(threshold = -21)
cens_prof_model = TrivialPhaseSpace()
sats_occ_model =  AssembiasZheng07Sats()
sats_prof_model = NFWPhaseSpace()

model= HodModelFactory(
        centrals_occupation = cens_occ_model,
        centrals_profile = cens_prof_model,
        satellites_occupation = sats_occ_model,
        satellites_profile = sats_prof_model)
baseline_model = PrebuiltHodModelFactory("Zheng07" , threshold = -21)


############### Setting the model parameters to those of Guo 15#######
model.param_dict['logM0'] =  12.59
model.param_dict['sigma_logM'] =  0.49
model.param_dict['logMmin'] =  12.78
model.param_dict['alpha'] =  1.14
model.param_dict['logM1'] =  13.99
baseline_model.param_dict['logM0'] =  12.59
baseline_model.param_dict['sigma_logM'] =  0.49
baseline_model.param_dict['logMmin'] =  12.78
baseline_model.param_dict['alpha'] =  1.14
baseline_model.param_dict['logM1'] =  13.99

############### Iniitializing the mass range##############

npts = 1e3
mass = np.logspace(11, 14, npts)

############## Maximum allowable assembly bias ############

model.param_dict['mean_occupation_satellites_assembias_param1'] = 1.0
nsat_a1_old = model.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=0)
nsat_a1_young = model.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=1)
############# half of maximum allowable assembly bias ############

model.param_dict['mean_occupation_satellites_assembias_param1'] = 0.5
nsat_a05_old = model.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=0)
nsat_a05_young = model.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=1)

############# Model with no assemly bias ########################
nsat_baseline = baseline_model.mean_occupation_satellites(prim_haloprop = mass)

fig2_sats = plt.figure()
ax = fig2_sats.add_subplot(111)

xlabel = ax.set_xlabel(r'$M_{\rm vir} [M_{\odot}]$', fontsize=25)
ylabel = ax.set_ylabel(r'$\langle N_{\rm s}\rangle$', fontsize=25)
title = ax.set_title('Satellite Galaxies with positive assembly bias', fontsize=20)

ax.plot(mass, nsat_a1_young, color='blue', linewidth=3.5)
ax.plot(mass, nsat_a1_old, color='red', linewidth=3.5)
ax.plot(mass, nsat_a05_old, color='r', linewidth=1.5)
ax.plot(mass, nsat_a05_young, color='b', linewidth=1.5)
ax.plot(mass, nsat_baseline, '--', color='k', linewidth=2.5)
ax.plot(mass, .5*(nsat_a1_young + nsat_a1_old) , '-' , color = "green", linewidth =.5)
ax.plot(mass, .5*(nsat_a05_young + nsat_a05_old) , '-' , color = "magenta" , linewidth =.5)

plt.loglog()
plt.xlim(xmin=1e11, xmax=1e14)
plt.ylim(ymin=5e-3, ymax=100)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

blue_line = mlines.Line2D([], [], ls = '-', c = 'b', linewidth=3, label = 'high-concentration halos')
red_line = mlines.Line2D([], [], ls = '-', c = 'r', linewidth=3, label = 'low-concentration halos')
first_legend = plt.legend(handles=[blue_line, red_line], frameon=False, loc='best', fontsize=17)

fig2_sats.savefig('nsat.pdf', bbox_extra_artists=[xlabel, ylabel], bbox_inches='tight')


############## Minimum allowable assembly bias ############

model.param_dict['mean_occupation_satellites_assembias_param1'] = -1.
nsat_a1_old = model.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=0)
nsat_a1_young = model.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=1)
model.populate_mock(halocat , enforce_PBC = True)
model.param_dict['mean_occupation_satellites_assembias_param1'] = -0.5
nsat_a05_old = model.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=0)
nsat_a05_young = model.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=1)

fig2_sats = plt.figure()
ax = fig2_sats.add_subplot(111)

xlabel = ax.set_xlabel(r'$M_{\rm vir} [M_{\odot}]$', fontsize=25)
ylabel = ax.set_ylabel(r'$\langle N_{\rm s}\rangle$', fontsize=25)
title = ax.set_title('Satellite Galaxies with negative assembly bias', fontsize=20)

ax.plot(mass, nsat_a1_young, color='blue', linewidth=3.5)
ax.plot(mass, nsat_a1_old, color='red', linewidth=3.5)
ax.plot(mass, nsat_a05_old, color='r', linewidth=1.5)
ax.plot(mass, nsat_a05_young, color='b', linewidth=1.5)
ax.plot(mass, nsat_baseline, '--', color='k', linewidth=2.5)
ax.plot(mass, .5*(nsat_a1_young + nsat_a1_old) , '-' , color = "green", linewidth =.5)
ax.plot(mass, .5*(nsat_a05_young + nsat_a05_old) , '-' , color = "magenta" , linewidth =.5)

plt.loglog()
plt.xlim(xmin=1e11, xmax=1e14)
plt.ylim(ymin=5e-3, ymax=100)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

blue_line = mlines.Line2D([], [], ls = '-', c = 'b', linewidth=3, label = 'high-concentration halos')
red_line = mlines.Line2D([], [], ls = '-', c = 'r', linewidth=3, label = 'low-concentration halos')
first_legend = plt.legend(handles=[blue_line, red_line], frameon=False, loc='best', fontsize=17)

fig2_sats.savefig('nsat2.pdf', bbox_extra_artists=[xlabel, ylabel], bbox_inches='tight')
