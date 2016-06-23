'''
Plotting modules
'''
from halotools.sim_manager import CachedHaloCatalog
import os.path as path

from Corrfunc import _countpairs
from Corrfunc.utils import read_catalog
from halotools.empirical_models.factories.mock_helpers import three_dim_pos_bundle
from halotools.mock_observables.catalog_analysis_helpers import return_xyz_formatted_array
from halotools.empirical_models import enforce_periodicity_of_box
from halotools.empirical_models import HodModelFactory
from halotools.empirical_models import TrivialPhaseSpace, AssembiasZheng07Cens
from halotools.empirical_models import AssembiasZheng07Sats , HeavisideAssembias
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

from matplotlib import lines as mlines
from matplotlib import gridspec
halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')
import os
import os.path as path 
import util as ut
import h5py
import corner
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import colorConverter
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors
plt.switch_backend("Agg")

from hod_group import MCMC_model as hod_gmf 
from hod import MCMC_model as hod_wp
from biased_hod_group import MCMC_model as dec_gmf
from biased_hod import MCMC_model as dec_wp 

from biased_hod import composite_model as dec_model
from hod import single_model as hod_model


def occupation_predictions(filename, Mr, nburnins, nchains , obs = "wp", model = "dec"):

    npts = 2e3
    mass = np.logspace(10, 14, npts)
    
    if model == "dec":
          mod = dec_model(Mr)
    if model == "hod":
          mod = hod_model(Mr)
    
    sample = h5py.File(filename , "r")["mcmc"][nburnins:nchains]
    sample = sample.reshape(sample.shape[0]*sample.shape[1] , sample.shape[2]) 
    
    if model == "dec":
        
    	ncen_old = []
    	ncen_young = []
        #ncen_all = []
        nsat_old = []
   	nsat_young = []
        #nsat_all = []
        for i in xrange(len(sample)):
         
          print i
       
          mod.param_dict['logM0'] =  sample[i][0]
          mod.param_dict['sigma_logM'] =  sample[i][1]
          mod.param_dict['logMmin'] =  sample[i][2]
	  mod.param_dict['alpha'] =  sample[i][3]
	  mod.param_dict['logM1'] =  sample[i][4]
          mod.param_dict['mean_occupation_centrals_assembias_param1'] = sample[i][5]
          mod.param_dict['mean_occupation_satellites_assembias_param1'] = sample[i][6]
          
          ncen_old.append(mod.mean_occupation_centrals(prim_haloprop = mass, sec_haloprop_percentile=0))
          ncen_young.append(mod.mean_occupation_centrals(prim_haloprop = mass, sec_haloprop_percentile=1))
          #ncen_all.append(mod.mean_occupation_centrals(prim_haloprop = mass))
          nsat_old.append(mod.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=0))
          nsat_young.append(mod.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=1))
          #nsat_all.append(mod.mean_occupation_satellites(prim_haloprop = mass))

        #np.savetxt("results/ncen_all_"+obs+"_"+model+"_"+str(Mr)+".dat" , np.array(ncen_all)) 
        np.savetxt("results/ncen_old_"+obs+"_"+model+"_"+str(Mr)+".dat" , np.array(ncen_old)) 
        np.savetxt("results/ncen_young_"+obs+"_"+model+"_"+str(Mr)+".dat" , np.array(ncen_young)) 
        #np.savetxt("results/nsat_all_"+obs+"_"+model+"_"+str(Mr)+".dat" , np.array(nsat_all)) 
        np.savetxt("results/nsat_old_"+obs+"_"+model+"_"+str(Mr)+".dat" , np.array(nsat_old)) 
        np.savetxt("results/nsat_young_"+obs+"_"+model+"_"+str(Mr)+".dat" , np.array(nsat_young)) 

    if model == "hod":

        nsat = []
        ncen = []
  
        for i in xrange(len(sample)):
         
          print i
       
          mod.param_dict['logM0'] =  sample[i][0]
          mod.param_dict['sigma_logM'] =  sample[i][1]
          mod.param_dict['logMmin'] =  sample[i][2]
	  mod.param_dict['alpha'] =  sample[i][3]
	  mod.param_dict['logM1'] =  sample[i][4]

          ncen.append(mod.mean_occupation_centrals(prim_haloprop = mass))  
          nsat.append(mod.mean_occupation_satellites(prim_haloprop = mass))  

        np.savetxt("ncen_"+obs+"_"+model+"_"+str(Mr)+".dat" , np.array(ncen)) 
        np.savetxt("nsat_"+obs+"_"+model+"_"+str(Mr)+".dat" , np.array(nsat)) 
 
    return None 


def plot_occupations_decs(obs = "gmf", galtype = "ncen", clotter = True):

    npts = 2e3
    mass = np.logspace(10, 14, npts)
    pretty_color = prettycolors()
    if clotter == False:

        print "running compute_occupation_prediction first"

    else:

       if galtype == "ncen":

           ylabel = r'$\langle N_{\rm cen}(M_{h})\rangle$'
           xlabel = r'$M_{h}\; [h^{-1} \; M_{\odot}]$'
       if galtype == "nsat":

           ylabel = r'$\langle N_{\rm sat}(M_{h})\rangle$'
           xlabel = r'$M_{h}\; [h^{-1} \; M_{\odot}]$'

       old_18 = np.loadtxt("results/"+galtype+"_old_"+obs+"_dec_18.0.dat")
       young_18 = np.loadtxt("results/"+galtype+"_young_"+obs+"_dec_18.0.dat")
       old_19 = np.loadtxt("results/"+galtype+"_old_"+obs+"_dec_19.0.dat")
       young_19 = np.loadtxt("results/"+galtype+"_young_"+obs+"_dec_19.0.dat")
       old_20 = np.loadtxt("results/"+galtype+"_old_"+obs+"_dec_20.0.dat")
       young_20 = np.loadtxt("results/"+galtype+"_young_"+obs+"_dec_20.0.dat") 
         
       fig = plt.figure(1, figsize=(15,5))
       gs = gridspec.GridSpec(1, 3)# height_ratios=[2.5, 1], width_ratios=[1,1])  
       pretty_colors=prettycolors()
      
       ##### MR18 DEC ######

       ax = plt.subplot(gs[0])
       a_old, b_old, c_old = np.percentile(old_18, [16, 50, 84], axis=0) 
       a_young, b_young, c_young = np.percentile(young_18, [16, 50, 84], axis=0) 

       ax.plot(mass, b_young, color=pretty_colors[1], linewidth=2.) 
       ax.plot(mass, b_old, color=pretty_colors[4], linewidth=2.)
       ax.fill_between(mass, a_young, c_young, color=pretty_colors[1], alpha=0.5, edgecolor=pretty_colors[1] ) 
       ax.fill_between(mass, a_old, c_old, color=pretty_colors[4], alpha=0.5, edgecolor=pretty_colors[4]) 
 
       blue_line = mlines.Line2D([], [], ls = '-', c = pretty_color[1], linewidth=2, label = '$type-1$ $halos$')
       red_line = mlines.Line2D([], [], ls = '-', c = pretty_color[4], linewidth=2, label = '$type-2$ $halos$')
       import matplotlib as mpl
       #mpl.rc('font',family='serif')
       plt.legend(handles=[blue_line, red_line], frameon=False, loc= 2 , fontsize=15) 
       ax.set_ylabel(ylabel, fontsize=27)
       ax.set_xlabel(xlabel, fontsize=27)
       ax.set_yscale('log') 
       ax.set_xscale('log')
       ax.set_xlim([10**10, 8. * 10**13]) 
       if galtype == "ncen":
           ax.set_ylim([0.1, 2]) 
           ax.text(4 * 10**12. , 1.3, r'$M_{r}<-18$', fontsize=18) 
       if galtype == "nsat":
           ax.set_ylim([0.1,100])
           ax.text(4 * 10**12. , 50, r'$M_{r}<-18$', fontsize=18) 
       ##### Mr19 DEC ######

       ax = plt.subplot(gs[1])
       a_old, b_old, c_old = np.percentile(old_19, [16, 50, 84], axis=0) 
       a_young, b_young, c_young = np.percentile(young_19, [16, 50, 84], axis=0) 

       ax.plot(mass, b_young, color=pretty_colors[1], linewidth=2.) 
       ax.plot(mass, b_old, color=pretty_colors[4], linewidth=2.)
       ax.fill_between(mass, a_young, c_young, color=pretty_colors[1], alpha=0.5, edgecolor=pretty_colors[1]) 
       ax.fill_between(mass, a_old, c_old, color=pretty_colors[4], alpha=0.5, edgecolor=pretty_colors[4]) 
 
       blue_line = mlines.Line2D([], [], ls = '-', c = pretty_color[3], linewidth=2, label = 'high-concentration halos')
       red_line = mlines.Line2D([], [], ls = '-', c = pretty_color[4], linewidth=2, label = 'low-concentration halos')
       ax.set_yscale('log') 
       ax.set_xscale('log')
       ax.set_yticklabels([])
       ax.set_xlabel(xlabel, fontsize=27)
       ax.set_xlim([10**10, 8.*10**13]) 
       if galtype == "ncen":
           ax.set_ylim([0.1, 2]) 
           ax.text(4 * 10**12. , 1.3, r'$M_{r}<-19$', fontsize=18) 
       if galtype == "nsat":
           ax.set_ylim([0.1,100])
           ax.text(4 * 10**12. , 50, r'$M_{r}<-19$', fontsize=18) 
       ##### Mr20 DEC #######      
 
       ax = plt.subplot(gs[2])
       a_old, b_old, c_old = np.percentile(old_20, [16, 50, 84], axis=0) 
       a_young, b_young, c_young = np.percentile(young_20, [16, 50, 84], axis=0) 

       ax.plot(mass, b_young, color=pretty_colors[1], linewidth=2.) 
       ax.plot(mass, b_old, color=pretty_colors[4], linewidth=2.)
       ax.fill_between(mass, a_young, c_young, color=pretty_colors[1], alpha=0.5, edgecolor=pretty_colors[1]) 
       ax.fill_between(mass, a_old, c_old, color=pretty_colors[4], alpha=0.5, edgecolor=pretty_colors[4]) 
 
       blue_line = mlines.Line2D([], [], ls = '-', c = pretty_color[3], linewidth=2, label = 'high-concentration halos')
       red_line = mlines.Line2D([], [], ls = '-', c = pretty_color[4], linewidth=2, label = 'low-concentration halos')
       
       ax.set_yscale('log') 
       ax.set_xscale('log')
       ax.set_yticklabels([])
       ax.set_xlabel(xlabel, fontsize=27)
       ax.set_xlim([10**10, 8. * 10**13]) 
       if galtype == "ncen":
           ax.set_ylim([0.1, 2]) 
           ax.text(4 * 10**12. , 1.3, r'$M_{r}<-20$', fontsize=18) 
       if galtype == "nsat":
           ax.set_ylim([0.1,100])
           ax.text(4 * 10**12. , 50, r'$M_{r}<-20$', fontsize=18) 
       fig.subplots_adjust(wspace=0.0, hspace=0.0)
       fig_name = ''.join([ut.fig_dir(), 
        'paper', 
        '.abhod.',
        galtype,
        obs,
        '.pdf'])
       fig.savefig(fig_name, bbox_inches='tight')
       plt.close()
       
    return None 

if __name__=='__main__':

   directory = "/export/bbq2/mj/chains/"
   filename = "group_mcmc_chain_Mr18.0.hdf5"
   filename = directory+filename
   Mr = 18.0
   obs , model = "gmf" , "dec"
   nchains = 2050
   nburnins = 2040
   #occupation_predictions(filename, Mr, nburnins, nchains , obs , model) 
   plot_occupations_decs(obs = "gmf", galtype = "nsat", clotter = True)
