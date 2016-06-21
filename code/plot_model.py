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

halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')
import os
import os.path as path 

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

def model_predictions(filename, Mr, nburnins, nchains , obs = "wp", model = "dec"):

    if obs == "wp":
       if model == "dec":
          mod = dec_wp(Mr)
       if model == "hod":
          mod = dec_gmf(Mr)
    
    if obs == "gmf":
       if model == "dec":
          mod = dec_wp(Mr)
       if model == "hod":
          mod = hod_gmf(Mr)

    sample = h5py.File(filename , "r")["mcmc"][nburnins:nchains]
    sample = sample.reshape(sample.shape[0]*sample.shape[1] , sample.shape[2]) 
    model_obs = []
 
    for i in xrange(len(sample)):
        print i
        model_obs.append(mod._sum_stat(sample[i] , prior_range = None))
    np.savetxt(obs+"_"+model+"_"+str(Mr)+".dat" , np.array(model_obs)) 
 
    return None 

def plot_model_prediction(obs = "wp", model= "dec", clotter = True):

    if clotter == False:

        print "running compute_model_prediction first"

    else:

        #### loading the model predictions for the observables ######         

        file18 = np.loadtxt(obs+"_"+model+"_"+"18.dat")
        file19 = np.loadtxt(obs+"_"+model+"_"+"19.dat")
        file20 = np.loadtxt(obs+"_"+model+"_"+"20.dat")

        #### loading the observables in Mr 18,19,20

        if obs == "wp":

            rbin18 = np.loadtxt(bin)    
            rbin19 = np.loadtxt(bin)    
            rbin20 = np.loadtxt(bin)

	    wbin18 = 1.
	    wbin19 = 1.
	    wbin20 = 1.
            
	    data18 = data.load_wp(18.)    
            data19 = data.load_wp(19.)    
            data20 = data.load_wp(20.)
            
            err18 = np.diag(data.load_wp_covariance(18.)) 
            err19 = np.diag(data.load_wp_covariance(19.)) 
            err20 = np.diag(data.load_wp_covariance(20.)) 
    
        if obs == "gmf":
 
            cat18 = np.loadtxt("../dat/gmf_mr18.0.dat")
            cat19 = np.loadtxt("../dat/gmf_mr19.0.dat")   
            cat20 = np.loadtxt("../dat/gmf_mr20.0.dat")

            rbin18 = np.hstack([cat18[:,0], cat18[-1,1]])
            rbin19 = np.hstack([cat19[:,0], cat19[-1,1]])
            rbin20 = np.hstack([cat20[:,0], cat20[-1,1]])
           
	    wbin18 = rbin18[1:] - rbin18[:-1]
	    wbin19 = rbin19[1:] - rbin19[:-1]
	    wbin20 = rbin20[1:] - rbin20[:-1]
 
            data18 = data_group.load_gmf(18.)    
            data19 = data_group.load_gmf(19.)    
            data20 = data_group.load_gmf(20.)
            
            err18 = data_group.load_gmf_covariance(18., pois=True) 
            err19 = data_group.load_gmf_covariance(19., pois=True) 
            err20 = data_group.load_gmf_covariance(20., pois=True)

	prettyplot()
    	pretty_colors=prettycolors()
    	fig = plt.figure(1, figsize=(16,12))
    	gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1], width_ratios=[1,1])  

        ax = plt.subplot(gs[0])
        a, b, c, d, e = np.percentile(file18, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin18, a, e, color=pretty_colors[3], alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin18, b, d, color=pretty_colors[3], alpha=0.5, edgecolor="none")

        ax.errorbar(rbin18, data18, err18, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)

        ax.scatter(rbin18, data18, c='k', s=10, lw=0)
        ax.set_ylabel(r'$w_{\rm p}(r_{\rm p}) \; [rm{Mpc}h^{-1}]$', fontsize=27)
        ax.set_yscale('log') 
        ax.set_xscale('log')
        ax.set_xticklabels([])
        ax.set_xlim([0.05, 25.])
        ax.set_ylim([0.09, 1000.])


        ax = plt.subplot(gs[1])
        a, b, c, d, e = np.percentile(file19, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin19, a, e, color=pretty_colors[3], alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin19, b, d, color=pretty_colors[3], alpha=0.5, edgecolor="none")

        ax.errorbar(rbin19, data19, err19, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)

        ax.scatter(rbin19, data19, c='k', s=10, lw=0)
        ax.set_yscale('log') 
        ax.set_xscale('log')
        ax.set_xticklabels([])
        ax.set_xlim([0.05, 25.])
        ax.set_ylim([0.09, 1000.])
        
        ax = plt.subplot(gs[2])
        a, b, c, d, e = np.percentile(file20, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin20, a, e, color=pretty_colors[3], alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin20, b, d, color=pretty_colors[3], alpha=0.5, edgecolor="none")

        ax.errorbar(rbin20, data20, err20, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)

        ax.scatter(r_bin, xi_data, c='k', s=10, lw=0)
        ax.set_yscale('log') 
        ax.set_xscale('log')
        ax.set_xticklabels([])
        ax.set_xlim([0.05, 25.])
        ax.set_ylim([0.09, 1000.])

        fig.subplots_adjust(wspace=0.05, hspace=0.0)
        fig_name = ''.join([ut.fig_dir(), 
        'paper', 
        '.model.',
        model,
        obs,
        '.pdf'])
    fig.savefig(fig_name, bbox_inches='tight')
    plt.close()
    return None 


def plot_occupations(filename, Mr, nburnins, nchains, assembly = True , clotter = False):

    model = composite_model(Mr)
    npts = 1e3
    mass = np.logspace(11, 14, npts)
    prettyplot()
    pretty_colors=prettycolors()
    fig = plt.figure(1, figsize=(16,12))

    if (assembly == True):
        filename = 'Mr'+str(Mr)+'.hdf5'
    else:
        filename = 'adhoc_Mr'+str(Mr)+'.hdf5'

    sample = h5py.File(filename , "r")["mcmc"]
    npars = sample.shape[2]
    nwalkers = sample.shape[1]
    sample = sample[nchains-2:nchains, : , :]
    print sample.shape 
    sample = sample.reshape(sample.shape[0]*sample.shape[1] , sample.shape[2])
    if (clotter == False):
    	nsat_old = []
   	nsat_young = []
    	ncen_old = []
    	ncen_young = []

    	for i in xrange(len(sample)):
        	print i     
        	model.param_dict['logM0'] =  sample[i][0]
		model.param_dict['sigma_logM'] =  sample[i][1]
		model.param_dict['logMmin'] =  sample[i][2]
		model.param_dict['alpha'] =  sample[i][3]
		model.param_dict['logM1'] =  sample[i][4]
        	model.param_dict['mean_occupation_centrals_assembias_param1'] = sample[i][5]
        	model.param_dict['mean_occupation_satellites_assembias_param1'] = sample[i][6]
        	ncen_old.append(model.mean_occupation_centrals(prim_haloprop = mass, sec_haloprop_percentile=0))
        	ncen_young.append(model.mean_occupation_centrals(prim_haloprop = mass, sec_haloprop_percentile=1))
        	nsat_old.append(model.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=0))
        	nsat_young.append(model.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=1))
    	nsat_old = np.array(nsat_old)
    	nsat_young = np.array(nsat_young)
    	ncen_old = np.array(ncen_old)
    	ncen_young = np.array(ncen_young)
    	np.savetxt("nsat_old"+str(Mr)+".dat" , nsat_old)
    	np.savetxt("nsat_young"+str(Mr)+".dat" , nsat_young)
    	np.savetxt("ncen_old"+str(Mr)+".dat" , ncen_old)
    	np.savetxt("ncen_young"+str(Mr)+".dat" , ncen_young)

    else:

       nsat_old = np.loadtxt("nsat_old"+str(Mr)+".dat")
       nsat_young = np.loadtxt("nsat_young"+str(Mr)+".dat")
       ncen_old = np.loadtxt("ncen_old"+str(Mr)+".dat")
       ncen_young = np.loadtxt("ncen_young"+str(Mr)+".dat")
       
       a1, b1, c1 = np.percentile(nsat_young, [16, 50, 84], axis=0)
       a2, b2, c2 = np.percentile(nsat_old, [16, 50, 84], axis=0)
       a3, b3, c3 = np.percentile(ncen_young, [16, 50, 84], axis=0)
       a4, b4, c4 = np.percentile(ncen_old, [16, 50, 84], axis=0)
        
       fig = plt.figure()
       ax = fig.add_subplot(111)

       xlabel = ax.set_xlabel(r'$M_{\rm vir} [M_{\odot}]$', fontsize=25)
       ylabel = ax.set_ylabel(r'$\langle N_{\rm s}\rangle$', fontsize=25)
       #title = ax.set_title(r'$\langle N_{\rm s} \rangle$ for $\mathrm{M_{r}}$ < -'+str(Mr), fontsize=20)

       
       ax.plot(mass, b1, color='blue', linewidth=3.5)
       ax.fill_between(mass, a1, c1 , color='blue', alpha = 0.1)
       ax.plot(mass, b2, color='red', linewidth=3.5)
       ax.fill_between(mass, a2, c2 , color='red', alpha = 0.1)
       ax.plot(mass, 0.5*(b1 + b2), '--', color='k', linewidth=2.5)

       plt.loglog()
       plt.xlim(xmin=1e11, xmax=1e14)
       plt.ylim(ymin=5e-3, ymax=100)
       plt.xticks(fontsize=20)
       plt.yticks(fontsize=20)

       blue_line = mlines.Line2D([], [], ls = '-', c = 'b', linewidth=3, label = 'high-concentration halos')
       red_line = mlines.Line2D([], [], ls = '-', c = 'r', linewidth=3, label = 'low-concentration halos')
       black_line = mlines.Line2D([], [], ls = '--', c = 'k', linewidth=3, label = 'all halos')
       first_legend = plt.legend(handles=[blue_line, red_line, black_line], frameon=False, loc='best', fontsize=17)
       fig.savefig("nsats"+str(Mr)+".pdf", bbox_extra_artists=[xlabel, ylabel], bbox_inches='tight')
       
       fig = plt.figure()
       ax = fig.add_subplot(111)

       xlabel = ax.set_xlabel(r'$M_{\rm vir} [M_{\odot}]$', fontsize=25)
       ylabel = ax.set_ylabel(r'$\langle N_{\rm c}\rangle$', fontsize=25)
       #title = ax.set_title(r'$\langle N_{\rm c} \rangle$ for $\mathrm{M_{r}}$ < -'+str(Mr), fontsize=20)

       
       ax.plot(mass, b3, color='blue', linewidth=3.5)
       ax.fill_between(mass, a3, c3 , color='blue', alpha = 0.1)
       ax.plot(mass, b4, color='red', linewidth=3.5)
       ax.fill_between(mass, a4, c4 , color='red', alpha = 0.1)
       ax.plot(mass, 0.5*(b3 + b4), '--', color='k', linewidth=2.5)

       plt.loglog()
       plt.xlim(xmin=1e11, xmax=1e14)
       plt.ylim(ymin=5e-3, ymax=100)
       plt.xticks(fontsize=20)
       plt.yticks(fontsize=20)

       blue_line = mlines.Line2D([], [], ls = '-', c = 'b', linewidth=3, label = 'high-concentration halos')
       red_line = mlines.Line2D([], [], ls = '-', c = 'r', linewidth=3, label = 'low-concentration halos')
       black_line = mlines.Line2D([], [], ls = '--', c = 'k', linewidth=3, label = 'all halos')
       first_legend = plt.legend(handles=[blue_line, red_line, black_line], frameon=False, loc='best', fontsize=17)
       fig.savefig("ncens"+str(Mr)+".pdf", bbox_extra_artists=[xlabel, ylabel], bbox_inches='tight')
    return None 

if __name__=='__main__':

   directory = "/export/bbq2/mj/chains/"
   filename = "mcmc_chain_Mr19.0.hdf5"
   filename = directory+filename
   Mr = 19.0
   obs , model = "wp" , "dec"
   nchains = 22000
   nburnins = 21998
   model_predictions(filename, Mr, nburnins, nchains, obs, model)   

