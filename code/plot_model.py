'''
Genrating posterior model prediction figures for the paper

author: MJ

pretty_colors[1] for decorated model and 
pretty_colors[5] for standard hod models

'''
from halotools.sim_manager import CachedHaloCatalog
import os.path as path
import util as ut
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
import matplotlib.pyplot as plt
from matplotlib import lines as mlines
from matplotlib import gridspec
from scipy.stats import norm
from matplotlib.colors import colorConverter
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors

halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')
import os
import os.path as path 
from numpy.linalg import solve
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
import data
import data_group
def model_predictions(filename, Mr, nburnins, nchains , obs = "wp", model = "dec"):

    if obs == "wp":
       if model == "dec":
          mod = dec_wp(Mr)
       if model == "hod":
          mod = hod_wp(Mr)
    
    if obs == "gmf":
       if model == "dec":
          mod = dec_gmf(Mr)
       if model == "hod":
          mod = hod_gmf(Mr)

    sample = h5py.File(filename , "r")["mcmc"][nburnins:nchains]
    sample = sample.reshape(sample.shape[0]*sample.shape[1] , sample.shape[2]) 
    model_obs = []
 
    for i in xrange(len(sample)):
        print i
        model_obs.append(mod._sum_stat(sample[i] , prior_range = None)[1])
    np.savetxt("results/"+obs+"_"+model+"_"+str(Mr)+".dat" , np.array(model_obs)) 
 
    return None 




def plot_model_prediction(obs = "wp", model= "dec", clotter = True):

    if clotter == False:

        print "running compute_model_prediction first"

    else:

        #### loading the model predictions for the observables ######         

        file18 = np.loadtxt("results/"+obs+"_"+model+"_"+"18.0.dat")
        file19 = np.loadtxt("results/"+obs+"_"+model+"_"+"19.0.dat")
        file20 = np.loadtxt("results/"+obs+"_"+model+"_"+"20.0.dat")

        #### loading the observables in Mr 18,19,20

        if obs == "wp":

            rbin18 = np.loadtxt("bin.dat")    
            rbin19 = np.loadtxt("bin.dat")    
            rbin20 = np.loadtxt("bin.dat")

            rbin18 = np.mean(rbin18, axis = 1)             
            rbin19 = np.mean(rbin19, axis = 1)             
            rbin20 = np.mean(rbin20, axis = 1)             

	    wbin18 = 1.
	    wbin19 = 1.
	    wbin20 = 1.
            
	    data18 = data.load_wp(18.)    
            data19 = data.load_wp(19.)    
            data20 = data.load_wp(20.)
            
            err18 = np.diag(data.load_wp_covariance(18.))**.5 
            err19 = np.diag(data.load_wp_covariance(19.))**.5 
            err20 = np.diag(data.load_wp_covariance(20.))**.5 
    
        if obs == "gmf":
 
            cat18 = np.loadtxt("../dat/gmf_mr18.0.dat")
            cat19 = np.loadtxt("../dat/gmf_mr19.0.dat")   
            cat20 = np.loadtxt("../dat/gmf_mr20.0.dat")

            rbin18 = np.hstack([cat18[:,0], cat18[-1,1]])
            rbin19 = np.hstack([cat19[:,0], cat19[-1,1]])
            rbin20 = np.hstack([cat20[:,0], cat20[-1,1]])
           
            rbin18 = 0.5 * (rbin18[1:] + rbin18[:-1])	    
            rbin19 = 0.5 * (rbin19[1:] + rbin19[:-1])	    
            rbin20 = 0.5 * (rbin20[1:] + rbin20[:-1])	    

            wbin18 = rbin18[1:] - rbin18[:-1]
	    wbin19 = rbin19[1:] - rbin19[:-1]
	    wbin20 = rbin20[1:] - rbin20[:-1]
 
            data18 = data_group.load_gmf(18.)    
            data19 = data_group.load_gmf(19.)    
            data20 = data_group.load_gmf(20.)
            
            err18 = data_group.load_gmf_covariance(18., pois=True)**.5 
            err19 = data_group.load_gmf_covariance(19., pois=True)**.5 
            err20 = data_group.load_gmf_covariance(20., pois=True)**.5

	#prettyplot()
    	pretty_colors=prettycolors()
    	fig = plt.figure(1, figsize=(24,8))
    	gs = gridspec.GridSpec(1,3)#, height_ratios=[1, 1], width_ratios=[1,1])  

        ax = plt.subplot(gs[0])
        a, b, c, d, e = np.percentile(file18, [2.5, 16, 50, 84, 97.5], axis=0) 
        
        print a.shape
        ax.fill_between(rbin18, a, e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin18, b, d, color=pretty_colors[1], alpha=0.7, edgecolor="none")

        ax.errorbar(rbin18, data18, yerr=err18, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin18, data18, c='k', s=10, lw=0)
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          #ax.set_xticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([0.5, 1000.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(5, 500, r'$M_{r}<-18$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1.5, 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks([2 , 10, 50],[2,10,50],fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(20, .002, r'$M_{r}<-18$', fontsize=25) 
        ax = plt.subplot(gs[1])
        a, b, c, d, e = np.percentile(file19, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin19, a, e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin19, b, d, color=pretty_colors[1], alpha=0.7, edgecolor="none")

        ax.errorbar(rbin19, data19, yerr=err19, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)

        ax.scatter(rbin19, data19, c='k', s=10, lw=0)
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          plt.xticks(fontsize=20)
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([0.5, 1000.])
          ax.text(5, 500, r'$M_{r}<-19$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          plt.xticks([2 , 10, 50],[2,10,50],fontsize=20)
          plt.xticks(fontsize=20)
          ax.set_xlim([1.5, 99.])
          ax.set_ylim([10**-8., 10**-2.])
          ax.text(20, .002, r'$M_{r}<-19$', fontsize=25) 
        ax = plt.subplot(gs[2])
        a, b, c, d, e = np.percentile(file20, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin20, a, e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin20, b, d, color=pretty_colors[1], alpha=0.7, edgecolor="none")

        ax.errorbar(rbin20, data20, yerr=err20, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)

        ax.scatter(rbin20, data20, c='k', s=10, lw=0)
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          plt.xticks(fontsize=20)
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([0.5, 1000.])
          ax.text(5, 500, r'$M_{r}<-20$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          plt.xticks([2 , 10, 50],[2,10,50],fontsize=20)
          plt.xticks(fontsize=20)
          ax.set_xlim([1.5, 99.])
          ax.set_ylim([10**-8., 10**-2.])
          ax.text(20, .002, r'$M_{r}<-20$', fontsize=25) 
        fig.subplots_adjust(wspace=0.0, hspace=0.0)
        fig_name = ''.join([ut.fig_dir(), 
        'paper', 
        '.model.',
        model,
        obs,
        '.pdf'])
    fig.savefig(fig_name, bbox_inches='tight')
    plt.close()
    return None 

def plot_wp_prediction(obs = "wp", model= "dec", clotter = True):

    if clotter == False:

        print "running compute_model_prediction first"

    else:

        #### loading the model predictions for the observables ######         

        file18 = np.loadtxt("results/"+obs+"_"+model+"_"+"18.0.dat")
        file185 = np.loadtxt("results/"+obs+"_"+model+"_"+"18.5.dat")
        file19 = np.loadtxt("results/"+obs+"_"+model+"_"+"19.0.dat")
        file195 = np.loadtxt("results/"+obs+"_"+model+"_"+"19.5.dat")
        file20 = np.loadtxt("results/"+obs+"_"+model+"_"+"20.0.dat")
        file205 = np.loadtxt("results/"+obs+"_"+model+"_"+"20.5.dat")

        #### loading the observables in Mr 18,19,20

        if obs == "wp":

            rbin18 = np.loadtxt("bin.dat")    
            rbin185 = np.loadtxt("bin.dat")    
            rbin19 = np.loadtxt("bin.dat")    
            rbin195 = np.loadtxt("bin.dat")    
            rbin20 = np.loadtxt("bin.dat")
            rbin205 = np.loadtxt("bin.dat")

            rbin18 = np.mean(rbin18, axis = 1)             
            rbin185 = np.mean(rbin185, axis = 1)             
            rbin19 = np.mean(rbin19, axis = 1)             
            rbin195 = np.mean(rbin195, axis = 1)             
            rbin20 = np.mean(rbin20, axis = 1)             
            rbin205 = np.mean(rbin205, axis = 1)             

	    wbin18 = 1.
	    wbin185 = 1.
	    wbin19 = 1.
	    wbin195 = 1.
	    wbin20 = 1.
	    wbin205 = 1.
            
	    data18 = data.load_wp(18.)    
	    data185 = data.load_wp(18.5)    
            data19 = data.load_wp(19.)    
	    data195 = data.load_wp(19.5)    
            data20 = data.load_wp(20.)
            data205 = data.load_wp(20.5)
            
            err18 = np.diag(data.load_wp_covariance(18.))**.5 
            err185 = np.diag(data.load_wp_covariance(18.5))**.5 
            err19 = np.diag(data.load_wp_covariance(19.))**.5 
            err195 = np.diag(data.load_wp_covariance(19.5))**.5 
            err20 = np.diag(data.load_wp_covariance(20.))**.5 
            err205 = np.diag(data.load_wp_covariance(20.5))**.5 
    
        if obs == "gmf":
 
            cat18 = np.loadtxt("../dat/gmf_mr18.0.dat")
            cat19 = np.loadtxt("../dat/gmf_mr19.0.dat")   
            cat20 = np.loadtxt("../dat/gmf_mr20.0.dat")

            rbin18 = np.hstack([cat18[:,0], cat18[-1,1]])
            rbin19 = np.hstack([cat19[:,0], cat19[-1,1]])
            rbin20 = np.hstack([cat20[:,0], cat20[-1,1]])
           
            rbin18 = 0.5 * (rbin18[1:] + rbin18[:-1])	    
            rbin19 = 0.5 * (rbin19[1:] + rbin19[:-1])	    
            rbin20 = 0.5 * (rbin20[1:] + rbin20[:-1])	    

            wbin18 = rbin18[1:] - rbin18[:-1]
	    wbin19 = rbin19[1:] - rbin19[:-1]
	    wbin20 = rbin20[1:] - rbin20[:-1]
 
            data18 = data_group.load_gmf(18.)    
            data19 = data_group.load_gmf(19.)    
            data20 = data_group.load_gmf(20.)
            
            err18 = data_group.load_gmf_covariance(18., pois=True)**.5 
            err19 = data_group.load_gmf_covariance(19., pois=True)**.5 
            err20 = data_group.load_gmf_covariance(20., pois=True)**.5

	#prettyplot()
    	pretty_colors=prettycolors()
    	fig = plt.figure(1, figsize=(24,12))
    	gs = gridspec.GridSpec(2,3)#, height_ratios=[1, 1], width_ratios=[1,1])  

        ax = plt.subplot(gs[0,0]) #Mr18

        a, b, c, d, e = np.percentile(file18, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin18, a, e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin18, b, d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        ax.errorbar(rbin18, data18, yerr=err18, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin18, data18, c='k', s=10, lw=0)
        if obs == "wp":
          #ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([0.5, 900.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(3.5, 500, r'$M_{r}<-18$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25) 

        ax = plt.subplot(gs[0,1]) #Mr18.5

        a, b, c, d, e = np.percentile(file185, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin185, a, e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin185, b, d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        ax.errorbar(rbin185, data185, yerr=err185, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin185, data185, c='k', s=10, lw=0)
        if obs == "wp":
          #ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xticklabels([])
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([0.5, 900.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(3.5, 500, r'$M_{r}<-18.5$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25) 

        ax = plt.subplot(gs[0,2]) #Mr19.

        a, b, c, d, e = np.percentile(file19, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin19, a, e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin19, b, d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        ax.errorbar(rbin19, data19, yerr=err19, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin19, data19, c='k', s=10, lw=0)
        if obs == "wp":
          #ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([0.5, 900.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(3.5, 500, r'$M_{r}<-19$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([.5, 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25) 

        ax = plt.subplot(gs[1,0]) #Mr19.5

        a, b, c, d, e = np.percentile(file195, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin195, a, e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin195, b, d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        ax.errorbar(rbin195, data195, yerr=err195, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin195, data195, c='k', s=10, lw=0)
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          #ax.set_xticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([0.5, 900.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(3.5, 500, r'$M_{r}<-19.5$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25)


        ax = plt.subplot(gs[1,1]) #Mr20.0

        a, b, c, d, e = np.percentile(file20, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin20, a, e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin20, b, d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        ax.errorbar(rbin20, data20, yerr=err20, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin20, data20, c='k', s=10, lw=0)
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([0.5, 900.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(3.5, 500, r'$M_{r}<-20$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25)


        ax = plt.subplot(gs[1,2]) #Mr20.5
        #plt.axis('off')
 
        a, b, c, d, e = np.percentile(file205, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin205, a, e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin205, b, d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        ax.errorbar(rbin205, data205, yerr=err205, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin205, data205, c='k', s=10, lw=0)
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([0.5, 900.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(3.5, 500, r'$M_{r}<-20.5$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-20.5$', fontsize=25)
        fig.subplots_adjust(wspace=0.0, hspace=0.0)
        fig_name = ''.join([ut.fig_dir(), 
        'paper', 
        '.model.',
        model,
        obs,
        '.pdf'])

    fig.savefig(fig_name, bbox_inches='tight')
    plt.close()
    return None 

def plot_wpres_prediction(obs = "wp", clotter = True):

    if clotter == False:

        print "running compute_model_prediction first"

    else:

        #### loading the model predictions for the observables ######         

        file18dec = np.loadtxt("results/"+obs+"_dec_"+"18.0.dat")
        file185dec = np.loadtxt("results/"+obs+"_dec_"+"18.5.dat")
        file19dec = np.loadtxt("results/"+obs+"_dec_"+"19.0.dat")
        file195dec = np.loadtxt("results/"+obs+"_dec_"+"19.5.dat")
        file20dec = np.loadtxt("results/"+obs+"_dec_"+"20.0.dat")
        file205dec = np.loadtxt("results/"+obs+"_dec_"+"20.5.dat")
        file21dec = np.loadtxt("results/"+obs+"_dec_"+"21.0.dat")
        file215dec = np.loadtxt("results/"+obs+"_dec_"+"21.5.dat")
        file18hod = np.loadtxt("results/"+obs+"_hod_"+"18.0.dat")
        file185hod = np.loadtxt("results/"+obs+"_hod_"+"18.5.dat")
        file19hod = np.loadtxt("results/"+obs+"_hod_"+"19.0.dat")
        file195hod = np.loadtxt("results/"+obs+"_hod_"+"19.5.dat")
        file20hod = np.loadtxt("results/"+obs+"_hod_"+"20.0.dat")
        file205hod = np.loadtxt("results/"+obs+"_hod_"+"20.5.dat")
        file21hod = np.loadtxt("results/"+obs+"_hod_"+"21.0.dat")
        file215hod = np.loadtxt("results/"+obs+"_hod_"+"21.5.dat")

        #### loading the observables in Mr 18,19,20

        if obs == "wp":

            rbin18 = np.loadtxt("bin.dat")    
            rbin185 = np.loadtxt("bin.dat")    
            rbin19 = np.loadtxt("bin.dat")    
            rbin195 = np.loadtxt("bin.dat")    
            rbin20 = np.loadtxt("bin.dat")
            rbin205 = np.loadtxt("bin.dat")
            rbin21 = np.loadtxt("bin.dat")
            rbin215 = np.loadtxt("bin.dat")

            rbin18 = np.mean(rbin18, axis = 1)             
            rbin185 = np.mean(rbin185, axis = 1)             
            rbin19 = np.mean(rbin19, axis = 1)             
            rbin195 = np.mean(rbin195, axis = 1)             
            rbin20 = np.mean(rbin20, axis = 1)             
            rbin205 = np.mean(rbin205, axis = 1)             
            rbin21 = np.mean(rbin21, axis = 1)             
            rbin215 = np.mean(rbin215, axis = 1)             

	    wbin18 = 1.
	    wbin185 = 1.
	    wbin19 = 1.
	    wbin195 = 1.
	    wbin20 = 1.
	    wbin205 = 1.
	    wbin21 = 1.
	    wbin215 = 1.
            
	    data18 = data.load_wp(18.)    
	    data185 = data.load_wp(18.5)    
            data19 = data.load_wp(19.)    
	    data195 = data.load_wp(19.5)    
            data20 = data.load_wp(20.)
            data205 = data.load_wp(20.5)
            data21 = data.load_wp(21.0)
            data215 = data.load_wp(21.5)
            
            err18 = np.diag(data.load_wp_covariance(18.))**.5 
            err185 = np.diag(data.load_wp_covariance(18.5))**.5 
            err19 = np.diag(data.load_wp_covariance(19.))**.5 
            err195 = np.diag(data.load_wp_covariance(19.5))**.5 
            err20 = np.diag(data.load_wp_covariance(20.))**.5 
            err205 = np.diag(data.load_wp_covariance(20.5))**.5 
            err21 = np.diag(data.load_wp_covariance(21.0))**.5 
            err215 = np.diag(data.load_wp_covariance(21.5))**.5 
    
        if obs == "gmf":
 
            cat18 = np.loadtxt("../dat/gmf_mr18.0.dat")
            cat19 = np.loadtxt("../dat/gmf_mr19.0.dat")   
            cat20 = np.loadtxt("../dat/gmf_mr20.0.dat")

            rbin18 = np.hstack([cat18[:,0], cat18[-1,1]])
            rbin19 = np.hstack([cat19[:,0], cat19[-1,1]])
            rbin20 = np.hstack([cat20[:,0], cat20[-1,1]])
           
            rbin18 = 0.5 * (rbin18[1:] + rbin18[:-1])	    
            rbin19 = 0.5 * (rbin19[1:] + rbin19[:-1])	    
            rbin20 = 0.5 * (rbin20[1:] + rbin20[:-1])	    

            wbin18 = rbin18[1:] - rbin18[:-1]
	    wbin19 = rbin19[1:] - rbin19[:-1]
	    wbin20 = rbin20[1:] - rbin20[:-1]
 
            data18 = data_group.load_gmf(18.)    
            data19 = data_group.load_gmf(19.)    
            data20 = data_group.load_gmf(20.)
            
            err18 = data_group.load_gmf_covariance(18., pois=True)**.5 
            err19 = data_group.load_gmf_covariance(19., pois=True)**.5 
            err20 = data_group.load_gmf_covariance(20., pois=True)**.5

	#prettyplot()
    	pretty_colors=prettycolors()
    	fig = plt.figure(1, figsize=(24,16))
    	gs = gridspec.GridSpec(2,4)#, height_ratios=[1, 1], width_ratios=[1,1])  

        ax = plt.subplot(gs[0,0]) #Mr18

        a, b, c, d, e = np.percentile(file18dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin18, a/data18 -1 , e/data18 -1, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin18, b/data18 -1, d/data18 -1, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file18hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin18, a/data18-1, e/data18-1, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin18, b/data18-1, d/data18-1, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin18, np.zeros_like(data18), yerr=err18/data18, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin18, np.zeros_like(data18), c='k', s=10, lw=0)

        if obs == "wp":
          #ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          ax.set_ylabel(r'$w_{p}^{\mathrm{model}}/w_{p}^{\mathrm{data}} \; -1 $', fontsize=27)
          #ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([-1., 1.2])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(2.0, -0.8, r'$M_{r}<-18$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25) 

        ax = plt.subplot(gs[0,1]) #Mr18.5

        a, b, c, d, e = np.percentile(file185dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin185, a/data185-1, e/data185-1, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin185, b/data185-1, d/data185-1, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file185hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin185, a/data185-1, e/data185-1, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin185, b/data185-1, d/data185-1, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin185, np.zeros_like(data185), yerr=err185/data185, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin185, np.zeros_like(data185), c='k', s=10, lw=0)

        if obs == "wp":
          #ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontize=27)
          #ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xticklabels([])
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([-1., 1.2])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(2.0, -0.8, r'$M_{r}<-18.5$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25) 

        ax = plt.subplot(gs[0,2]) #Mr19.

        a, b, c, d, e = np.percentile(file19dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin19, a/data19-1, e/data19-1, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin19, b/data19-1, d/data19-1, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file19hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin19, a/data19-1, e/data19-1, color='#ee6a50', alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin19, b/data19-1, d/data19-1, color='#ee6a50', alpha=0.7, edgecolor="none")
        ax.errorbar(rbin19, np.zeros_like(data19), yerr=err19/data19, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin19, np.zeros_like(data19), c='k', s=10, lw=0)
        if obs == "wp":
          #ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          #ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([-1, 1.2])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(2.0, -0.8, r'$M_{r}<-19$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([.5, 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25) 

        ax = plt.subplot(gs[0,3]) #Mr19.5

        a, b, c, d, e = np.percentile(file195dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin195, a/data195-1, e/data195-1, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin195, b/data195-1, d/data195-1, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file195hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin195, a/data195-1, e/data195-1, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin195, b/data195-1, d/data195-1, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin195, np.zeros_like(data195), yerr=err195/data195, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin195, np.zeros_like(data195), c='k', s=10, lw=0)
        if obs == "wp":
          #ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          #ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([-1., 1.2])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(2.0, -0.8, r'$M_{r}<-19.5$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25)


        ax = plt.subplot(gs[1,0]) #Mr20.0

        a, b, c, d, e = np.percentile(file20dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin20, a/data20-1, e/data20-1, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin20, b/data20-1, d/data20-1, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file20hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin20, a/data20-1, e/data20-1, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin20, b/data20-1, d/data20-1, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin20, np.zeros_like(data20), yerr=err20/data20, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin20, np.zeros_like(data20), c='k', s=10, lw=0)
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          ax.set_ylabel(r'$w_{p}^{\mathrm{model}}/ w_{p}^{\mathrm{data}} \; -1 $', fontsize=27)
          #ax.set_yscale('log') 
          ax.set_xscale('log')
          #ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([-1., 1.2])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(2.0, -0.8, r'$M_{r}<-20$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25)


        ax = plt.subplot(gs[1,1]) #Mr20.5
        #plt.axis('off')
 
        a, b, c, d, e = np.percentile(file205dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin205, a/data205-1, e/data205-1, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin205, b/data205-1, d/data205-1, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file205hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin205, a/data205-1, e/data205-1, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin205, b/data205-1, d/data205-1, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin205, np.zeros_like(data205), yerr=err205/data205, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin205, np.zeros_like(data205), c='k', s=10, lw=0)
        
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          #ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([-1., 1.2])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(2.0, -0.8, r'$M_{r}<-20.5$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-20.5$', fontsize=25)
        
        ax = plt.subplot(gs[1,2]) #Mr21.0
        #plt.axis('off')
 
        a, b, c, d, e = np.percentile(file21dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin21, a/data21-1, e/data21-1, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin21, b/data21-1, d/data21-1, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file21hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin21, a/data21-1, e/data21-1, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin21, b/data21-1, d/data21-1, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin21, np.zeros_like(data21), yerr=err21/data21, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin21, np.zeros_like(data21), c='k', s=10, lw=0)
        
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          #ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([-1., 1.2])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(2.0, -0.8, r'$M_{r}<-21$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-21$', fontsize=25)
        
        ax = plt.subplot(gs[1,3]) #Mr21.5
        #plt.axis('off')
        
        a, b, c, d, e = np.percentile(file215dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin215, a/data215-1, e/data215-1, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin215, b/data215-1, d/data215-1, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file215hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin215, a/data215-1, e/data215-1, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin215, b/data215-1, d/data215-1, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin215, np.zeros_like(data215), yerr=err215/data215, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin215, np.zeros_like(data215), c='k', s=10, lw=0)
        
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          #ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([-1., 1.2])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(2.0, -0.8, r'$M_{r}<-21.5$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-21$', fontsize=25)
        
        fig.subplots_adjust(wspace=0.0, hspace=0.0)
        fig_name = ''.join([ut.fig_dir(), 
        'paper', 
        '.wpres',
        '.pdf'])
        
    fig.savefig(fig_name, bbox_inches='tight')
    plt.close()
    return None 

def plot_wpmodel_prediction(obs = "wp", clotter = True):

    if clotter == False:

        print "running compute_model_prediction first"

    else:

        #### loading the model predictions for the observables ######         

        file18dec = np.loadtxt("results/"+obs+"_dec_"+"18.0.dat")
        file185dec = np.loadtxt("results/"+obs+"_dec_"+"18.5.dat")
        file19dec = np.loadtxt("results/"+obs+"_dec_"+"19.0.dat")
        file195dec = np.loadtxt("results/"+obs+"_dec_"+"19.5.dat")
        file20dec = np.loadtxt("results/"+obs+"_dec_"+"20.0.dat")
        file205dec = np.loadtxt("results/"+obs+"_dec_"+"20.5.dat")
        file21dec = np.loadtxt("results/"+obs+"_dec_"+"21.0.dat")
        file215dec = np.loadtxt("results/"+obs+"_dec_"+"21.5.dat")
        file18hod = np.loadtxt("results/"+obs+"_hod_"+"18.0.dat")
        file185hod = np.loadtxt("results/"+obs+"_hod_"+"18.5.dat")
        file19hod = np.loadtxt("results/"+obs+"_hod_"+"19.0.dat")
        file195hod = np.loadtxt("results/"+obs+"_hod_"+"19.5.dat")
        file20hod = np.loadtxt("results/"+obs+"_hod_"+"20.0.dat")
        file205hod = np.loadtxt("results/"+obs+"_hod_"+"20.5.dat")
        file21hod = np.loadtxt("results/"+obs+"_hod_"+"21.0.dat")
        file215hod = np.loadtxt("results/"+obs+"_hod_"+"21.5.dat")

        #### loading the observables in Mr 18,19,20

        if obs == "wp":

            rbin18 = np.loadtxt("bin.dat")    
            rbin185 = np.loadtxt("bin.dat")    
            rbin19 = np.loadtxt("bin.dat")    
            rbin195 = np.loadtxt("bin.dat")    
            rbin20 = np.loadtxt("bin.dat")
            rbin205 = np.loadtxt("bin.dat")
            rbin21 = np.loadtxt("bin.dat")
            rbin215 = np.loadtxt("bin.dat")

            rbin18 = np.mean(rbin18, axis = 1)             
            rbin185 = np.mean(rbin185, axis = 1)             
            rbin19 = np.mean(rbin19, axis = 1)             
            rbin195 = np.mean(rbin195, axis = 1)             
            rbin20 = np.mean(rbin20, axis = 1)             
            rbin205 = np.mean(rbin205, axis = 1)             
            rbin21 = np.mean(rbin21, axis = 1)             
            rbin215 = np.mean(rbin215, axis = 1)             

	    wbin18 = 1.
	    wbin185 = 1.
	    wbin19 = 1.
	    wbin195 = 1.
	    wbin20 = 1.
	    wbin205 = 1.
	    wbin21 = 1.
	    wbin215 = 1.
            
	    data18 = data.load_wp(18.)    
	    data185 = data.load_wp(18.5)    
            data19 = data.load_wp(19.)    
	    data195 = data.load_wp(19.5)    
            data20 = data.load_wp(20.)
            data205 = data.load_wp(20.5)
            data21 = data.load_wp(21.0)
            data215 = data.load_wp(21.5)
            
            err18 = np.diag(data.load_wp_covariance(18.))**.5 
            err185 = np.diag(data.load_wp_covariance(18.5))**.5 
            err19 = np.diag(data.load_wp_covariance(19.))**.5 
            err195 = np.diag(data.load_wp_covariance(19.5))**.5 
            err20 = np.diag(data.load_wp_covariance(20.))**.5 
            err205 = np.diag(data.load_wp_covariance(20.5))**.5 
            err21 = np.diag(data.load_wp_covariance(21.0))**.5 
            err215 = np.diag(data.load_wp_covariance(21.5))**.5 
    
        if obs == "gmf":
 
            cat18 = np.loadtxt("../dat/gmf_mr18.0.dat")
            cat19 = np.loadtxt("../dat/gmf_mr19.0.dat")   
            cat20 = np.loadtxt("../dat/gmf_mr20.0.dat")

            rbin18 = np.hstack([cat18[:,0], cat18[-1,1]])
            rbin19 = np.hstack([cat19[:,0], cat19[-1,1]])
            rbin20 = np.hstack([cat20[:,0], cat20[-1,1]])
           
            rbin18 = 0.5 * (rbin18[1:] + rbin18[:-1])	    
            rbin19 = 0.5 * (rbin19[1:] + rbin19[:-1])	    
            rbin20 = 0.5 * (rbin20[1:] + rbin20[:-1])	    

            wbin18 = rbin18[1:] - rbin18[:-1]
	    wbin19 = rbin19[1:] - rbin19[:-1]
	    wbin20 = rbin20[1:] - rbin20[:-1]
 
            data18 = data_group.load_gmf(18.)    
            data19 = data_group.load_gmf(19.)    
            data20 = data_group.load_gmf(20.)
            
            err18 = data_group.load_gmf_covariance(18., pois=True)**.5 
            err19 = data_group.load_gmf_covariance(19., pois=True)**.5 
            err20 = data_group.load_gmf_covariance(20., pois=True)**.5

	#prettyplot()
    	pretty_colors=prettycolors()
    	fig = plt.figure(1, figsize=(24,16))
    	gs = gridspec.GridSpec(2,4)#, height_ratios=[1, 1], width_ratios=[1,1])  

        ax = plt.subplot(gs[0,0]) #Mr18

        a, b, c, d, e = np.percentile(file18dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin18, a , e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin18, b , d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file18hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin18, a, e, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin18, b, d, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin18, data18, yerr=err18, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin18, data18, c='k', s=10, lw=0)

        if obs == "wp":
          #ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([.3, 2000.0])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(.1, 0.8, r'$M_{r}<-18$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25) 

        ax = plt.subplot(gs[0,1]) #Mr18.5

        a, b, c, d, e = np.percentile(file185dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin185, a , e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin185, b , d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file185hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin185, a, e, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin185, b, d, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin185, data185, yerr=err185, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin185, data185, c='k', s=10, lw=0)

        if obs == "wp":
          #ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xticklabels([])
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([0.3, 2000.0])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(0.1, 0.8, r'$M_{r}<-18.5$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25) 

        ax = plt.subplot(gs[0,2]) #Mr19.

        a, b, c, d, e = np.percentile(file19dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin19, a , e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin19, b , d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file19hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin19, a, e, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin19, b, d, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin19, data19, yerr=err19, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin19, data19, c='k', s=10, lw=0)
        if obs == "wp":
          #ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([0.3, 2000.0])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(0.1, 0.8, r'$M_{r}<-19$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([.5, 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25) 

        ax = plt.subplot(gs[0,3]) #Mr19.5

        a, b, c, d, e = np.percentile(file195dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin195, a , e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin195, b , d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file195hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin195, a, e, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin195, b, d, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin195, data195, yerr=err195, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin195, data195, c='k', s=10, lw=0)
        if obs == "wp":
          #ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([0.3, 2000.0])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(0.1, 0.8, r'$M_{r}<-19.5$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25)


        ax = plt.subplot(gs[1,0]) #Mr20.0

        a, b, c, d, e = np.percentile(file20dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin20, a , e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin20, b , d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file20hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin20, a, e, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin20, b, d, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin20, data20, yerr=err20, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin20, data20, c='k', s=10, lw=0)
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          ax.set_ylabel(r'$w_{p} (r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          #ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([.3, 2000.0])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(.1, 0.8, r'$M_{r}<-20$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-18$', fontsize=25)


        ax = plt.subplot(gs[1,1]) #Mr20.5
        #plt.axis('off')
 
        a, b, c, d, e = np.percentile(file205dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin205, a , e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin205, b , d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file205hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin205, a, e, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin205, b, d, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin205, data205, yerr=err205, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin205, data205, c='k', s=10, lw=0)
        
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([.3, 2000.0])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(0.1, 0.8, r'$M_{r}<-20.5$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-20.5$', fontsize=25)
        
        ax = plt.subplot(gs[1,2]) #Mr21.0
        #plt.axis('off')
 
        a, b, c, d, e = np.percentile(file21dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin21, a , e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin21, b , d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file21hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin21, a, e, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin21, b, d, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin21, data21, yerr=err21, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin21, data21, c='k', s=10, lw=0)
        
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([.3, 2000.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(0.1, 0.8, r'$M_{r}<-21$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-21$', fontsize=25)
        
        ax = plt.subplot(gs[1,3]) #Mr21.5
        #plt.axis('off')
        
        a, b, c, d, e = np.percentile(file215dec, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin215, a , e, color=pretty_colors[1], alpha=0.4, edgecolor="none") 
        ax.fill_between(rbin215, b , d, color=pretty_colors[1], alpha=0.7, edgecolor="none")
        a, b, c, d, e = np.percentile(file215hod, [2.5, 16, 50, 84, 97.5], axis=0) 
        ax.fill_between(rbin215, a, e, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin215, b, d, color='#ee6a50', alpha=0.6, edgecolor="none")
        ax.errorbar(rbin215, data215, yerr=err215, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin215, data215, c='k', s=10, lw=0)
        ax.fill_between(rbin215, a, e, color='#ee6a50', alpha=0.3, edgecolor="none") 
        ax.fill_between(rbin215, b, d, color='#ee6a50', alpha=0.6, edgecolor="none")
        
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([.3, 2000.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(.1, 0.8, r'$M_{r}<-21.5$', fontsize=25) 
        if obs == "gmf":
          ax.set_xlabel(r'$N$', fontsize=27)
          ax.set_ylabel(r'$g(N) \; [\mathrm{Mpc}^{-3}\; h^{3}]$', fontsize=27)
          ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_xlim([1., 99.])
          ax.set_ylim([10**-8., 10**-2.])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(35, .002, r'$M_{r}<-21$', fontsize=25)
        
        fig.subplots_adjust(wspace=0.0, hspace=0.0)
        fig_name = ''.join([ut.fig_dir(), 
        'paper', 
        '.wpmodel',
        '.pdf'])
        
    fig.savefig(fig_name, bbox_inches='tight')
    plt.close()
    return None 

def compute_chisq_aic_bic(filename, nchains, nburnins, Mr, obs, model):

    import data as data_wp
    #import data_group as data_group

    sample = h5py.File(filename , "r")["mcmc"][nburnins:nchains]
    sample = sample.reshape(sample.shape[0]*sample.shape[1] , sample.shape[2])
    bestfit = np.median(sample, axis = 0)
    print bestfit
    bestfit[1]  = 0.01
    #for i in range(sample.shape[1]):
    #    n, b, patches = plt.hist(sample[:,i], 20)
    #    bin_max = np.argmax(n)
    #    print bin_max
    #    print b[bin_max]
    #    #bin_max = np.where(n == n.max())
    #    #print 'maxbin', b[bin_max][0]
    #    bestfit[i] = b[bin_max]
    #bestfit = np.mean(sample, axis = 0)
    #print bestfit
    
    bestfit = np.median(sample, axis = 0) 
    if obs == "wp":

       data_obs = data_wp.load_data(Mr)
       cov_obs = data_wp.load_covariance(Mr)

       if model == "dec":
          mod = dec_wp(Mr)
          npar = 7
       if model == "hod":
          mod = hod_wp(Mr)
          npar = 5
 
    if obs == "gmf":

       data_obs = data_group.load_data(Mr)
       cov_obs = data_group.load_covariance(Mr, pois = True)

       if model == "dec":
          mod = dec_gmf(Mr)
          npar = 7
       if model == "hod":
          mod = hod_gmf(Mr)
          npar = 5

    n_ens , w_ens = [] , []

    #for i in range(5):

    #    a , b = mod._sum_stat(bestfit , prior_range = None)
    #    n_ens.append(a)
    #    w_ens.append(b)
    #n_opt , w_opt = np.median(np.array(n_ens) , axis=0) , np.mean(np.array(w_ens) , axis=0)
    #model_obs = [n_opt , w_opt]

    if obs == "gmf":
  
       model_nbar , model_gmf = model_obs[0] , model_obs[1]
       data_nbar , data_gmf = data_obs[0] , data_obs[1]
       nbar_var , gmf_cov = cov_obs[0] , cov_obs[1]
       res_nbar = model_nbar - data_nbar
       res_gmf = model_gmf - data_gmf
       chisq_nbar = (res_nbar**2.)/(nbar_var)
       chisq_gmf = np.sum((res_gmf**2./gmf_cov))
       chisq = chisq_nbar + chisq_gmf
       ndata = 1. + len(res_gmf)      
 
    if obs == "wp":

       model_nbar , model_wp = mod._sum_stat(bestfit , prior_range = None)
       #model_nbar , model_wp = model_obs[0] , model_obs[1]
       data_nbar , data_wp = data_obs[0] , data_obs[1]
       nbar_var , wp_cov = cov_obs[0] , cov_obs[1]
       res_nbar = model_nbar - data_nbar
       res_wp = model_wp - data_wp
       f_bias = (400. - len(res_wp) -2.)/(400. - 1.)
       chisq_nbar = (res_nbar**2.)/(nbar_var)
       chisq_wp = f_bias * np.sum(np.dot(res_wp , solve(wp_cov , res_wp)))
       chisq = chisq_nbar + chisq_wp
       ndata = 1. + len(res_wp)

    aic1 = chisq + 2. * npar
    aic2 = chisq + (2. * npar * ndata)/(ndata - npar -1.)    
    bic = chisq + npar * np.log(ndata)
    print chisq , aic1 , aic2 , bic
 
    goodness = np.array([chisq, aic1, aic2, bic])
    np.savetxt("results/goodness_"+obs+"_"+model+"_"+str(Mr)+".dat" , goodness) 

    return None

def compute(filename1, filename2, nchains, nburnins, Mr, obs):

    import data as data_wp
    #import data_group as data_group

    fig = plt.figure(figsize=(6,6))
    #gs = gridspec.GridSpec(2,3)#, height_ratios=[1, 1], width_ratios=[1,1])  
    ax = fig.add_subplot(111) #Mr18
    
    sample_dec = h5py.File(filename1 , "r")["mcmc"][nburnins:nchains]
    sample_hod = h5py.File(filename2 , "r")["mcmc"][nburnins:nchains]

    sample_dec = sample_dec.reshape(sample_dec.shape[0]*sample_dec.shape[1] , sample_dec.shape[2])
    sample_hod = sample_hod.reshape(sample_hod.shape[0]*sample_hod.shape[1] , sample_hod.shape[2])

    bestfit_dec = np.median(sample_dec, axis = 0)
    bestfit_hod = np.median(sample_hod, axis = 0)
   
    bestfit_hod2 = bestfit_hod[:5]
 
    print bestfit_dec
    print bestfit_hod
    print bestfit_hod2
 
    if obs == "wp":

       data = data_wp.load_data(Mr)[1]
       cov = data_wp.load_covariance(Mr)[1]
       mod_dec = dec_wp(Mr)._sum_stat(bestfit_dec , prior_range = None)[1]
       mod_hod = hod_wp(Mr)._sum_stat(bestfit_hod , prior_range = None)[1]
       mod_hod2 = hod_wp(Mr)._sum_stat(bestfit_hod2 , prior_range = None)[1]
       #print mod_dec/data
       #print mod_hod/data 
       
    rbin = np.loadtxt("bin.dat")
    rbin = np.mean(rbin , axis = 1) 
    
    err = np.diag(cov)**0.5 / data
    pretty_colors = prettycolors()    
    #pretty_colors[1] for decorated model and 
    #pretty_colors[5] for standard hod models
    ax.errorbar(rbin , np.ones_like(rbin) , yerr = err , color = "k")    
    ax.plot(rbin , mod_dec/data, color = pretty_colors[1])
    ax.plot(rbin , mod_hod/data, color = "green") # pretty_colors[5])
    ax.plot(rbin , mod_hod2/data , color = "red")

    ax.set_xscale('log')
    ax.set_xlim([0.05, 30.])
    ax.set_ylim([0., 2.])
    plt.xticks(fontsize=20)
    plt.savefig("res.pdf")
       
    return None

 
if __name__=='__main__':

   directory = "/export/bbq2/mj/chains/"
   filename1 = "mcmc_chain_Mr19.0.hdf5"
   filename1 = directory+filename1
   filename2 = "adhoc_mcmc_chain_Mr19.0.hdf5"
   filename2 = directory+filename2
   Mr = 19.0
   obs , model = "wp" , "dec"
   nchains = 3000
   nburnins = 2500
   #compute(filename1, filename2, nchains, nburnins, Mr, obs)

   #compute_chisq_aic_bic(filename , nchains, nburnins, Mr, obs, model)
   #model_predictions(filename, Mr, nburnins, nchains, obs, model)   
   #plot_model_prediction(obs = "wp", model= "dec", clotter = True)
   plot_wpmodel_prediction(obs = "wp", clotter = True)
