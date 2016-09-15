import util as ut
import data
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



def plot_wprandom(obs = "wp", clotter = True):

    if clotter == False:

        print "running compute_model_prediction first"

    else:

        #### loading the model predictions for the observables ######         

        file20 = np.loadtxt("../dat/wp_AGEMATCHING_Mr-20.0.dat")
        file205 = np.loadtxt("../dat/wp_AGEMATCHING_Mr-20.5.dat")
        file21 = np.loadtxt("../dat/wp_AGEMATCHING_Mr-21.0.dat")

        rand20 = np.loadtxt("../dat/wp_AGEMATCHING_randomized_Mr-20.0.dat")
        rand205 = np.loadtxt("../dat/wp_AGEMATCHING_randomized_Mr-20.5.dat")
        rand21 = np.loadtxt("../dat/wp_AGEMATCHING_randomized_Mr-21.0.dat")
        #### loading the observables in Mr 18,19,20


        if obs == "wp":

            rbin20 = np.loadtxt("bin.dat")
            rbin205 = np.loadtxt("bin.dat")
            rbin21 = np.loadtxt("bin.dat")

            rbin20 = np.mean(rbin20, axis = 1)             
            rbin205 = np.mean(rbin205, axis = 1)             
            rbin21 = np.mean(rbin21, axis = 1)             

	    wbin20 = 1.
	    wbin205 = 1.
	    wbin21 = 1.
            
            data20 = data.load_wp(20.)
            data205 = data.load_wp(20.5)
            data21 = data.load_wp(21.0)
            
            err20 = np.diag(data.load_wp_covariance(20.))**.5 
            err205 = np.diag(data.load_wp_covariance(20.5))**.5 
            err21 = np.diag(data.load_wp_covariance(21.0))**.5 

	#prettyplot()
    	pretty_colors=prettycolors()
    	fig = plt.figure(1, figsize=(21,7))
    	gs = gridspec.GridSpec(1,3)#, height_ratios=[1, 1], width_ratios=[1,1])  

        ax = plt.subplot(gs[0,0]) #Mr20

        ax.plot(rbin20, rand20/file20 - 1 , color='#ee6a50', alpha=1.0 , lw = 3) 
        wp = np.loadtxt("results/dec_random_20.0.dat") 
        nsamples = wp.shape[0]/2
        mock , random = wp[:nsamples,:] , wp[nsamples:,:]

        red_line = mlines.Line2D([], [], ls = '-', c = '#ee6a50', linewidth=3, 
                           label=r'HW+13 abundance matching')
        blue_line = mlines.Line2D([], [], ls = '-', c = 'blue', alpha = 1.0, 
                           label=r'HOD with assembly bias, this work')

        
        ax.errorbar(rbin20, np.zeros_like(data20), yerr=err20/data20, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin20, np.zeros_like(data20), c='k', s=10, lw=0)
        for i in range(nsamples):
            ax.plot(rbin20 , random[i,:]/mock[i,:] - 1. , alpha = 0.1 , color = 'blue')

        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          ax.set_ylabel(r'$w_{p}^{\mathrm{ranomized}}/w_{p}^{\mathrm{mock}} \; -1 $', fontsize=27)
          #ax.set_yscale('log') 
          ax.set_xscale('log')
          #ax.set_xticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([-0.5, 0.5])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(2.0, -0.4, r'$M_{r}<-20$', fontsize=25) 
          plt.legend(handles=[blue_line, red_line], frameon=False, loc='best', fontsize=15)
        ax = plt.subplot(gs[0,1]) #Mr20.5
        ax.plot(rbin20, rand205/file205 - 1 , color='#ee6a50', alpha=1.0 , lw = 3) 

        wp = np.loadtxt("results/dec_random_20.5.dat") 
        nsamples = wp.shape[0]/2
        mock , random = wp[:nsamples,:] , wp[nsamples:,:]
        for i in range(nsamples):
            ax.plot(rbin20 , random[i,:]/mock[i,:] - 1. , alpha = 0.1 , color = 'blue')
        
        ax.errorbar(rbin205, np.zeros_like(data205), yerr=err205/data205, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin20, np.zeros_like(data20), c='k', s=10, lw=0)
        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          ax.set_xscale('log')
          #ax.set_xticklabels([])
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([-0.5, 0.5])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(2.0, -0.4, r'$M_{r}<-20.5$', fontsize=25) 
        
        ax = plt.subplot(gs[0,2]) #Mr21.
        ax.plot(rbin20, rand21/file21 - 1 , color='#ee6a50', alpha=1.0 , lw = 3) 

        ax.errorbar(rbin21, np.zeros_like(data21), yerr=err21/data21, fmt="o", color='k', markersize=0, lw=0, capsize=3, elinewidth=1.5)
        ax.scatter(rbin21, np.zeros_like(data21), c='k', s=10, lw=0)
        wp = np.loadtxt("results/dec_random_21.0.dat") 
        nsamples = wp.shape[0]/2
        mock , random = wp[:nsamples,:] , wp[nsamples:,:]
        for i in range(nsamples):
            ax.plot(rbin20 , random[i,:]/mock[i,:] - 1. , alpha = 0.1 , color = 'blue')

        if obs == "wp":
          ax.set_xlabel(r'$r_{p} \; [\mathrm{Mpc}\; h^{-1}]$', fontsize=27)
          #ax.set_ylabel(r'$w_{p}(r_{p}) \; [\mathrm{Mpc} \; h^{-1}]$', fontsize=27)
          #ax.set_yscale('log') 
          ax.set_xscale('log')
          ax.set_yticklabels([])
          ax.set_xlim([0.05, 30.])
          ax.set_ylim([-0.5, 0.5])
          plt.xticks(fontsize=20)
          plt.yticks(fontsize=20)
          ax.text(2.0, -0.4, r'$M_{r}<-21$', fontsize=25) 

        
        fig.subplots_adjust(wspace=0.0, hspace=0.0)
        fig_name = ''.join([ut.fig_dir(), 
        'paper', 
        '.wprandom',
        '.pdf'])
        
    fig.savefig(fig_name, bbox_inches='tight')
    plt.close()
    return None 


if __name__=='__main__':

    plot_wprandom(obs = "wp", clotter = True)
