'''

Plotting modules

'''
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

prior_min = [10., 0.05, 10., 0.85, 12., -1. , -1.]
prior_max = [14.5, 0.7, 14., 1.45, 15., 1. , 1.]  

def plot_time_mcmc(Nwalkers, Nburns, filename):

    sample = h5py.File(filename , "r")["mcmc"]
    npars = sample.shape[2]
    fig , axes = plt.subplots(npars , 1 , sharex=True, figsize=(10, 12))

    for i in xrange(npars):
        axes[i].plot(sample[:, :, i], color="b", alpha=.4 , lw = .5)
	axes[i].yaxis.set_major_locator(MaxNLocator(5))
        #axes[i].set_ylim([plot_range[i,0], plot_range[i,1]])
        axes[i].set_xlim(0, 20000)
        #axes[i].set_ylabel(labels[i], fontsize=25)

    axes[-1].set_xlabel("Step Number", fontsize=25)
    fig.tight_layout(h_pad=0.0)
    fig_file = "mcmc_time.pdf"
    plt.savefig(fig_file)
    plt.close()

def plot_corner_mcmc(Nchains , Nburns, Mr , filename):

    sample = h5py.File(filename , "r")["mcmc"]
    npars = sample.shape[2]
    nwalkers = sample.shape[1]
    sample = sample[Nburns:Nchains, : , :]
    sample = sample.reshape(sample.shape[0]*sample.shape[1] , sample.shape[2])
    print sample.shape    
    #prior_min, prior_max = PriorRange('first_try' , Mr)
    prior_range = np.zeros((len(prior_min),2))
    prior_range[:,0] = np.array(prior_min)
    prior_range[:,1] = np.array(prior_max) 

    fig = corner.corner(
            sample,
            #labels=[
            #    r'$\log\;\mathcal{M}_{0}}$',
            #    r'$\log\;\sigma_\mathtt{log\;M}}$',
            #    r'$\log\;\mathcal{M}_\mathtt{min}}$',
            #    r'$\alpha$',
            #    r'$\log\;\mathcal{M}_{1}}$'
            #    r'$\mathcal{A}_{1}}$'
            #    r'$\mathcal{A}_{2}}$'
            #    ],
            #label_kwargs={'fontsize': 25},
            range=prior_range,
            quantiles=[0.16,0.5,0.84],
            show_titles=True,
            title_args={"fontsize": 12},
            plot_datapoints=True,
            fill_contours=True,
            levels=[0.68, 0.95],
            color='#ee6a50',
            bins=20,
            smooth=1.0)
    fig_name = ''.join(['post',
         str(Mr),  
        '.pdf'])
    fig.savefig(fig_name, bbox_inches='tight', dpi=150) 
    plt.close()
    return None 

if __name__=='__main__':
   filename = "mcmc_chain_Mr20.0.hdf5"
   #filename = "../../dat/mcmc/mcmc_chain_Mr20.0.hdf5"
   #plot_mcmc(Nwalkers=24, Nburns=3000, filename=filename)
   plot_corner_mcmc(Nchains = 12000 , Nburns=6000, Mr = 20, filename=filename)
