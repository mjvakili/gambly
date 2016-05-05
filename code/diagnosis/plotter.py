'''

Plotting modules

'''
import os
import h5py
import corner
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator

plt.switch_backend("Agg")


def plot_mcmc(Nwalkers, Nburns, filename):

    sample = h5py.File(filename , "r")["mcmc"]
    npars = sample.shape[2]
    fig , axes = plt.subplots(npars , 1 , sharex=True, figsize=(10, 12))

    for i in xrange(npars):
        axes[i].plot(sample[:, :, i], color="b", alpha=.4 , lw = .5)
	axes[i].yaxis.set_major_locator(MaxNLocator(5))
        #axes[i].set_ylim([plot_range[i,0], plot_range[i,1]])
        axes[i].set_xlim(0, 7750)
        #axes[i].set_ylabel(labels[i], fontsize=25)

    axes[-1].set_xlabel("Step Number", fontsize=25)
    fig.tight_layout(h_pad=0.0)
    fig_file = "mcmc_time.pdf"
    plt.savefig(fig_file)
    plt.close()

if __name__=='__main__':
   plot_mcmc(Nwalkers=24, Nburns=3000, filename = "../../dat/mcmc/mcmc_chain2.hdf5")
