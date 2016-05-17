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

def composite_model(Mr):

    cens_occ_model =  AssembiasZheng07Cens(threshold = -1.*Mr)
    cens_prof_model = TrivialPhaseSpace()
    sats_occ_model =  AssembiasZheng07Sats(threshold = -1.*Mr)
    sats_prof_model = NFWPhaseSpace()
    
    return HodModelFactory(
               centrals_occupation = cens_occ_model,
               centrals_profile = cens_prof_model,
               satellites_occupation = sats_occ_model,
               satellites_profile = sats_prof_model)

class MCMC_model(object):

    def __init__(self, Mr):
        
        self.Mr = Mr
        self.model = composite_model(Mr)
        self.halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')

        ###pair counter settings ###

        self.boxsize = self.halocat.Lbox
        self.nthreads = 1
        self.pimax = 40.0
        self.binfile = path.join(path.dirname(path.abspath(__file__)),
                        "../", "bin")
        self.autocorr = 1

    def __call__(self, theta, prior_range):
        return self._sum_stat(theta, prior_range=prior_range)

    def _sum_stat(self, theta, prior_range=None):
        
        self.model.param_dict['logM0'] = theta[0]
        self.model.param_dict['sigma_logM'] = theta[1]
        self.model.param_dict['logMmin'] = theta[2]
        self.model.param_dict['alpha'] = theta[3]
        self.model.param_dict['logM1'] = theta[4]
        self.model.param_dict['mean_occupation_centrals_assembias_param1']= theta[5]
        self.model.param_dict['mean_occupation_satellites_assembias_param1']= theta[6]

        self.model.populate_mock(self.halocat) 
        x = self.model.mock.galaxy_table['x']
        y = self.model.mock.galaxy_table['y']
        z = self.model.mock.galaxy_table['z']
        vz = self.model.mock.galaxy_table['vz']
        # applying RSD
        pos = return_xyz_formatted_array(x, y, z, velocity = vz, velocity_distortion_dimension = 'z')
        # enforcing PBC
        pos = enforce_periodicity_of_box(pos, self.boxsize)
        pos = pos.astype(np.float32)
        x, y, z = pos[:,0] , pos[:,1] , pos[:,2]
        results_wp = _countpairs.countpairs_wp(self.boxsize, self.pimax, 
                                           self.nthreads,
                                           self.binfile, x, y, z)
        wp = np.array(results_wp)[:,3]
        return wp

prior_min = [10., 0.05, 10., 0.85, 12., -1. , -1.]
prior_max = [14.5, 0.7, 14., 1.45, 15., 1. , 1.]  

def plot_time_mcmc(Nwalkers, Nchains, filename):

    sample = h5py.File(filename , "r")["mcmc"]
    npars = sample.shape[2]
    fig , axes = plt.subplots(npars , 1 , sharex=True, figsize=(10, 12))

    for i in xrange(npars):
        axes[i].plot(sample[:, :, i], color="b", alpha=.4 , lw = .5)
	axes[i].yaxis.set_major_locator(MaxNLocator(5))
        #axes[i].set_ylim([plot_range[i,0], plot_range[i,1]])
        axes[i].set_xlim(0, Nchains)
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
            plot_datapoints=False,
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

def plot_predictions(Mr, nburnins, assembly = True):


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
    sample = sample[-nburnins:, : , :]
    sample = sample.reshape(sample.shape[0]*sample.shape[1] , sample.shape[2])
    model = MCMC_model(Mr) 

    model_wp = []
    for i in xrange(len(sample)):
        print i
        model_wp.append(model._sum_stat(sample[i] , prior_range = None))
    np.savetxt("wp_assembly_"+str(Mr)+".dat" , np.array(model_wp)) 
 
    return None 

if __name__=='__main__':

   filename = "Mr20.0.hdf5"
   #filename = "mcmc_chain_Mr20.0.hdf5"
   #filename = "../../dat/mcmc/mcmc_chain_Mr20.0.hdf5"
   #plot_time_mcmc(Nwalkers = 24, Nchains = 30000, filename=filename)
   plot_predictions(20.0 , 1000 , True)
   #plot_corner_mcmc(Nchains = 21000 , Nburns=17000, Mr = 20, filename=filename)
