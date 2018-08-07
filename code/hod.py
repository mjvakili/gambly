'''
Central/Satellite galaxy without assembly biased model
'''
import numpy as np
import os.path as path
from halo_utils import load_project_halocat 
from Corrfunc.theory.wp import wp

from halotools.sim_manager import CachedHaloCatalog
from halotools.empirical_models import HodModelFactory
from halotools.empirical_models import TrivialPhaseSpace, AssembiasZheng07Cens
from halotools.empirical_models import NFWPhaseSpace, AssembiasZheng07Sats
from halotools.empirical_models import PrebuiltHodModelFactory
from halotools.empirical_models.factories.mock_helpers import three_dim_pos_bundle
from halotools.mock_observables.catalog_analysis_helpers import return_xyz_formatted_array
from halotools.empirical_models import enforce_periodicity_of_box

def single_model(Mr):

    model = PrebuiltHodModelFactory("zheng07" , threshold = -1.*Mr)

    return model

class MCMC_model(object):

    def __init__(self, Mr, box, halocat):
        
        self.Mr = Mr
        self.model = single_model(Mr)
        #self.halocat = load_project_halocat(box)
        self.halocat = halocat
        self.boxsize = self.halocat.Lbox[0]
        self.nthreads = 1
        self.pimax = 40.0
        self.binfile = path.join(path.dirname(path.abspath(__file__)),
                        "../", "bin")
        self.rbins = np.loadtxt("rbins.dat") 
        self.autocorr = 1

    def __call__(self, theta, prior_range):
        return self._sum_stat(theta, prior_range=prior_range)

    def _sum_stat(self, theta, prior_range=None):
        
        self.model.param_dict['logM0'] = theta[0]
        self.model.param_dict['sigma_logM'] = theta[1]
        self.model.param_dict['logMmin'] = theta[2]
        self.model.param_dict['alpha'] = theta[3]
        self.model.param_dict['logM1'] = theta[4]

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
	#print pos
        x, y, z = pos[:,0], pos[:,1], pos[:,2] 

        wp_result = wp(self.boxsize, self.pimax, 
                       self.nthreads, self.rbins, 
		       x, y, z)
        wp_result = np.array([wp_bin[3] for wp_bin in wp_result])
        
        nbar = 1.*len(pos)/(self.boxsize)**3.
        
	return nbar , wp_result

def edge(boxsize, index , nsub):
    '''returns edges of a sub-box of 
       a given index
    '''
    subbox_size = 1.*boxsize / nsub

    zi = (index / (nsub**2)) * subbox_size
    i2 = index % (nsub**2)
    yi = (i2 / nsub) * subbox_size
    i3 = i2 % nsub
    xi = (i3) * subbox_size

    return xi , yi , zi

def mask_positions(pos , boxsize, subvol_index , nsub):

    '''masks the positions of galaxies in
       model to compute jk covariance'''

    
    subbox_size = 1.*boxsize / nsub
    
    xi , yi , zi  = edge(boxsize, subvol_index, nsub)

    submask = np.where((xi <pos[:, 0]) * \
                       (pos[:, 0] < xi + subbox_size) * \
                       (yi <pos[:, 1]) * \
                       (pos[:, 1] < yi + subbox_size) * \
                       (zi <pos[:, 2]) *  \
                       (pos[:, 2] < zi + subbox_size))
    #mask  = np.where(np.arange(len(pos))!=submask) 
    return pos[submask,:][0]
