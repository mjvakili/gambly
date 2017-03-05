'''
Central/Satellite galaxy assembly biased model
'''
import numpy as np
import os.path as path

from Corrfunc import _countpairs
from Corrfunc.utils import read_catalog

from halotools.sim_manager import CachedHaloCatalog
from halotools.empirical_models import HodModelFactory
from halotools.empirical_models import TrivialPhaseSpace, AssembiasZheng07Cens
from halotools.empirical_models import NFWPhaseSpace, AssembiasZheng07Sats
from halotools.empirical_models import PrebuiltHodModelFactory
from halotools.empirical_models.factories.mock_helpers import three_dim_pos_bundle
from halotools.mock_observables.catalog_analysis_helpers import return_xyz_formatted_array
from halotools.empirical_models import enforce_periodicity_of_box

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
        self.rbins = np.loadtxt(self.binfile) 
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
        
        ########################## rewrite for model covariance ######################
        
        #applying RSD
        pos = return_xyz_formatted_array(x, y, z, velocity = vz, velocity_distortion_dimension = 'z')
        # enforcing PBC
        pos = enforce_periodicity_of_box(pos, self.boxsize)
        pos = pos.astype(np.float32)
        #x, y, z = pos[:,0] , pos[:,1] , pos[:,2]

        wp , wp_cov = helpers.projected_correlation(pos, self.rbins, self.pimax, self.boxsize, 
                                           jackknife_nside = 8 , bias_correction = True) 

        nbar = 1.*len(pos)/(self.boxsize)**3.
        
        nsub = 8. #FIXME
        number_of_subboxes = nsub ** 3
        subbox_size = self.boxsize / nsub
        
        nbars = np.zeros((number_of_subboxes , 1)) #placeholder for nbars of jk subsamples

        for subvol_index in xrange(number_of_subboxes):

            sub_rsd_pos = mask_positions(rsd_pos , subvol_index , nsub)
            sub_rsd_pos = sub_rsd_pos.astype(np.float32)
            nbars[subvol_index] = 1.*len(sub_rsd_pos)/subbox_size**3.
           
        nbar_var = np.array([np.var(nbars)])    
        
        ######################## End of rewriting for model covariance ############## 
        
        """
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
        nbar = 1.*len(pos)/(self.boxsize)**3.
        """
        return nbar , wp , nbar_var , wp_cov

def mask_positions(pos , boxdize, subvol_index , nsub):

    '''masks the positions of galaxies in
       model to compute jk covariance'''

    
    subbox_size = 1.*box_size / nsub
    
    xi , yi , zi  = edge(subvol_index, nsub)
    submask = np.where((xi <pos[:, 0]) * \
                       (pos[:, 0] < xi + subbox_size) * \
                       (yi <pos[:, 1]) * \
                       (pos[:, 1] < yi + subbox_size) * \
                       (zi <pos[:, 2]) *  \
                       (pos[:, 2] < zi + subbox_size))
    
    return pos[submask] - edge(subvol_index,nsub)
