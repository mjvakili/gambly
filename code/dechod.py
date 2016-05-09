'''
Model Evaluation
'''
import numpy as np
import os.path as path

from Corrfunc import _countpairs
from Corrfunc.utils import read_catalog

from halotools.sim_manager import CachedHaloCatalog
from halotools.empirical_models import HodModelFactory
from halotools.empirical_models import TrivialPhaseSpace, Zheng07Cens
from halotools.empirical_models import Zheng07Sats , HeavisideAssembias
from halotools.empirical_models import PrebuiltHodModelFactory
from halotools.empirical_models.factories.mock_helpers import three_dim_pos_bundle
from halotools.mock_observables.catalog_analysis_helpers import return_xyz_formatted_array
from halotools.empirical_models import NFWPhaseSpace
from halotools.empirical_models import enforce_periodicity_of_box

class AssembiasZheng07Sats(Zheng07Sats, HeavisideAssembias):
   
    def __init__(self, Mr = 21 , **kwargs):
        
        self.Mr = Mr
        Zheng07Sats.__init__(self, threshold = -1.*self.Mr)

        HeavisideAssembias.__init__(self,
            method_name_to_decorate = 'mean_occupation',
            lower_assembias_bound = 0.,
            upper_assembias_bound = np.inf,
            **kwargs)

def composite_model(Mr):

    cens_occ_model =  Zheng07Cens(threshold = -1.*Mr)
    cens_prof_model = TrivialPhaseSpace()
    sats_occ_model =  AssembiasZheng07Sats(Mr = 21)
    sats_prof_model = NFWPhaseSpace()
    
    return HodModelFactory(
               centrals_occupation = cens_occ_model,
               centrals_profile = cens_prof_model,
               satellites_occupation = sats_occ_model,
               satellites_profile = sats_prof_model)


class MCMC_model(object):

    def __init__(self, Mr=21):
        
        self.Mr = Mr
        thr = -1. * np.float(Mr)
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
        self.model.param_dict['mean_occupation_satellites_assembias_param1']= theta[5]
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
        nbar = 1.*len(pos)/(self.boxsize)**3.

        return nbar , wp
