'''
model without galaxy assembly bias
'''
import numpy as np
import os.path as path

from Corrfunc import _countpairs
from Corrfunc.utils import read_catalog

from astropy.table import Table
from halotools.mock_observables import FoFGroups
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

def richness(group_id):
    gals = Table()
    gals['groupid'] = group_id
    gals['dummy'] = 1
    grouped_table = gals.group_by('groupid')
    grp_richness = grouped_table['dummy'].groups.aggregate(np.sum)
    return grp_richness

class MCMC_model(object):

    def __init__(self, Mr):
        
        self.Mr = Mr
        self.model = single_model(Mr)
        self.halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')

        #GMF binning settings

        self.boxsize = self.halocat.Lbox
        if self.Mr == 18:
           gmf_cat = np.loadtxt("../dat/gmf_mr18.0.dat")
        if self.Mr == 19:
           gmf_cat = np.loadtxt("../dat/gmf_mr19.0.dat")
        if self.Mr == 20:
           gmf_cat = np.loadtxt("../dat/gmf_mr20.0.dat")

        self.data_gmf_bin = np.hstack([gmf_cat[:,0],gmf_cat[-1,1]])
        self.data_gmf_bin_width=(self.data_gmf_bin[1:]-self.data_gmf_bin[:-1])        
    def __call__(self, theta, prior_range):
        return self._sum_stat(theta, prior_range=prior_range)

    def _sum_stat(self, theta, prior_range=None):
        
        self.model.param_dict['logM0'] = theta[0]
        self.model.param_dict['sigma_logM'] = theta[1]
        self.model.param_dict['logMmin'] = theta[2]
        self.model.param_dict['alpha'] = theta[3]
        self.model.param_dict['logM1'] = theta[4]
        gmff = []
        for i in xrange(1):
          self.model.populate_mock(self.halocat) 
          x = self.model.mock.galaxy_table['x']
          y = self.model.mock.galaxy_table['y']
          z = self.model.mock.galaxy_table['z']
          vz = self.model.mock.galaxy_table['vz']
          # applying RSD
          pos = return_xyz_formatted_array(x, y, z, velocity = vz, velocity_distortion_dimension = 'z')
          # enforcing PBC
          pos = enforce_periodicity_of_box(pos, self.boxsize)
          bperp = 0.14
          bpar = 0.75
          Lbox = np.array([self.boxsize, self.boxsize, self.boxsize])
          period = Lbox
          groups = FoFGroups(pos, b_perp=bperp, b_para=bpar, Lbox = Lbox , period=Lbox) 
  
          group_ids = groups.group_ids
          group_richness = richness(group_ids)
          gmff.append(np.histogram(np.array(group_richness) ,self.data_gmf_bin)[0] / (self.data_gmf_bin_width * self.boxsize**3.))
        gmf = np.mean(np.array(gmff) , axis = 0)
        nbar = 1.*len(pos)/(self.boxsize)**3.

        return nbar , gmf

