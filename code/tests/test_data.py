'''
testing how the model fits the data
'''
from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
import matplotlib.pyplot as plt
import os.path as path
import time
from Corrfunc import _countpairs
from Corrfunc.utils import read_catalog
# --- Local ---
# --- halotools ---
from halotools.sim_manager import CachedHaloCatalog
from halotools.empirical_models import HodModelFactory
from halotools.empirical_models import TrivialPhaseSpace, Zheng07Cens
from halotools.empirical_models import Zheng07Sats , HeavisideAssembias
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
from halotools.empirical_models import enforce_periodicity_of_box
prettyplot()

class AssembiasZheng07Sats(Zheng07Sats, HeavisideAssembias):

    def __init__(self, **kwargs):

        Zheng07Sats.__init__(self, threshold = -21)

        HeavisideAssembias.__init__(self,
            method_name_to_decorate = 'mean_occupation',
            lower_assembias_bound = 0.,
            upper_assembias_bound = np.inf,
            **kwargs)

cens_occ_model =  Zheng07Cens(threshold = -21)
cens_prof_model = TrivialPhaseSpace()
sats_occ_model =  AssembiasZheng07Sats()
sats_prof_model = NFWPhaseSpace()

model= HodModelFactory(
        centrals_occupation = cens_occ_model,
        centrals_profile = cens_prof_model,
        satellites_occupation = sats_occ_model,
        satellites_profile = sats_prof_model)
def main():
     
    cov = np.loadtxt("../../data/wpxicov_dr72_bright0_mr21.0_z0.159_nj400")
    f_MD = (1. + 71.74*10**6. / (1000.)**3.)
    f_bol = (1. + 71.74*10.**6. / (250.)**3.)
      
    print("covariance correction factor=" , f_bol/f_MD)
    cov = cov*f_bol/f_MD
   
    model.param_dict['logM0'] =  12.59
    model.param_dict['sigma_logM'] =  0.49
    model.param_dict['logMmin'] =  12.78
    model.param_dict['alpha'] =  1.14
    model.param_dict['logM1'] =  13.99

    model.param_dict['mean_occupation_satellites_assembias_param1'] = 0.0
    halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')
    model.populate_mock(halocat, enforce_PBC = True)
    pos = three_dim_pos_bundle(model.mock.galaxy_table, 'x', 'y', 'z')
    
    x = model.mock.galaxy_table['x']
    y = model.mock.galaxy_table['y']
    z = model.mock.galaxy_table['z']
    vz = model.mock.galaxy_table['vz']

    # applying RSD
    pos = return_xyz_formatted_array(x, y, z, velocity = vz, velocity_distortion_dimension = 'z')
    # enfprcing PBC
    pos =  enforce_periodicity_of_box(pos, halocat.Lbox)

    tstart = time.time()
    t0 = tstart
    pos = pos.astype(np.float32)
    x, y, z = pos[:,0] , pos[:,1] , pos[:,2]
    t1 = time.time()
    print("Done reading the data - time taken = {0:10.1f} seconds"
          .format(t1 - t0))
    print("Beginning Correlation functions calculations")
    boxsize = halocat.Lbox
    nthreads = 4
    pimax = 40.0
    binfile = path.join(path.dirname(path.abspath(__file__)),
                        "../../", "bin")
    autocorr = 1
    numbins_to_print = 12
    
    print("\nRunning 2-D projected correlation function wp(rp)")
    results_wp = _countpairs.countpairs_wp(boxsize, pimax, nthreads,
                                           binfile, x, y, z)
    print("\n#            ******    wp: first {0} bins  *******         "
          .format(numbins_to_print))
    print("#      rmin        rmax       rpavg        wp       npairs")
    print("##########################################################")
    for ibin in range(numbins_to_print):
        items = results_wp[ibin]
        print("{0:12.4f} {1:12.4f} {2:10.4f} {3:10.1f} {4:10d}"
              .format(items[0], items[1], items[2], items[3], items[4]))
    print("-----------------------------------------------------------")
    
    data_wp = np.loadtxt("../../data/wpxi_dr72_bright0_mr21.0_z0.159_nj400")[:,1]
    print(data_wp.shape)
    data_wp_error = np.sqrt(np.diag(cov)[:12])
    print(data_wp_error.shape) 
    rbins = np.loadtxt(binfile)
    rs = np.mean(rbins , axis = 1)
    plt.figure(figsize=(10,10))
    plt.errorbar(rs , data_wp , data_wp_error , fmt=".k" , capsize = 2)
    plt.plot(rs , np.array(results_wp)[:,3])
    plt.loglog()
    plt.savefig("wp.pdf")
    
 
if __name__ == "__main__":
    main()
