
from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import os.path as path
import time
from Corrfunc import _countpairs
from Corrfunc.utils import read_catalog
import numpy as np 
# --- Local ---
# --- halotools ---
from halotools.sim_manager import CachedHaloCatalog
from halotools.empirical_models import PrebuiltHodModelFactory
from halotools.mock_observables import tpcf
from halotools.empirical_models.factories.mock_helpers import three_dim_pos_bundle
from halotools.mock_observables import FoFGroups
from halotools.mock_observables.pair_counters import npairs_3d

import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors 
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors
from Corrfunc.utils import read_catalog

def main():
     
    # Entire MultiDark Volume (Analytic xi) 
    cov = np.loadtxt("../data/wpxicov_dr72_bright0_mr21.0_z0.159_nj400")
    print(cov.shape)
    model = PrebuiltHodModelFactory('zheng07', threshold=-21)
    print(model.param_dict)
    model.param_dict['logM0'] =  12.59
    model.param_dict['sigma_logM'] =  0.49
    model.param_dict['logMmin'] =  12.78
    model.param_dict['alpha'] =  1.14
    model.param_dict['logM1'] =  13.99
    #, 'sigma_logM': 0.39, 'logMmin': 12.79, 'alpha': 1.15, 'logM1': 13.94}

    halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')
    model.populate_mock(halocat, enforce_PBC = True)
    pos = three_dim_pos_bundle(model.mock.galaxy_table, 'x', 'y', 'z')

    tstart = time.time()
    t0 = tstart
    pos = pos.astype(np.float32)
    x, y, z = pos[:,0] , pos[:,1] , pos[:,2]
    t1 = time.time()
    print("Done reading the data - time taken = {0:10.1f} seconds"
          .format(t1 - t0))
    print("Beginning Correlation functions calculations")
    boxsize = 250
    nthreads = 4
    pimax = 40.0
    binfile = path.join(path.dirname(path.abspath(__file__)),
                        "../", "bin")
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
    
    data_wp = np.loadtxt("../data/wpxi_dr72_bright0_mr21.0_z0.159_nj400")[:,1]
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
