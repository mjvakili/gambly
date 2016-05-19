'''
testing how the model fits the data
'''
from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import os.path as path
import time
from Corrfunc import _countpairs
from Corrfunc.utils import read_catalog
# --- Local ---
# --- halotools ---
from halotools.empirical_models import TrivialPhaseSpace, AssembiasZheng07Cens
from halotools.empirical_models import NFWPhaseSpace, AssembiasZheng07Sats
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

        Zheng07Sats.__init__(self, threshold = -20)

        HeavisideAssembias.__init__(self,
            method_name_to_decorate = 'mean_occupation',
            lower_assembias_bound = 0.,
            upper_assembias_bound = np.inf,
            **kwargs)

cens_occ_model =  AssembiasZheng07Cens(threshold = -20.)
cens_prof_model = TrivialPhaseSpace()
sats_occ_model =  AssembiasZheng07Sats(threshold = -20.)
sats_prof_model = NFWPhaseSpace()

model= HodModelFactory(
        centrals_occupation = cens_occ_model,
        centrals_profile = cens_prof_model,
        satellites_occupation = sats_occ_model,
        satellites_profile = sats_prof_model)

def richness(group_id):
    gals = Table()
    gals['groupid'] = group_id
    gals['dummy'] = 1
    grouped_table = gals.group_by('groupid')
    grp_richness = grouped_table['dummy'].groups.aggregate(np.sum)
    return grp_richness

def main():
     
    #cov = np.loadtxt("../../data/wpxicov_dr72_bright0_mr21.0_z0.159_nj400")
    #f_MD = (1. + 71.74*10**6. / (1000.)**3.)
    #f_bol = (1. + 71.74*10.**6. / (250.)**3.)
      
    #print("covariance correction factor=" , f_bol/f_MD)
    #cov = cov*f_bol/f_MD


    gmf_cat = np.loadtxt("gmf20.dat")
    data_gmf = gmf_cat[:,2]
    data_gmf_err = (gmf_cat[:,3]**2. + gmf_cat[:,4]**2.)**.5
    data_gmf_bin = np.hstack([gmf_cat[:,0],gmf_cat[-1,1]])
    print(data_gmf_bin) 
    model.param_dict['logM0'] =  11.16
    model.param_dict['sigma_logM'] =  0.61
    model.param_dict['logMmin'] =  12.10	
    model.param_dict['alpha'] =  1.08
    model.param_dict['logM1'] =  13.33
    model.param_dict['mean_occupation_centrals_assembias_param1'] = -1.
    
    halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')
    model.populate_mock(halocat)
    pos = three_dim_pos_bundle(model.mock.galaxy_table, 'x', 'y', 'z')
  
    print("modelnumber density=" , len(pos)/halocat.Lbox**3.)
    print("data number density=" , 1.16*10**-3.) 

    x = model.mock.galaxy_table['x']
    y = model.mock.galaxy_table['y']
    z = model.mock.galaxy_table['z']
    vz = model.mock.galaxy_table['vz']

    # applying RSD
    pos = return_xyz_formatted_array(x, y, z, velocity = vz, velocity_distortion_dimension = 'z')
    # enfprcing PBC
    #pos =  enforce_periodicity_of_box(pos, halocat.Lbox)

    from halotools.mock_observables import FoFGroups

    bperp = 0.14
    bpar = 0.75
    Lbox = np.array([halocat.Lbox,halocat.Lbox,halocat.Lbox])
    period = Lbox
    groups = FoFGroups(pos, b_perp=bperp, b_para=bpar, Lbox = Lbox , period=Lbox) 
  
    #model.mock.galaxy_table['fof_group_id'] = groups.group_ids
    #model.mock.galaxy_table.sort(['fof_group_id'])
    #from halotools.utils import group_member_generator


    #grouping_key = 'fof_group_id'
    #requested_columns = []

    #group_gen = group_member_generator(model.mock.galaxy_table, grouping_key, requested_columns)

    #group_richness = np.zeros(len(model.mock.galaxy_table), dtype=int)
    #for first, last, member_props in group_gen:
    #    group_richness[first:last] = last-first
    
    group_ids = groups.group_ids
    group_richness = richness(group_ids)        
    gmf = np.histogram(np.array(group_richness) , data_gmf_bin)[0] / halocat.Lbox**3.
   
    chi = np.sum(((gmf - data_gmf)**2. / data_gmf_err**2.)[:-1])
    print(chi/(len(gmf)-1)) 

    bin_centers = 0.5 * (data_gmf_bin[1:] + data_gmf_bin[:-1])

    plt.figure(figsize=(10,10))
    plt.errorbar(bin_centers , data_gmf , data_gmf_err , fmt=".k" , capsize = 2)
    plt.plot(bin_centers , gmf)
    plt.yscale("log")
    plt.xlim([0,60])
    plt.savefig("gmf.pdf")
    
 
if __name__ == "__main__":
    main()
