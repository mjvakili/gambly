import numpy as np
from matplotlib import lines as mlines
from matplotlib import pyplot as plt
import numpy as np
from halotools.empirical_models import PrebuiltHodModelFactory


def xigm_tabulated_generator(hod = 'leauthaud11' , threshold = 10.6 , \
                             rmin = 0.001 , rmax = 80 , nbins = 100 , \
                             scale = 'log' , sim = 'bolshoi', z = 0.0):


    """
    returns galaxy_matter_cross_correlation <delta_g  delta_m >
    
    arguments:
              hod : HOD model, default = Leauthaud
  	      threshold : log10 of stellar mass threshold in units of Msun/h, default = 10.6
              rmin : minimum radius of the r-bins in Mpc/h, default = .001 Mpc/h
              rmax : maximum radius of the r-bins in Mpc/h, default = 50 Mpc/h
              nbins : number of rbins, default = 100
              scale : log or linear, scaling of the rbins, default = log
              simulation name : bolshoi or multidark, default = bolshoi
              z = redshift

    """
    
    if z<0:
        
        raise ValueError("redshift cannot be negative")

    if scale == 'log':

        rbin = np.logspace(np.log10(rmin) , np.log10(rmax) , nbins)

    elif scale == 'linear':

        rbin = np.linspace(rmin , rmax , nbins)

    else:

        raise ValueError("invalid scaling for the rbins")

    model = PrebuiltHodModelFactory(hod , threshold)
    model.populate_mock(simname = sim , redshift = z)


    rr , xigm = model.mock.compute_galaxy_matter_cross_clustering(rbins = rbin)

    return rr , xigm
