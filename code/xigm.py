import numpy as np
from matplotlib import lines as mlines
from matplotlib import pyplot as plt
import numpy as np
from halotools.empirical_models import PrebuiltHodModelFactory
from scipy import interpolate 

###Constants necessary in calculation of Sigma(R)###

Omegam = 0.315
rho_crit = 2.77536627 * 1.e11   #units h^2M_sun/Mpc^-3
r_los_max = 50
r_los_min = .001

def xigm_tabulated_generator(hod = 'leauthaud11' , tr = 9.8 , \
                             rmin = 0.001 , rmax = 80 , nbins = 10 , \
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

    model = PrebuiltHodModelFactory(hod ,threshold = tr , redshift = z)
    model.populate_mock(simname = sim , redshift = z)


    rr , xigm = model.mock.compute_galaxy_matter_cross_clustering(rbins = rbin)

    return rr , xigm


def xigm_func(r , r1 , r2):
 
    tabulated_rr , tabulated_xigm = xigm_tabulated_generator(rmin = r1 , rmax = r2)
   
    if r.min() < tabulated_rr.min():
        if r.max() > tabulated_rr.max():
    
            raise ValueError("invalid input for radius")    
 
    f = interpolate.interp1d(tabulated_rr, tabulated_xigm , kind = 'cubic')
    
    return f(r)
    


if  __name__=='__main__':

    import pylab
    import time
    a = time.time()
    rr ,xigm = xigm_tabulated_generator(rmin = .1 , rmax = 10. , nbins = 50)
    print time.time() - a
    pylab.loglog(rr , xigm)
    #pylab.show()
    rr2 = np.logspace(np.log10(.2) , np.log10(8) , 10)
    a = time.time()
    xigm2 = xigm_func(rr2, .1 , 10)
    print time.time() - a
    pylab.loglog(rr2 , xigm2)
    pylab.show()
