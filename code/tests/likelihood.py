import wp_model
import numpy as np
from scipy import interpolate 
import pylab
"""
{'alphasat': 1.0,
 'bcut': 1.47,
 'betacut': -0.13,
 'betasat': 0.859,
 'bsat': 10.62,
 'mean_occupation_centrals_assembias_param1': 1,
 'mean_occupation_satellites_assembias_param1': 0.2,
 u'scatter_model_param1': 0.2,
 u'smhm_beta_0': 0.43,
 u'smhm_beta_a': 0.18,
 u'smhm_delta_0': 0.56,
 u'smhm_delta_a': 0.18,
 u'smhm_gamma_0': 1.54,
 u'smhm_gamma_a': 2.52,
 u'smhm_m0_0': 10.72,
 u'smhm_m0_a': 0.59,
 u'smhm_m1_0': 12.35,
 u'smhm_m1_a': 0.3}
"""
data_table = np.loadtxt("../data/wp_all.dat")
rp_data = data_table[:,0]
wp_data = data_table[:,-2]
wp_error = data_table[:,-1]
#initial theta

theta0 = np.array([1. , 1.47 , -.13 , 0.859 , 10.62 , 1. , 0.2 , 0.2 , 0.43 , .18 , 0.56 ,0.18 , 1.54 , 2.52 , 10.72 , 0.59 , 12.35 , 0.3])
#theta0 = np.array([0.5 , 0.5])


def lnlike(theta):

    rp , wp = wp_model.model(theta , Mstar = 10.6)

    
    f_interp = interpolate.interp1d(rp , 
                                    wp,
                                    kind = "cubic")

    wp_interp = f_interp(rp_data)
    
    pylab.errorbar(rp_data , wp_data , yerr = wp_error)
    pylab.xscale("log")
    pylab.yscale("log")
    pylab.loglog(rp_data , wp_interp , "r-")
    pylab.show()
    
    return -0.5 * np.sum((wp_data - wp_interp) ** 2. / (wp_error ** 2.))


import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.fmin_l_bfgs_b(nll, theta0, fprime = None, args=(), approx_grad = True, disp = 1 ) 
#bounds = [(-1. , 1.) , (-1. , 1.)], disp = 1)

