import numpy as np
#import matplotlib.pyplot as plt
#plt.switch_backend("Agg")
import pylab
from matplotlib import lines as mlines
from matplotlib import pyplot as plt
import numpy as np

"""
from halotools.empirical_models import PrebuiltHodModelFactory
model = PrebuiltHodModelFactory('hearin15', threshold = 10.6, redshift = 0., 
                        central_assembias_strength = 0, 
                        satellite_assembias_strength = 0)
#model.param_dict['scatter_model_param1'] = 0.4 
baseline_model = PrebuiltHodModelFactory('leauthaud11', threshold = 10.6, redshift = 0.)
#baseline_model.param_dict['scatter_model_param1'] = 0.4
baseline_model.populate_mock()
rbins = np.loadtxt("../data/gm_all.dat")[:,0]/1000
print rbins.shape
gg = np.loadtxt("../data/gm_all.dat")[:,1]
print rbins
baseline_model.populate_mock()
a = np.logspace(np.log10(.02) , np.log10(80) , 25)
print .5*(a[1:]+a[:-1])
rr , gg_model = baseline_model.mock.compute_galaxy_matter_cross_clustering(rbins = a)

#rr , gg_model = baseline_model.mock.compute_galaxy_matter_cross_clustering(rbin = a)
plt.loglog(rr , gg_model , "b.")
plt.loglog(rbins , gg , "r.")
"""


eps = 1e-13
def sigma(xigm , rbins):

    sig = np.zeros_like((xigm))

    r = .5*(rbins[1:]+rbins[:-1])
    dr = rbins[1:]-rbins[:-1]

    f = (1 + xigm)*r*dr

    for j in range(len(r)):

        sig[j] =  np.sum(f[j:] * np.sqrt(r[j:]**2. - rbins[j]**2. + eps)**-1.)

    return sig 


def sigmabar(xigm , rbins):

    sigbar = np.zeros_like((xigm))
    
    sig = sigma(xigm , rbins)

    r = .5*(rbins[1:]+rbins[:-1])
    dr = rbins[1:]-rbins[:-1]

    for j in range(len(r)):

        sigbar[j] = 2.* (r[j]**-2.) * np.sum((sig*r*dr)[:j])

    return sigbar


def dsigma(xigm , rbins):

    A = np.pi * (rbins[-1:]**2. - rbins[1:]**2.) 
    R = rbins
    delta = np.zeros_like((xigm))
    r = .5*(rbins[1:]+rbins[:-1])
    sig = sigma(xigm , rbins)
    for j in range(len(r)):

        (np.pi * rbins[j]**2.)**-1. * np.sum((A*sig)[:j] - sig[j])

    return sig


def dsigmap(xigm , rbins):


    return sigma(xigm , rbins) - sigmabar(xigm , rbins)   


if __name__ =='__main__':

     import pylab

     rbins = np.logspace(np.log10(.1) , np.log10(20) , 50)
     r = .5*(rbins[1:]+rbins[:-1])
     xigm =  r**-3.

     print sigma(xigm , rbins)
     
     pylab.plot( r , dsigma(xigm , rbins))
     pylab.xlabel(r"$R$")
     pylab.ylabel(r"$\Sigma(R)$")
     pylab.show()

    
