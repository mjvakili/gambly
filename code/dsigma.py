import numpy as np
import pylab
from matplotlib import lines as mlines
from matplotlib import pyplot as plt
import numpy as np


eps = 1.e-13

def get_sigma(xigm , r):

    """
    computes projected matter surface density as a function of 
    comoving radius (Mpc/h)
    from a tabulated \xi_{gm} and bins of radii
    """

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

     pylab.plot( r , dsigma(xigm , rbins))
     pylab.xlabel(r"$R$")
     pylab.ylabel(r"$\Sigma(R)$")
     pylab.show()
