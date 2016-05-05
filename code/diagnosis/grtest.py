'''
implementation of Gelman-Rubin convergence test
'''
import numpy as np

def single_parameter_gr_test(chains):
    """
    inputs:
      chains : MCMC samples for one parameter.
      shape = (nwalkers , niters)
    returns:
      potential scale reduction factor
      and the variance of the distribution      
    """
    nwalkers =  chains.shape[0] 
    niters =  chains.shape[1]
    

    #Discarding the first half of draws:
    chains = chains[: , niters/2:]
    nwalkers , niters = chains.shape[0] , chains.shape[1]

    #Calculating the within-chain variance:
    W = np.mean(np.var(chains, axis=1))

    #Calculating the between-chain variance:
    chains_means = np.mean(chains, axis=1)
    mean_of_chains_means = np.mean(chains_means)
    B = (niters/(nwalkers-1.0)) * np.sum((chains_means - mean_of_chains_means)**2.)

    # Estimating the variance of distribution:
    V = (1. - 1./niters) * W + (1./niters) * B

    # Calculating the potential scale reduction factor:
    R = np.sqrt(V/W)
    return R , V


def gr_test(sample , nwalkers , nburnins , npars):

    """
    inputs: 
      sample = an emcee sample
      nwalkers = number of walkers
    returns:
      Rs = npar-dimensional vector of the 
      potential scale reduction factors
      Vs = npar-dimensional vector of the 
      variances
    """


    #npar = len(sample[0])
    niters = len(sample)
    chain_ensemble = sample.reshape(niters , nwalkers , npars)
    chain_ensemble = chain_ensemble[nburnins: , :]
    Rs = np.zeros((npars))
    Vs = np.zeros((npars))
    for i in range(npars):
     
        chains = chain_ensemble[ : , : , i].T
        Rs[i] = single_parameter_gr_test(chains)[0]
        Vs[i] = single_parameter_gr_test(chains)[1]

    return Rs , Vs
      

if __name__ == "__main__":

  ###########NOTE: This part is provided by the user#########
  import h5py
  sample = h5py.File("emcee2.hdf5")["k"]
  nwalkers  = 6
  nburnins = 1700
  npars = 3
  print gr_test(sample , nwalkers , nburnins , npars)
