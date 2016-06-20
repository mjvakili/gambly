import numpy as np
import h5py

def gelman_rubin(chain):
    ssq = np.var(chain, axis=1, ddof=1)
    W = np.mean(ssq, axis=0)
    xb = np.mean(chain, axis=1)
    xbb = np.mean(xb, axis=0)
    m = chain.shape[0]
    n = chain.shape[1]
    B = n / (m - 1.) * np.sum((xbb - xb)**2., axis=0)
    var_x = (n - 1.) / n * W + 1. / n * B
    Rhat = np.sqrt(var_x / W)
    return Rhat

if __name__=='__main__':
   directory = "/export/bbq2/mj/chains/"
   filename = directory+"group_mcmc_chain_Mr18.0.hdf5"
   chain = h5py.File(filename)["mcmc"][500:2100]
   print gelman_rubin(chain)
