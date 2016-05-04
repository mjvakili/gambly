import os 
import sys
import numpy as np
import emcee
from numpy.linalg import solve
from emcee.utils import MPIPool
from numpy.linalg import solve
import h5py
# --- Local ---
import util
import data as Data
from dechod import MCMC_model
from prior import PriorRange


def lnPost(theta, **kwargs):
    def lnprior(theta, **kwargs):
        '''log prior 
        '''
        fake_obs = kwargs['data']
    	fake_obs_icov = kwargs['data_icov']
    	kwargs.pop('data', None)
    	kwargs.pop('data_icov', None)
    	prior_range = kwargs['prior_range']
    	prior_min = prior_range[:,0]
    	prior_max = prior_range[:,1]
    	# Prior 
        if prior_min[0] < theta[0] < prior_max[0] and \
       	    prior_min[1] < theta[1] < prior_max[1] and \
            prior_min[2] < theta[2] < prior_max[2] and \
            prior_min[3] < theta[3] < prior_max[3] and \
            prior_min[4] < theta[4] < prior_max[4] and \
            prior_min[5] < theta[5] < prior_max[5]:
                return 0
    
        else:
		return -np.inf

    def lnlike(theta, **kwargs):

    	fake_obs = kwargs['data']
    	fake_obs_icov = kwargs['data_icov']
    	kwargs.pop('data', None)
    	kwargs.pop('data_icov', None)
    	prior_range = kwargs['prior_range']
    	# Likelihood
    	model_obvs = generator(theta, prior_range)
        #print "model=" , model_obvs
        res = fake_obs - model_obvs
        f = 1.
        neg_chisq = -0.5*f*np.sum(np.dot(res , solve(fake_obs_icov , res)))
    	print "neg_chi_tot" , neg_chisq
        return neg_chisq

    lp = lnprior(theta , **kwargs)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, **kwargs)


def mcmc_mpi(Nwalkers, Niters, data_dict={'Mr':21}, prior_name = 'first_try'): 
    '''
    Parameters
    -----------
    - Nwalker : 
        Number of walkers
    - Nchains : 
        Number of MCMC chains   
    '''
    #data and covariance matrix
    fake_obs_icov = Data.load_covariance(**data_dict)
    fake_obs = Data.load_data(**data_dict)
        
    # True HOD parameters
    data_hod = Data.load_hod_random_guess(Mr=21)
    Ndim = len(data_hod)

    # Priors
    prior_min, prior_max = PriorRange(prior_name)
    prior_range = np.zeros((len(prior_min),2))
    prior_range[:,0] = prior_min
    prior_range[:,1] = prior_max
    
    # mcmc chain output file 
    chain_file_name = ''.join([util.mcmc_dir(),'.mcmc_chain.hdf5'])
 

    if os.path.isfile(chain_file_name) and continue_chain:   
        print 'Continuing previous MCMC chain!'
        sample = h5py.File(chain_file_name , "r") 
        Nchains = Niters - len(sample) # Number of chains left to finish 
        if Nchains > 0: 
            pass
        else: 
            raise ValueError
        print Nchains, ' iterations left to finish'

        # Initializing Walkers from the end of the chain 
        pos0 = sample[-Nwalkers:]
    else:
        # new chain 
        sample_file = h5py.File(chain_file_name , 'w')
        sample_file.create_dataset("mcmc",(Niters, Nwalkers, Ndim), data = np.zeros((Niters, Nwalkers , Ndim)))
        sample_file.close()
         
        # Initializing Walkers
        random_guess = data_hod
        pos0 = np.repeat(random_guess, Nwalkers).reshape(Ndim, Nwalkers).T + \
                         5.e-2 * np.random.randn(Ndim * Nwalkers).reshape(Nwalkers, Ndim)
    print "initial position of the walkers = " , pos0
    # Initializing MPIPool
    pool = MPIPool(loadbalance=True)
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    # Initializing the emcee sampler
    hod_kwargs = {
            'prior_range': prior_range, 
            'data': fake_obs, 
            'data_icov': fake_obs_icov, 
            'Mr': data_dict['Mr']
            }
    sampler = emcee.EnsembleSampler(Nwalkers, Ndim, lnPost, pool=pool, kwargs=hod_kwargs)

    cnt = 0

    # Initializing Walkers 
    for result in sampler.sample(pos0, iterations = Niters, storechain=False):
        position = result[0]
        sample_file = h5py.File(chain_file_name)
        sample_file["k"][cnt] = position
        sample_file.close()
        print cnt
        cnt += 1
        pass
    pool.close()


if __name__=="__main__": 
    generator = MCMC_model(Mr = 21)
    continue_chain = False
    Nwalkers = int(sys.argv[1])
    print 'N walkers = ', Nwalkers
    Niters = int(sys.argv[2])
    print 'N iterations = ', Niters
    mcmc_mpi(Nwalkers, Niters)
