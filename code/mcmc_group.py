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
import data_group as Data
from biased_hod_group import MCMC_model
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
            prior_min[5] < theta[5] < prior_max[5] and \
            prior_min[6] < theta[6] < prior_max[6]:
                return 0
    
        else:
		return -np.inf

    def lnlike(theta, **kwargs):

    	fake_obs = kwargs['data']
    	fake_obs_icov = kwargs['data_icov']

        data_nbar , data_gmf = fake_obs[0] , fake_obs[1]
        nbar_var , gmf_cov = fake_obs_icov[0] , fake_obs_icov[1]

    	kwargs.pop('data', None)
    	kwargs.pop('data_icov', None)
    	prior_range = kwargs['prior_range']
    	# Likelihood
    	model_obvs = generator(theta, prior_range)
        model_nbar , model_gmf = model_obvs[0] , model_obvs[1]
        
        res_nbar = model_nbar - data_nbar
        res_gmf = model_gmf - data_gmf

        f_vol = 1. #Data.load_Volume_corrector(Mr)**-1.

        f_bias = 1. #(400. - len(res_wp) -2.)/(400. - 1.)

        neg_chisq_nbar = -0.5*(res_nbar**2.)/(nbar_var)
        neg_chisq_gmf = -0.5*f_bias * f_vol * np.sum((res_gmf**2./gmf_cov))
    	neg_chisq = neg_chisq_nbar + neg_chisq_gmf

        print "neg_chi_tot" , neg_chisq
        return neg_chisq

    lp = lnprior(theta , **kwargs)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, **kwargs)


def mcmc_mpi(Nwalkers, Niters, Mr, prior_name = 'first_try', pois = False): 
    '''
    Parameters
    -----------
    - Nwalker : 
        Number of walkers
    - Nchains : 
        Number of MCMC chains   
    '''
    #data and covariance matrix
    fake_obs_icov = Data.load_covariance(Mr , pois = False)
    fake_obs = Data.load_data(Mr)
        
    # True HOD parameters
    data_hod = Data.load_dechod_random_guess(Mr)
    Ndim = len(data_hod)

    # Priors
    prior_min, prior_max = PriorRange(prior_name , Mr)
    prior_range = np.zeros((len(prior_min),2))
    prior_range[:,0] = prior_min
    prior_range[:,1] = prior_max
    
    # mcmc chain output file 
    chain_file_name = ''.join([util.mcmc_dir(),'group_nopoisson_mcmc_chain_Mr',str(Mr),'.hdf5'])
 

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
        print "chain_file_name=" , chain_file_name
 
        sample_file = h5py.File(chain_file_name , 'w')
        sample_file.create_dataset("mcmc",(Niters, Nwalkers, Ndim), data = np.zeros((Niters, Nwalkers , Ndim)))
        sample_file.close()
         
        # Initializing Walkers
        random_guess = data_hod
        pos0 = np.repeat(random_guess, Nwalkers).reshape(Ndim, Nwalkers).T + \
                         5.e-2 * np.random.randn(Ndim * Nwalkers).reshape(Nwalkers, Ndim)
    print "initial position of the walkers = " , pos0.shape
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
            'Mr': Mr
            }
    sampler = emcee.EnsembleSampler(Nwalkers, Ndim, lnPost, pool=pool, kwargs=hod_kwargs)

    cnt = 0

    # Initializing Walkers 
    for result in sampler.sample(pos0, iterations = Niters, storechain=False):
        position = result[0]
        sample_file = h5py.File(chain_file_name)
        sample_file["mcmc"][cnt] = position
        sample_file.close()
        print "iteration=" , cnt
        cnt += 1
        pass
    pool.close()

if __name__=="__main__": 
    continue_chain = False
    Nwalkers = int(sys.argv[1])
    print 'N walkers = ', Nwalkers
    Niters = int(sys.argv[2])
    print 'N iterations = ', Niters
    Mr = np.float(sys.argv[3])
    print 'Mr = ', np.float(Mr)
    generator = MCMC_model(Mr)
    mcmc_mpi(Nwalkers, Niters, Mr)
