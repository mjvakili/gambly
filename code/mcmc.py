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
import biased_hod
from biased_hod import MCMC_model
from prior import PriorRange
from halo_utils import load_project_halocat 


def lnPost(theta, box, **kwargs):
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

        data_nbar , data_wp = fake_obs[0] , fake_obs[1]
        nbar_var , wp_cov = fake_obs_icov[0] , fake_obs_icov[1]

    	kwargs.pop('data', None)
    	kwargs.pop('data_icov', None)
    	prior_range = kwargs['prior_range']
    	# Likelihood
    	model_obvs = generator(theta, prior_range)
        
        model_nbar , model_wp = model_obvs[0] , model_obvs[1]
       
        nbar_var = nbar_var
        
        res_nbar = model_nbar - data_nbar
        res_wp = model_wp - data_wp

        f_vol = Data.load_Volume_corrector(Mr, box)**-1.

        f_bias = (400. - len(res_wp) -2.)/(400. - 1.)

        neg_chisq_nbar = -0.5*(res_nbar**2.)/(nbar_var)
        neg_chisq_wp = -0.5 * f_bias * f_vol * np.sum(np.dot(res_wp , solve(wp_cov , res_wp)))
    	neg_chisq = neg_chisq_nbar + neg_chisq_wp

        print "neg_chi_tot" , neg_chisq
        return neg_chisq

    lp = lnprior(theta , **kwargs)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, **kwargs)


def mcmc_mpi(Nwalkers, Niters, Mr, box, it, start_chain, chain_file_name, prior_name = 'first_try'): 
    '''
    Parameters
    -----------
    - Nwalker : 
        Number of walkers
    - Nchains : 
        Number of MCMC chains   
    '''
    #data and covariance matrix
    fake_obs_icov = Data.load_covariance(Mr)
    fake_obs = Data.load_data(Mr)
        
    # True HOD parameters
    data_hod = Data.load_dechod_random_guess(Mr)
    Ndim = len(data_hod)

    # Priors
    prior_min, prior_max = PriorRange(prior_name, Mr, box)
    prior_range = np.zeros((len(prior_min),2))
    prior_range[:,0] = prior_min
    prior_range[:,1] = prior_max
    
    #Initializing Walkers
    
    if it == 0:

       random_guess = data_hod
       pos0 = np.repeat(random_guess, Nwalkers).reshape(Ndim, Nwalkers).T + \
                         5.e-3 * np.random.randn(Ndim * Nwalkers).reshape(Nwalkers, Ndim)
    else:
       
       pos0 = start_chain
       
   
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
            'Mr': Mr,
	    'box': box
            }
    sampler = emcee.EnsembleSampler(Nwalkers, Ndim, lnPost, pool=pool, kwargs=hod_kwargs)

    cnt = 0

    for result in sampler.sample(pos0, iterations = Niters, storechain=False):
        position = result[0]
        sample_file = h5py.File(chain_file_name)
        sample_file["mcmc"][cnt] = position
        sample_file.close()
        print cnt
        cnt += 1
        pass
    pool.close()

if __name__=="__main__": 
   
    continue_chain = True
    Nwalkers = int(sys.argv[1])
    print 'N walkers = ', Nwalkers
    Niters = int(sys.argv[2])
    print 'N iterations = ', Niters
    Mr = np.float(sys.argv[3])
    print 'Mr = ', np.float(Mr)
    box = sys.argv[4]
    print 'box = ', box
    it = np.int(sys.argv[5])
    print 'it = ', it
    Ndim = 7 # the number of hod parameters
  
    halocat = load_project_halocat(box)
    generator = MCMC_model(Mr, box, halocat)
    
    if it == 0: 

        chain_file_name = ''.join(['/disks/shear14/mj/mcmc/','mcmc_chain_Mr',str(Mr),'_box_',box,'.hdf5'])
    	sample_file = h5py.File(chain_file_name , 'w')
    	sample_file.create_dataset("mcmc",(Niters, Nwalkers, Ndim), data = np.zeros((Niters, Nwalkers , Ndim)))
    	sample_file.close()

    if it == 1:	
        
        old_chain_file_name = ''.join(['/disks/shear14/mj/mcmc/','mcmc_chain_Mr',str(Mr),'_box_',box,'.hdf5'])
	old_sample = h5py.File(old_chain_file_name, 'r')
        start_chain = old_sample["mcmc"][-1,:,:]
	old_sample.close()

        chain_file_name = ''.join(['/disks/shear14/mj/mcmc/','mcmc_chain_Mr',str(Mr),'_box_',box,'_1.hdf5'])
    	sample_file = h5py.File(chain_file_name , 'w')
    	sample_file.create_dataset("mcmc",(Niters, Nwalkers, Ndim), data = np.zeros((Niters, Nwalkers , Ndim)))
    	sample_file.close()

    else:

        old_chain_file_name = ''.join(['/disks/shear14/mj/mcmc/','mcmc_chain_Mr',str(Mr),'_box_',box,'_'+str(it-1)+'.hdf5'])
	old_sample = h5py.File(old_chain_file_name, 'r')
        start_chain = old_sample["mcmc"][-2,:,:]
	old_sample.close()
        print start_chain
        chain_file_name = ''.join(['/disks/shear14/mj/mcmc/','mcmc_chain_Mr',str(Mr),'_box_',box,'_'+str(it)+'.hdf5'])
    	sample_file = h5py.File(chain_file_name , 'w')
    	sample_file.create_dataset("mcmc",(Niters, Nwalkers, Ndim), data = np.zeros((Niters, Nwalkers , Ndim)))
    	sample_file.close()
    print start_chain
    mcmc_mpi(Nwalkers, Niters, Mr, box, it, start_chain, chain_file_name)
