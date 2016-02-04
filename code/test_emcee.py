import sys
import numpy as np
import emcee
from emcee.utils import MPIPool

# Choose the "true" parameters.
m_true = -0.9594
b_true = 4.294
f_true = 0.534

# Generate some synthetic data from the model.
N = 50
x = np.sort(10*np.random.rand(N))
yerr = 0.1+0.5*np.random.rand(N)
y = m_true*x+b_true
y += np.abs(f_true*y) * np.random.randn(N)
y += yerr * np.random.randn(N)

A = np.vstack((np.ones_like(x), x)).T
C = np.diag(yerr * yerr)
cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))
import time
def lnlike(theta, x, y, yerr):
    a = time.time()
    m, b, lnf = theta
    model = m * x + b
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    print time.time() - a
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))
m_ml, b_ml, lnf_ml = result["x"]

def lnprior(theta):
    m, b, lnf = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
        return 0.0
    return -np.inf

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)


ndim, nwalkers = 3, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

#pool = MPIPool()
#if not pool.is_master():
#    pool.wait()
#    sys.exit(0)

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
f = open("chain.dat", "w")
f.close()
for result in sampler.sample(pos, iterations=600, storechain=False):
        #print result
        position = result[0]
        print position.shape
        f = open("chain.dat", "a")
        for k in range(position.shape[0]):
	    output_str = '\t'.join(position[k].astype('str')) + '\n'
            f.write(output_str)
        f.close()

sampler.run_mcmc(pos, 600)

samples = sampler.chain[:, 100:, :].reshape((-1, ndim))

#pool.close()

#import corner
#fig = corner.corner(samples, labels=["$m$", "$b$", "$\ln\,f$"],
#                      truths=[m_true, b_true, np.log(f_true)])
#fig.savefig("/home/mj/public_html/triangle_final.png")
#
#samples2 = np.loadtxt("chain.dat")
#
#fig = corner.corner(samples2[100:], labels=["$m$", "$b$", "$\ln\,f$"],
#                      truths=[m_true, b_true, np.log(f_true)])
#
#fig.savefig("/home/mj/public_html/triangle_incremental.png")



