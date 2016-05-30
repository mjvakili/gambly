'''
Module for handling the data and coavariance matrices from Berlind 2006
'''
import numpy as np
import util

def load_data(Mr):
    '''loads wp and nbar
    '''
    return [load_nbar(Mr) , load_gmf(Mr)]

def load_covariance(Mr , pois = False):
    '''loads wp and nbar
    '''
    return [load_nbar_variance(Mr) , load_gmf_covariance(Mr , pois = False)]

def load_gmf(Mr): 
    ''' loads wp 
    
    Parameters
    ----------
    Mr : (int)
        Absolute r-band magnitude threshold. Default is M_r = -21
    '''
    if Mr == 20.:
        data_file = ''.join([util.dat_dir(),
            'gmf_mr20.0.dat']) 
    if Mr == 19.:
        data_file = ''.join([util.dat_dir(),
            'gmf_mr19.0.dat'])
    if Mr == 18.:
        data_file = ''.join([util.dat_dir(),
            'gmf_mr18.0.dat'])

    gmf = np.loadtxt(data_file)[:,2]

    return gmf

def load_nbar(Mr):
    '''load the number density of the data
    '''
    if Mr == 20.:
        nbar = 0.00673
    if Mr == 19.:
        nbar = 0.01396
    if Mr == 18.:
        nbar = 0.02434
    return nbar

def load_gmf_covariance(Mr , pois = False): 
    ''' loads the jackknife covariance matrix associated with gmf
    
    Parameters
    ----------
    Mr : (int)
        Absolute r-band magnitude threshold. Default is M_r = -21
    pois : True r false, whether we have poisson error in likelihood or not
    '''
    if Mr == 20.:
        data_file = ''.join([util.dat_dir(),
            'gmf_mr20.0.dat']) 
    if Mr == 19.:
        data_file = ''.join([util.dat_dir(),
            'gmf_mr19.0.dat'])
    if Mr == 18.:
        data_file = ''.join([util.dat_dir(),
            'gmf_mr18.0.dat'])
    if pois == True:
        gmf_err = (np.loadtxt(data_file)[:,3]**2. + np.loadtxt(data_file)[:,4]**2.)**.5 
    if pois == False:
        gmf_err = np.loadtxt(data_file)[:,3] 
        
    return gmf_err ** 2.

def load_nbar_variance(Mr):
    '''load the variance of the number density of the data
    '''
    if Mr == 20.:
        nbarerr = 0.75 * 10**-3.
    if Mr == 19.:
        nbarerr = 2.06 * 10**-3.
    if Mr == 18.:
        nbarerr = 4.03 * 10**-3.
    return nbarerr ** 2.

def load_hod_random_guess(Mr):
    '''random guess for the walkers'''

    if Mr == 20.:
        pos = [12.10 , 0.26 , 11.95 , 1.08 , 13.31]
    if Mr == 19.:
        pos = [11.6 , 0.26 , 11.58 , 1.12 , 13.04]
    if Mr == 18.:
        pos = [11.57 , 0.1 , 11.18 , 0.97 , 12.48]

    return pos

def load_dechod_random_guess(Mr):
    '''random guess for the walkers'''

    if Mr == 20.:
        pos = [12.10 , 0.26 , 11.95 , 1.08 , 13.31 , 0.5 , 0.5]
    if Mr == 19.:
        pos = [11.6 , 0.26 , 11.58 , 1.12 , 13.04 , 0.5 , 0.5]
    if Mr == 18.:
        pos = [11.57 , 0.1 , 11.18 , 0.97 , 12.48 , 0.5 , 0.5]
    return pos

def load_Volume_corrector(Mr):
    '''calculate the volume correction term of the 
       covariance matrix'''
    Vsim = 250. ** 3.
    
    if Mr == 20.:
        Vd = 22.0 * 10**6.
    if Mr == 19.:
        Vd = 4.87 * 10**6.

    return (1. + Vd / Vsim)

