'''
Module for handling the data anc coavariance matrices
'''

import numpy as np
import util

def load_data(Mr):
    '''loads wp and nbar
    '''
    return [load_nbar(Mr) , load_wp(Mr)]

def load_covariance(Mr):
    '''loads wp and nbar
    '''
    return [load_nbar_variance(Mr) , load_wp_covariance(Mr)]

def load_wp(Mr): 
    ''' loads wp 
    
    Parameters
    ----------
    Mr : (int)
        Absolute r-band magnitude threshold. Default is M_r = -21
    '''
    if Mr == 21.5:
        data_file = ''.join([util.dat_dir(),
            'wpxi_dr72_bright0_mr21.5_z0.198_nj400'])
    if Mr == 21.:
        data_file = ''.join([util.dat_dir(),
            'wpxi_dr72_bright0_mr21.0_z0.159_nj400']) 
    if Mr == 20.5:
        data_file = ''.join([util.dat_dir(),
            'wpxi_dr72_bright0_mr20.5_z0.132_nj400']) 
    if Mr == 20.:
        data_file = ''.join([util.dat_dir(),
            'wpxi_dr72_bright0_mr20.0_z0.106_nj400']) 
    if Mr == 19.5:
        data_file = ''.join([util.dat_dir(),
            'wpxi_dr72_bright0_mr19.5_z0.085_nj400']) 
    if Mr == 19.:
        data_file = ''.join([util.dat_dir(),
            'wpxi_dr72_bright0_mr19.0_z0.064_nj400'])
    if Mr == 18.5:
        data_file = ''.join([util.dat_dir(),
            'wpxi_dr72_bright0_mr18.5_z0.053_nj400'])
    if Mr == 18.:
        data_file = ''.join([util.dat_dir(),
            'wpxi_dr72_bright0_mr18.0_z0.041_nj400'])
    wp = np.loadtxt(data_file)[:,1]

    return wp

def load_nbar(Mr):
    '''load the number density of the data
    '''
    if Mr == 21.5:
        nbar = 0.29 * 10.**-3.
    if Mr == 21.:
        nbar = 1.16 * 10**-3.
    if Mr == 20.5:
        nbar = 3.13 * 10**-3.
    if Mr == 20.:
        nbar = 6.37 * 10**-3.
    if Mr == 19.5:
        nbar = 11.64 * 10**-3.
    if Mr == 19.:
        nbar = 15.66 * 10**-3.
    if Mr == 18.5:
        nbar = 22.25 * 10**-3.
    if Mr == 18.:
        nbar = 31.42 * 10**-3.
    return nbar

def load_wp_covariance(Mr): 
    ''' loads the jackknife covariance matrix associated with  wp
    
    Parameters
    ----------
    Mr : (int)
        Absolute r-band magnitude threshold. Default is M_r = -21
    '''
    if Mr == 21.5:
        data_file = ''.join([util.dat_dir(),
            'wpxicov_dr72_bright0_mr21.5_z0.198_nj400'])
    if Mr == 21.:
        data_file = ''.join([util.dat_dir(),
            'wpxicov_dr72_bright0_mr21.0_z0.159_nj400']) 
    if Mr == 20.5:
        data_file = ''.join([util.dat_dir(),
            'wpxicov_dr72_bright0_mr20.5_z0.132_nj400']) 
    if Mr == 20.:
        data_file = ''.join([util.dat_dir(),
            'wpxicov_dr72_bright0_mr20.0_z0.106_nj400']) 
    if Mr == 19.5:
        data_file = ''.join([util.dat_dir(),
            'wpxicov_dr72_bright0_mr19.5_z0.085_nj400']) 
    if Mr == 19.:
        data_file = ''.join([util.dat_dir(),
            'wpxicov_dr72_bright0_mr19.0_z0.064_nj400'])
    if Mr == 18.5:
        data_file = ''.join([util.dat_dir(),
            'wpxicov_dr72_bright0_mr18.5_z0.053_nj400'])
    if Mr == 18.:
        data_file = ''.join([util.dat_dir(),
            'wpxicov_dr72_bright0_mr18.0_z0.041_nj400'])
    wpcov = np.loadtxt(data_file)[:12 , :12]

    return wpcov

def load_nbar_variance(Mr):
    '''load the variance of the number density of the data
    '''
    if Mr == 21.5:
        nbarerr= 0.03 * 10**-3.
    if Mr == 21.:
        nbarerr= 0.12 * 10**-3.
    if Mr == 20.5:
        nbarerr = 0.3 * 10**-3.
    if Mr == 20.:
        nbarerr = 0.75 * 10**-3.
    if Mr == 19.5:
        nbarerr = 1.27 * 10**-3.
    if Mr == 19.:
        nbarerr = 2.06 * 10**-3.
    if Mr == 18.5:
        nbarerr = 2.70 * 10**-3.
    if Mr == 18.:
        nbarerr = 4.03 * 10**-3.
    return nbarerr ** 2.

def load_hod_random_guess(Mr):
    '''random guess for the walkers'''

    if Mr == 21.5:
        pos = [13.13 , 0.72 , 13.53 , 1.14 , 14.52]
    if Mr == 21.:
        pos = [12.78 , 0.49 , 12.59 , 1.14 , 13.99]
    if Mr == 20.5:
        pos = [12.23 , 0.18 , 12.42 , 1.06 , 13.58]
    if Mr == 19.5:
        pos = [11.95 , 0.26 , 12.10 , 1.08 , 13.31]
    if Mr == 20.:
        pos = [11.55 , 0.26 , 11.6 , 1.12 , 13.01]
    if Mr == 19.5:
        pos = [11.69 , 0.28 , 11.75 , 1.05 , 13.01]
    if Mr == 19.:
        pos = [11.49 , 0.26 , 11.6 , 1.02 , 12.83]
    if Mr == 18.5:
        pos = [11.38 , 0.26 , 11.73 , 1.02 , 12.71]
    if Mr == 18.:
        pos = [11.18 , 0.1 , 11.6 , 0.97 , 12.48]

    return pos

def load_dechod_random_guess(Mr):
    '''random guess for the walkers'''
    
    if Mr == 21.5:
        pos = [13.13 , 0.72 , 13.53 , 1.14 , 14.52 , 0.0 , 0.0]
    if Mr == 21.:
        pos = [12.59 , 0.49 , 12.78 , 1.14 , 13.99 , 0.5 , 0.5]
    if Mr == 20.5:
        pos = [12.42 , 0.18 , 12.23 , 1.06 , 13.58 , 0.5 , 0.5]
    if Mr == 20.:
        pos = [12.10 , 0.26 , 11.95 , 1.08 , 13.31 , 0.5 , 0.5]
    if Mr == 19.5:
        pos = [11.80 , 0.28 , 11.67 , 1.05 , 13.01 , 0.5 , 0.5]
    if Mr == 19.:
        pos = [11.6 , 0.26 , 11.58 , 1.12 , 13.04 , 0.5 , 0.5]
    if Mr == 18.5:
        pos = [11.73 , 0.26 , 11.38 , 1.02 , 12.71 , 0.5 , 0.5]
    if Mr == 18.:
        pos = [11.57 , 0.1 , 11.18 , 0.97 , 12.48 , 0.5 , 0.5]
    return pos

def load_Volume_corrector(Mr, box):
    '''calculate the volume correction term of the 
       covariance matrix'''

    if box == "smd":   
       Vsim = 400. ** 3.
    elif box == "bolshoi_planck":
       Vsim = 250 ** 3.
    
    if Mr == 21.5:
        Vd = 134.65 * 10**6.
    if Mr == 21.:
        Vd = 71.74 * 10**6.
    if Mr == 20.5:
        Vd = 41.82 * 10**6.
    if Mr == 20.:
        Vd = 22.0 * 10**6.
    if Mr == 19.5:
        Vd = 11.44 * 10**6.
    if Mr == 19.:
        Vd = 4.87 * 10**6.

    return (1. + Vd / Vsim)

