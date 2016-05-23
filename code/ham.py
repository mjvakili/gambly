'''
Module for handling the Abundance Matching catalog data and their
corresponding covariance matrices
'''

import numpy as np
import util

def load_data(Mr , style):
    '''loads wp and nbar
    '''
    return [load_nbar(Mr, style) , load_wp(Mr, style)]

def load_covariance(Mr, style):
    '''loads wp and nbar
    '''
    return [load_nbar_variance(Mr,style) , load_wp_covariance(Mr,style)]

def load_wp(Mr, style): 
    ''' loads wp 
    
    Parameters
    ----------
    Mr : (int)
        Absolute r-band magnitude threshold. Default is M_r = -21
    '''
    
    if style == "SHAM":
       style = "SHAM"
    elif style == "AM":
       style = "AM"

    data_file = ''.join([util.dat_dir(),
            "wp_"+str(style)+"_Mr"+str(Mr)+".dat"]) 
    wp = np.loadtxt(data_file)

    return wp

def load_wp_covariance(Mr , style): 
    ''' loads the jackknife covariance matrix associated with  wp
    Parameters
    ----------
    Mr : (int)
        Absolute r-band magnitude threshold. Default is M_r = -21
    '''
    if style == "SHAM":
       style = "SHAM"
    elif style == "AM":
       style = "AM"

    data_file = ''.join([util.dat_dir(),
            "wpcov_"+str(style)+"_Mr"+str(Mr)+".dat"]) 
    wpcov = np.loadtxt(data_file)

    return wpcov

def load_nbar(Mr, style):
    '''load the number density of the data
    '''

    if style == "SHAM":
       style = "SHAM"
    elif style == "AM":
       style = "AM"

    data_file = ''.join([util.dat_dir(),
            "nbar_"+str(style)+"_Mr"+str(Mr)+".dat"]) 
    nbar = np.loadtxt(data_file)

    return nbar

def load_nbar_variance(Mr, style):
    '''load the variance of the number density of the data
    '''
    if style == "SHAM":
       style = "SHAM"
    elif style == "AM":
       style = "AM"

    data_file = ''.join([util.dat_dir(),
            "nbarcov_"+str(style)+"_Mr"+str(Mr)+".dat"]) 
    nbarcov = np.loadtxt(data_file)

    return nbarcov


def load_hod_random_guess(Mr):
    '''random guess for the walkers'''

    if Mr == 21.:
        pos = [12.59 , 0.49 , 12.78 , 1.14 , 13.99]
    if Mr == 20.5:
        pos = [11.84 , 0.39 , 12.79 , 1.12 , 13.58]
    if Mr == 19.5:
        pos = [12.59 , 0.49 , 12.78 , 1.14 , 13.99]
    if Mr == 20.:
        pos = [11.38 , 0.26 , 12.02 , 1.06 , 13.31]
    if Mr == 19.5:
        pos = [11.69 , 0.28 , 11.75 , 1.05 , 13.01]
    if Mr == 19.:
        pos = [11.49 , 0.26 , 11.6 , 1.02 , 12.83]

    return pos

def load_dechod_random_guess(Mr):
    '''random guess for the walkers'''

    if Mr == 21.:
        pos = [12.59 , 0.49 , 12.78 , 1.14 , 13.99 , 0.5 , 0.5]
    if Mr == 20.5:
        pos = [11.84 , 0.39 , 12.79 , 1.12 , 13.58 , 0.5 , 0.5]
    if Mr == 19.5:
        pos = [12.59 , 0.49 , 12.78 , 1.14 , 13.99 , 0.5 , 0.5]
    if Mr == 20.:
        pos = [11.38 , 0.26 , 12.02 , 1.06 , 13.31 , 0.5 , 0.5]
    if Mr == 19.5:
        pos = [11.69 , 0.28 , 11.75 , 1.05 , 13.01 , 0.5 , 0.5]
    if Mr == 19.:
        pos = [11.49 , 0.26 , 11.6 , 1.02 , 12.83 , 0.5 , 0.5]

    return pos

def load_Volume_corrector(Mr):
    '''calculate the volume correction term of the 
       covariance matrix'''

    return 2.

