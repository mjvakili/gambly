import numpy as np
import os

def fig_dir():
    '''
    figure directory
    '''
    fig_dir = os.path.dirname(os.path.realpath(__file__)).split('code')[0]+'fig/'
    return fig_dir

def dat_dir():
    '''
    Dat directory
    '''
    return os.path.dirname(os.path.realpath(__file__)).split('code')[0]+'dat/'
