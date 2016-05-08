'''
Module to deal with prior range
'''
import numpy as np

def PriorRange(prior_name): 
    ''' Given prior dictionary name, return the prior range. 
    '''
    if prior_name is None: 
        prior_name = 'first_try'
    
    dict_table = prior_dict_table() 

    prior_dict = dict_table[prior_name]

    prior_min = prior_dict['prior_min']
    prior_max = prior_dict['prior_max']
    return [prior_min, prior_max]

def prior_dict_table(Mr): 
    ''' dictionary table of priors 
    '''
    if Mr ==21 :
      dict_table = { 
            'first_try': {
                'prior_min': [10.5, 0.1, 11.8, 0.85, 13., -1. , -1.],
                'prior_max': [14.5, 0.8, 13.8, 1.45, 15., 1. , 1.]
                }
              }
            
    if Mr ==20.5 :
      dict_table = { 
            'first_try': {
                'prior_min': [10.5, 0.1, 11.8, 0.85, 13., -1. , -1.],
                'prior_max': [14.5, 0.8, 13.8, 1.45, 15., 1. , 1.]
                }
              }
    if Mr ==20. :
      dict_table = { 
            'first_try': {
                'prior_min': [10.5, 0.1, 11.8, 0.85, 13., -1. , -1.],
                'prior_max': [14.5, 0.8, 13.8, 1.45, 15., 1. , 1.]
                }
              }
    if Mr ==19.5 :
      dict_table = { 
            'first_try': {
                'prior_min': [10.5, 0.1, 11.8, 0.85, 13., -1. , -1.],
                'prior_max': [14.5, 0.8, 13.8, 1.45, 15., 1. , 1.]
                }
              }
    if Mr ==19. :
      dict_table = { 
            'first_try': {
                'prior_min': [10.5, 0.1, 11.8, 0.85, 13., -1. , -1.],
                'prior_max': [14.5, 0.8, 13.8, 1.45, 15., 1. , 1.]
                }
              }

    return dict_table 
