
from halotools.mock_observables import wp
from halotools.empirical_models.factories.mock_helpers import three_dim_pos_bundle
import numpy as np
from halotools.empirical_models import PrebuiltHodModelFactory
#from halotools.mock_observables import return_xyz_formatted_array

rbins = np.logspace(-1, 1.25, 15)
rmax = rbins.max()
rbin_centers = (rbins[1:] + rbins[0:-1])/2.
pbins = np.linspace(.1 , 40.1, 2) #los radial bins
rbins = np.logspace(np.log10(.1) , np.log10(20) , 17) #projected radial bins
rbins_center =.5*(rbins[1:] + rbins[:-1]) 


pi_max = 40.
rp_bins = np.logspace(np.log10(.01),np.log10(20),18)
pi_bins = pbins #np.logspace(-1,np.log10(pi_max),17)
rp_bin_centers = (rp_bins[:-1] + rp_bins[1:])/2.
#print rp_bins
#print rp_bin_centers
def proj_clustering(pos , rbins, pbins , cellsize):

    return wp(sample1 = pos, rp_bins = rbins, 
              pi_bins = pbins, sample2 = None, 
              period = np.array([250,250,250]), approx_cell1_size = cellsize)

def model(theta, Mstar, prior_range = None):

    '''
    Given theta (HOD parameters) and stellar mass threshold Mstar, 
    compute wp
    not sure prior_range is necessary here
    '''

    model = PrebuiltHodModelFactory('hearin15', threshold = Mstar, 
                                     redshift = 0.0 ,
                                     central_assembias_strength = 0.,
                                     satellite_assembias_strength = 0.)

    model.param_dict['alphasat'] = theta[0]#1.0
    model.param_dict['bcut'] = theta[1]#1.47, 
    model.param_dict['betacut'] = theta[2]#0.13
    model.param_dict['betasat'] =  theta[3]#0.859
    model.param_dict['bsat'] = theta[4]#0.62
    model.param_dict['mean_occupation_centrals_assembias_param1'] = theta[5] # 0
    model.param_dict['mean_occupation_satellites_assembias_param1'] = theta[6] #.75
    model.param_dict['scatter_model_param1'] = theta[7]#0.2
    model.param_dict['smhm_beta_0'] = theta[8]# 0.43 
    model.param_dict['smhm_beta_a'] = theta[9]#0.18 
    model.param_dict['smhm_delta_0'] = theta[10]#0.56
    model.param_dict['smhm_delta_a'] = theta[11]#0.18, 
    model.param_dict['smhm_gamma_0']= theta[12]#1.54 
    model.param_dict['smhm_gamma_a'] = theta[13]#2.52, 
    model.param_dict['smhm_m0_0'] = theta[14]# 10.72, 
    model.param_dict['smhm_m0_a'] = theta[15]# 0.59 
    model.param_dict['smhm_m1_0'] = theta[16]# 12.35, 
    model.param_dict['smhm_m1_a'] = theta[17]#0.3
 
    
    model.populate_mock(simname = 'bolshoi', redshift = 0, halo_finder = 'rockstar')
    #x = model.mock.galaxy_table['x']
    #y = model.mock.galaxy_table['y']
    #z = model.mock.galaxy_table['z']
    #all_positions = return_xyz_formatted_array(x, y, z)    
    all_positions = three_dim_pos_bundle(model.mock.galaxy_table, 'x', 'y', 'z')
    
    wp_all = wp(all_positions, rp_bins, pi_bins,
                period=model.mock.Lbox, estimator = 'Landy-Szalay', num_threads = 1)
    #wp = proj_clustering(pos , rbins , pbins , cellsize = [rmax, rmax, rmax])
    print 1
    return rp_bin_centers, wp_all
