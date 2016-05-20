import numpy as np
import util
from Corrfunc import _countpairs
from Corrfunc.utils import read_catalog
import os.path as path
from halotools.mock_observables.catalog_analysis_helpers import return_xyz_formatted_array
from halotools.empirical_models import enforce_periodicity_of_box
from halotools.mock_observables import jackknife_covariance_matrix

def read_catalog():
    '''read Age-matching catalog
    '''

    filename =  util.dat_dir()+'Mr19_AM.dat'
    cat = np.loadtxt(filename)

    return cat

def impose_luminosity_threshold(Mr):
    '''divid up the galaxy catalog to different 
       luminosity thresholds, Mr = -19, -20, -21'''
    
    cat = read_catalog()
    luminosity = cat[:,9]
    reduced_cat = cat[luminosity>-1.*Mr]

    return reduced_cat

def save_luminosity_threshold(Mr):
    '''save a galaxy catalog with the imposed 
       luminosity threshold'''

    reduced_cat = impose_luminosity_threshold(Mr)
    filename = util.dat_dir+'AM_Mr'+str(Mr)+'dat'
    np.savetxt(filename , reduced_cat)

    return None

def measure_nbar_clustering(Mr):
    '''measure wp for the galaxy catalog
       with the luminosity threshold Mr'''

    #Corrfunc settings for wp measurements:

    boxsize = 250
    nthreads = 4
    pimax = 40.0
    binfile = path.join(path.dirname(path.abspath(__file__)),
                        "../", "bin")
    autocorr = 1

    filename = util.dat_dir+'AM_Mr'+str(Mr)+'dat'
    cat = np.loadtxt(filename)
    x, y, z = cat[:,1:4]
    vx, vy, vz = cat[:,4:7]
    #applying RSD:
    pos = return_xyz_formatted_array(x, y, z, velocity = vz, velocity_distortion_dimension = 'z')
    #enforcing periodic boundary conditions:
    pos = enforce_periodicity_of_box(pos, self.boxsize)
    pos = pos.astype(np.float32)
    x, y, z = pos[:,0] , pos[:,1] , pos[:,2]
    wp_result = _countpairs.countpairs_wp(boxsize, pimax, 
                                   nthreads, binfile, 
				   x, y, z)
    nbar = 1.*len(pos)/boxsize**3.
    wp = np.array(results_wp)[:,3]

    return nbar , wp    
    
def save_nbar_clustering(Mr):
    '''save nbar and wp measured for 
       the luminosity threshold Mr'''

    wp_filename = util.dat_dir()+'wp_AM_Mr'+str(Mr)+'.dat'
    nbar_filename = util.dat_dir()+'nbar_AM_Mr'+str(Mr)+'.dat'
    nbar , wp = measure_nbar_clustering(Mr)
    np.savetxt(nbar_filename ,  nbar)
    np.savetxt(wp_filename ,  wp)

    return None

def divid_box(nsub):
    '''divide the box for jackknife 
       stuff
       Nsub : the number of subboxes along
       each dimension of the box'''

    return None 
     
def edge(index , nsub):
    '''returns edges of a sub-box of 
       a given index
    '''
    box_size = 250.
    subbox_size = 1.*box_size / nsub

    zi = (index / (nsub**2)) * subbox_size
    i2 = index % (nsub**2)
    yi = (i2 / nsub) * subbox_size
    i3 = i2 % nsub
    xi = (i3) * subbox_size

    return xi , yi , zi
    
         
def mask_galaxy_table(pos , subvol_index , nsub):
    '''masks the positions of galaxies
    '''
    box_size = 250.
    subbox_size = 1.*box_size / nsub
    
    xi , yi , zi  = edge(subvol_index, nsub)
    submask = np.where((xi <pos[:, 0]) * \
                       (pos[:, 0] < xi + subbox_size) * \
                       (yi <pos[:, 1]) * \
                       (pos[:, 1] < yi + subbox_size) * \
                       (zi <pos[:, 2]) *  \
                       (pos[:, 2] < zi + subbox_size))
    
    return submask

def compute_jackknife_covariance(Mr , nsub):
    
    filename = util.dat_dir+'AM_Mr'+str(Mr)+'dat'
    cat = np.loadtxt(filename)
    pos = cat[:,1:4]  #positions of all galaxies in the box
    number_of_subboxes = nsub ** 3
   
    #Corrfunc settings for wp measurements:
    
    box_size = 250.
    subbox_size = box_size / nsub
    nthreads = 4
    pimax = 40.0
    binfile = path.join(path.dirname(path.abspath(__file__)),
                        "../", "bin")
    len_wp = len(np.loadtxt(binfile))
    wps = np.zeros((number_of_subboxes , len_wp)) 
    nbars = np.zeros((number_of_subboxes , 1))

    for subvol_index in xrange(number_of_subboxes):

        mask = mask_galaxy_table(pos , subvol_index , nsub)
        sub_cat = cat[mask]
        
        x, y, z = cat[:,1:4]
        vx, vy, vz = cat[:,4:7]
        #applying RSD:
        rsd_pos = return_xyz_formatted_array(x, y, z, velocity = vz, velocity_distortion_dimension = 'z')
        #enforcing periodic boundary conditions:
        rsd_pos = enforce_periodicity_of_box(rsd_pos, self.boxsize)
        rsd_pos = rsd_pos.astype(np.float32)
        x, y, z = rsd_pos[:,0] , rsd_pos[:,1] , rsd_pos[:,2]
        wp_result = _countpairs.countpairs_wp(boxsize, pimax, 
                                   nthreads, binfile, 
				   x, y, z)
        nbars[subvol_index] = 1.*len(rsd_pos)/subbox_size**3.
        wps[subvol_index] = np.array(results_wp)[:,3]
           
    nbar_covar = jackknife_covariance_matrix(nbars)    
    wp_covar = jackknife_covariance_matrix(wp)    

    return nbar_covar , wp_covar

def save_jackknife_covariance(Mr , nsub):

    wpcov_filename = util.dat_dir()+'wpcov_AM_Mr'+str(Mr)+'.dat'
    nbarcov_filename =  util.dat_dir()+'nbarcov_AM_Mr'+str(Mr)+'.dat'

    nbarcov , wpcov = compute_jackknife_covariance(Mr , nsub)

    np.savetxt(nbarcov_filename ,  nbarcov)
    np.savetxt(wpcov_filename ,  wpcov)
   
    return None 

if __name__ == "__main__":

     impose_luminosity_threshold(19)   
