'''
Module for handling the SHAM catalogs 
made by ChangHoon Hahn
'''
import h5py
import numpy as np
import util
from Corrfunc import _countpairs
from Corrfunc.utils import read_catalog
import os.path as path
from halotools.mock_observables.catalog_analysis_helpers import return_xyz_formatted_array
from halotools.empirical_models import enforce_periodicity_of_box


def read_catalog():
    '''read SHAM catalog
    '''

    filename =  util.dat_dir()+'bolshoi_a1.00231.mag_r.source_blanton.scatter0.2.Vpeak.hdf5'
    cat = h5py.File(filename, "r")['data']
    length_file = len(cat['x'][:])
    cat_file = np.zeros((length_file , 7))
    cat_file[:,0] = cat['x'][:]
    cat_file[:,1] = cat['y'][:]
    cat_file[:,2] = cat['z'][:]
    cat_file[:,3] = cat['vx'][:]
    cat_file[:,4] = cat['vy'][:]
    cat_file[:,5] = cat['vz'][:]
    cat_file[:,6] = cat['mag_r'][:]
    
    return cat_file

def impose_luminosity_threshold(Mr):
    '''divid up the galaxy catalog to different 
       luminosity thresholds, Mr = -19, -20, -21'''
    
    cat = read_catalog()
    luminosity = cat[:,-1]
    reduced_cat = cat[luminosity > Mr]

    return reduced_cat

def save_luminosity_threshold(Mr):
    '''save a galaxy catalog with the imposed 
       luminosity threshold'''

    reduced_cat = impose_luminosity_threshold(Mr)
    filename = util.dat_dir()+'SHAM_Mr'+str(Mr)+'dat'
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

    filename = util.dat_dir()+'SHAM_Mr'+str(Mr)+'dat'
    cat = np.loadtxt(filename)
    pos = cat[:,0:3]
    x, y, z = pos[:,0], pos[:,1], pos[:,2] 
    vel = cat[:,3:6]
    vx, vy, vz =  vel[:,0], vel[:,1], vel[:,2] 
    #applying RSD:
    pos = return_xyz_formatted_array(x, y, z, velocity = vz, velocity_distortion_dimension = 'z')
    #enforcing periodic boundary conditions:
    pos = enforce_periodicity_of_box(pos, boxsize)
    pos = pos.astype(np.float32)
    x, y, z = pos[:,0] , pos[:,1] , pos[:,2]
    wp_result = _countpairs.countpairs_wp(boxsize, pimax, 
                                   nthreads, binfile, 
				   x, y, z)
    nbar = 1.*len(pos)/boxsize**3.
    wp = np.array(wp_result)[:,3]

    return nbar , wp    
    
def save_nbar_clustering(Mr):
    '''save nbar and wp measured for 
       the luminosity threshold Mr'''

    wp_filename = util.dat_dir()+'wp_SHAM_Mr'+str(Mr)+'.dat'
    nbar_filename = util.dat_dir()+'nbar_SHAM_Mr'+str(Mr)+'.dat'
    nbar , wp = measure_nbar_clustering(Mr)
    np.savetxt(nbar_filename ,  np.array([nbar]))
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
    
         
def mask_catalog(cat , subvol_index , nsub):
    '''masks the SHAM catalogs
    '''
    box_size = 250.
    subbox_size = 1.*box_size / nsub
    
    xi , yi , zi  = edge(subvol_index, nsub)
    submask = np.where((xi <cat[:, 1]) * \
                       (cat[:, 1] < xi + subbox_size) * \
                       (yi <cat[:, 2]) * \
                       (cat[:, 2] < yi + subbox_size) * \
                       (zi <cat[:, 3]) *  \
                       (cat[:, 3] < zi + subbox_size))
    
    return cat[submask]

def mask_positions(pos , subvol_index , nsub):
    '''masks the positions of galaxies in
       SHAM catalogs'''

    box_size = 250.
    subbox_size = 1.*box_size / nsub
    
    xi , yi , zi  = edge(subvol_index, nsub)
    submask = np.where((xi <pos[:, 0]) * \
                       (pos[:, 0] < xi + subbox_size) * \
                       (yi <pos[:, 1]) * \
                       (pos[:, 1] < yi + subbox_size) * \
                       (zi <pos[:, 2]) *  \
                       (pos[:, 2] < zi + subbox_size))
    
    return pos[submask] - edge(subvol_index,nsub)

def compute_jackknife_covariance(Mr , nsub):
    
    box_size = 250.
    filename = util.dat_dir()+'SHAM_Mr'+str(Mr)+'dat'
    cat = np.loadtxt(filename)
    pos = cat[:,0:3]  #positions of all galaxies in the box
    vel = cat[:,3:6]  #velocities of all galaxies
    
    x , y , z = pos[:,0], pos[:,1], pos[:,2]
    vx, vy, vz = vel[:,0], vel[:,1], vel[:,2]
    rsd_pos = return_xyz_formatted_array(x, y, z, velocity = vz, velocity_distortion_dimension = 'z')
    rsd_pos = enforce_periodicity_of_box(rsd_pos, box_size)
    
    number_of_subboxes = nsub ** 3
   
    #Corrfunc settings for wp measurements:
    
    subbox_size = box_size / nsub
    nthreads = 4
    pimax = 25.0
    binfile = path.join(path.dirname(path.abspath(__file__)),
                        "../", "bin")
    len_wp = len(np.loadtxt(binfile))
    wps = np.zeros((number_of_subboxes , len_wp)) 
    nbars = np.zeros((number_of_subboxes , 1))

    for subvol_index in xrange(number_of_subboxes):

        sub_rsd_pos = mask_positions(rsd_pos , subvol_index , nsub)
        sub_rsd_pos = sub_rsd_pos.astype(np.float32)
        sub_x, sub_y, sub_z = sub_rsd_pos[:,0] , sub_rsd_pos[:,1] , sub_rsd_pos[:,2]
        wp_result = _countpairs.countpairs_wp(subbox_size, pimax, 
                                   nthreads, binfile, 
				   sub_x, sub_y, sub_z)
        nbars[subvol_index] = 1.*len(sub_rsd_pos)/subbox_size**3.
        wps[subvol_index] = np.array(wp_result)[:,3]
           
    nbar_covar = np.array([np.var(nbars)])    
    wp_covar = np.cov(wps.T)    

    return nbar_covar , wp_covar

def save_jackknife_covariance(Mr , nsub):

    wpcov_filename = util.dat_dir()+'wpcov_SHAM_Mr'+str(Mr)+'.dat'
    nbarcov_filename =  util.dat_dir()+'nbarcov_SHAM_Mr'+str(Mr)+'.dat'

    nbarcov , wpcov = compute_jackknife_covariance(Mr , nsub)

    np.savetxt(nbarcov_filename ,  nbarcov)
    np.savetxt(wpcov_filename ,  wpcov)
   
    return None 

if __name__ == "__main__":

     Mr = 19
     nsub = 3
     save_luminosity_threshold(Mr)    
     save_nbar_clustering(Mr)
     save_jackknife_covariance(Mr , nsub)   
