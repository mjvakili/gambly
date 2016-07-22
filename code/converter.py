import numpy as np
import h5py


filepath = '/export/bbq2/mj/'
#################################################################
catalogues =filepath+ 'bolshoi_new_a1.00231.mag_r.source_blanton.scatter0.17.tailored'
#catalogues =filepath+ 'bolshoi_new_a1.00231.mag_r.source_blanton.scatter0.15.Vpeak'

sham_columns = ['id', 'host_mass', 'upid' , 'x' , 'y' , 'z' , 'vx' , 'vy' , 'vz' , 'mag_r' , 'mvir' , 'rvir']

def reduce_cat(Mr):

    """
    impose luminosity cut-off 
    """

    catfile = h5py.File(catalogues+'.hdf5', 'r')
    cat = catfile['data']
    length_file = len(cat['x'][:]) 
    cat_file = np.zeros((length_file , len(sham_columns)))
     
    cat_file[:,0] = cat['id'][:]
    cat_file[:,1] = cat['host_mass'][:]
    cat_file[:,2] = cat['upid'][:]
    cat_file[:,3] = cat['x'][:]
    cat_file[:,4] = cat['y'][:]
    cat_file[:,5] = cat['z'][:]
    cat_file[:,6] = cat['vx'][:]
    cat_file[:,7] = cat['vy'][:]
    cat_file[:,8] = cat['vz'][:]
    cat_file[:,9] = cat['mag_r'][:]
    cat_file[:,10] = cat['mvir'][:]
    cat_file[:,11] = cat['rvir'][:]

    #for i, col in enumerate(sham_columns): 
    #        cat_file[:,i] = f[col][:]
    #        print cat_file[:,i]
    luminosity = cat_file[:,9]
    cat_file = cat_file[luminosity > Mr] 
    catfile.close()
    return cat_file

def format_cat(Mr):

    cat_file = reduce_cat(Mr)
    ngal = cat_file.shape[0]
    dtype = [(col , float) for col in sham_columns]
    cat_array = np.recarray((ngal,)  , dtype)
    cat_array['id'] = cat_file[:,0]   
    cat_array['host_mass'] = cat_file[:,1]   
    cat_array['upid'] = cat_file[:,2]   
    cat_array['x'] = cat_file[:,3]   
    cat_array['y'] = cat_file[:,4]   
    cat_array['z'] = cat_file[:,5]   
    cat_array['vx'] = cat_file[:,6]   
    cat_array['vy'] = cat_file[:,7]   
    cat_array['vz'] = cat_file[:,8]   
    cat_array['mag_r'] = cat_file[:,9]   
    cat_array['mvir'] = cat_file[:,10]   
    cat_array['rvir'] = cat_file[:,11]   
    print cat_array
    return cat_array 
       
def save_cat(Mr):

    cat_array = format_cat(Mr)
    #print cat_array
    name = catalogues+"_Mr"+str(Mr)+".hdf5"
    filename = h5py.File(name , "w")
    filename.create_dataset("data" , data = cat_array)
    filename.close()
    
    return None        

if __name__ == "__main__":

   for mr in [18.0,18.5,19.0,19.5,20.0]:
       print "writing out reduced catalog for the Mr=" , mr 
       save_cat(mr)
        


    
     
