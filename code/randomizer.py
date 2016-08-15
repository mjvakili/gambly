import yao_shuffler
import h5py
import numpy as np

def load_galcat(Mr , filename):
    
    cat = h5py.File(filename , "r")["data"]
    mrfilter = np.where(cat["mag_r"][:]>Mr)[0]
    mock_ids = cat["id"][:][mrfilter]
    return mock_ids


def load_halocat(filename):

    halocat = h5py.File(filename , "r")["data"]
    length_file = len(halocat['x'][:]) 
    yao_columns = ["id" , "upid" , "pid" , "x" , "y" , "z", "vz" , "mvir"]
    cat_array = np.recarray((length_file,) ,  [("id" , "i8"), ("upid" , "i8") , ("pid" , "i8") , ("x" , "f4") , ("y" ,"f4" ) , ("z" , "f4"), ("vz" , "f4") , ("mvir" , "f4")])
    for col in yao_columns:
        cat_array[col] = halocat[col][:].copy()
    generator = yao_shuffler.generate_upid
    upid = generator(halocat['pid'] , halocat['id'])
    cat_array["upid"] = upid
    return cat_array


if __name__ == '__main__':

   shuffler = yao_shuffler.shuffleMockCatalog
   Mr = 18.
   gal_file = "/export/bbq2/mj/bolshoi_a1.00231.mag_r.source_blanton.scatter0.17.tailored.hdf5"
   mock_ids = load_galcat(Mr , gal_file)
   halo_file = "/export/bbq2/mj/bolshoi_a1.00231.hdf5"
   halocat = load_halocat(halo_file)
   xyz = shuffler(mock_ids , halocat , bin_width = 0.1 , proxy = "mvir" , box_size = 250)

    
