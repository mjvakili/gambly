from halotools.empirical_models import Zheng07Sats , HeavisideAssembias
import numpy as np
from halotools.mock_observables.catalog_analysis_helpers import return_xyz_formatted_array
from halotools.empirical_models import HodModelFactory
from halotools.empirical_models import TrivialPhaseSpace, Zheng07Cens
from halotools.sim_manager import CachedHaloCatalog
halocat = CachedHaloCatalog(simname = 'multidark', redshift = 0)

class AssembiasZheng07Sats(Zheng07Sats, HeavisideAssembias):

    def __init__(self, **kwargs):

        Zheng07Sats.__init__(self, threshold = -21)

        HeavisideAssembias.__init__(self,
            method_name_to_decorate = 'mean_occupation',
            lower_assembias_bound = 0,
            upper_assembias_bound = np.inf,
            **kwargs)



cens_occ_model =  Zheng07Cens(threshold = -21)
cens_prof_model = TrivialPhaseSpace()

from halotools.empirical_models import NFWPhaseSpace
sats_occ_model =  AssembiasZheng07Sats()
sats_prof_model = NFWPhaseSpace()

model_instance = HodModelFactory(
        centrals_occupation = cens_occ_model,
        centrals_profile = cens_prof_model,
        satellites_occupation = sats_occ_model,
        satellites_profile = sats_prof_model)

model_instance.param_dict['mean_occupation_satellites_assembias_param1'] = -1.
model_instance.populate_mock(simname = 'multidark')


x = model_instance.mock.galaxy_table['x']
y = model_instance.mock.galaxy_table['y']
z = model_instance.mock.galaxy_table['z']
vz = model_instance.mock.galaxy_table['vz']

pos = return_xyz_formatted_array(x, y, z, velocity = vz, velocity_distortion_dimension = 'z')

