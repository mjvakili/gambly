import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend("Agg")
from halotools.empirical_models import PrebuiltModelFactory
from scipy import interpolate




model = PrebuiltModelFactory("zheng07")
model.populate_mock()



