import corner
import matplotlib.pyplot as plt
import numpy as np

ndim, nsamples = 3, 50000

# Generate some fake data.
np.random.seed(42)
data1 = np.random.randn(ndim * 4 * nsamples // 5).reshape([4 * nsamples // 5, ndim])
data2 = (4*np.random.rand(ndim)[None, :] + np.random.randn(ndim * nsamples // 5).reshape([nsamples // 5, ndim]))
data = np.vstack([data1, data2])

import matplotlib.lines as mlines

#blue_line = mlines.Line2D([], [], color='blue', label='SGM')
#red_line = mlines.Line2D([], [], color='red', label='GMM')

data2 = data.copy()
#data2[:,2] = data2[:,2] + 10.
data[:,2] = data[:,2] * 0.
print data[:,:2] - data2[:,:2]
#plt.legend(handles=[blue_line,red_line], bbox_to_anchor=(0., 1.0, 1., .0), loc=4)

prior_range = [[-4. , 10],[-3.,3.],[-3.,3.]]

# Plot it.
figure = corner.corner(data, labels=[r"$x$", r"$y$", r"$\log \alpha$", r"$\Gamma \, [\mathrm{parsec}]$"],
                       title_kwargs={"fontsize": 12}, 
                       range=prior_range,
                       quantiles=[0.16,0.5,0.84],
                       show_titles=False,
                       title_args={"fontsize": 12},
                       plot_datapoints=False,
                       fill_contours=True,
                       levels=[0.68, 0.95],
                       color='#ee6a50',
                       bins=50,
                       smooth=1.0)


corner.corner(data2, labels=[r"$x$", r"$y$", r"$\log \alpha$", r"$\Gamma \, [\mathrm{parsec}]$"],
                       title_kwargs={"fontsize": 12}, 
                       range=prior_range,
                       quantiles=[0.16,0.5,0.84],
                       show_titles=False,
                       title_args={"fontsize": 12},
                       plot_datapoints=False,
                       fill_contours=True,
                       levels=[0.68, 0.95],
                       color='blue',
                       bins=50,
                       smooth=1.0 , fig = figure)
plt.show()
