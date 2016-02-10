
from halotools.mock_observables import wp
from halotools.empirical_models.factories.mock_helpers import three_dim_pos_bundle

from matplotlib import lines as mlines
from matplotlib import pyplot as plt
plt.switch_backend("Agg")
import numpy as np

def fracdiff(alt, fid):
    return (alt-fid)/fid


rbins = np.logspace(-1, 1.25, 15)
rmax = rbins.max()
rbin_centers = (rbins[1:] + rbins[0:-1])/2.
pbins = np.logspace(-1, 1.25, 15)



def proj_clustering(pos , rbins, pbins , cellsize):

    return wp(sample1 = pos, rp_bins = rbins, 
              pi_bins = pbins, sample2 = None, 
              period = np.array([250,250,250]), approx_cell1_size = cellsize)


from halotools.empirical_models import PrebuiltHodModelFactory

# The "fiducial" stellar mass threshold is 10**10.5
model = PrebuiltHodModelFactory('hearin15', threshold = 10.5, redshift = 0., 
                        central_assembias_strength = 1, 
                        satellite_assembias_strength = 1)
model.param_dict['scatter_model_param1'] = 0.4 # This is the "fiducial" scatter used throughout the paper

baseline_model = PrebuiltHodModelFactory('leauthaud11', threshold = 10.5, redshift = 0.)
baseline_model.param_dict['scatter_model_param1'] = 0.4

#base_model

baseline_model.populate_mock()
pos = three_dim_pos_bundle(baseline_model.mock.galaxy_table, 'x', 'y', 'z')
wp_base = proj_clustering(pos , rbins , pbins , cellsize = [rmax, rmax, rmax])


#central assembly bias only

model.param_dict['mean_occupation_satellites_assembias_param1'] = 0
model.param_dict['mean_occupation_centrals_assembias_param1'] = 1
model.populate_mock()
pos = three_dim_pos_bundle(model.mock.galaxy_table, 'x', 'y', 'z')
wp_cen_only = proj_clustering(pos , rbins , pbins , cellsize = [rmax, rmax, rmax])

# Satellite assembias only

model.param_dict['mean_occupation_satellites_assembias_param1'] = 1
model.param_dict['mean_occupation_centrals_assembias_param1'] = 0
model.populate_mock()
pos = three_dim_pos_bundle(model.mock.galaxy_table, 'x', 'y', 'z')
wp_sat_only = proj_clustering(pos , rbins , pbins , cellsize = [rmax, rmax, rmax])



model.param_dict['mean_occupation_satellites_assembias_param1'] = .75
model.param_dict['mean_occupation_centrals_assembias_param1'] = 0
model.populate_mock()
pos = three_dim_pos_bundle(model.mock.galaxy_table, 'x', 'y', 'z')
wp_sat_only_75 = proj_clustering(pos , rbins , pbins , cellsize = [rmax, rmax, rmax])



model.param_dict['mean_occupation_satellites_assembias_param1'] = .5
model.param_dict['mean_occupation_centrals_assembias_param1'] = 0
model.populate_mock()
pos = three_dim_pos_bundle(model.mock.galaxy_table, 'x', 'y', 'z')
wp_sat_only_50 = proj_clustering(pos , rbins , pbins , cellsize = [rmax, rmax, rmax])


model.param_dict['mean_occupation_satellites_assembias_param1'] = .25
model.param_dict['mean_occupation_centrals_assembias_param1'] = 0
model.populate_mock()
pos = three_dim_pos_bundle(model.mock.galaxy_table, 'x', 'y', 'z')
wp_sat_only_25 = proj_clustering(pos , rbins , pbins , cellsize = [rmax, rmax, rmax])

dy_censonly = fracdiff(wp_cen_only, wp_base)
dy_satsonly = fracdiff(wp_sat_only, wp_base)
dy_satsonly75 = fracdiff(wp_sat_only_75, wp_base)
dy_satsonly50 = fracdiff(wp_sat_only_50, wp_base)
dy_satsonly25 = fracdiff(wp_sat_only_25, wp_base)

#Plotting 

fig3 = plt.figure()
ax = fig3.add_subplot(111)

xlabel = ax.set_xlabel(r'$R [\mathrm{Mpch^{-1}}]$', fontsize=20)
ylabel = ax.set_ylabel(r'$\Delta w_{\rm p} / w_{\rm p}$', fontsize=25)
#title = ax.set_title(r'$M_{\ast} > 10^{10.5} M_{\odot}$', fontsize=25)

ax.plot(rbin_centers, dy_satsonly75, '--', color='blue', linewidth=3, 
        label = r'$\mathcal{A}_{\rm bias}^{\rm cens} = 1, \mathcal{A}_{\rm bias}^{\rm sats} = 0$')
ax.plot(rbin_centers, dy_satsonly25, '--', color='black', linewidth=3, 
        label = r'$\mathcal{A}_{\rm bias}^{\rm cens} = 1, \mathcal{A}_{\rm bias}^{\rm sats} = 0$')
ax.plot(rbin_centers, dy_satsonly, '--', color='red', linewidth=3, 
        label = r'$\mathcal{A}_{\rm bias}^{\rm cens} = 0, \mathcal{A}_{\rm bias}^{\rm sats} = 1$')

ax.plot(np.logspace(-2, 2, 100), np.zeros(100), color='gray')


black_line = mlines.Line2D([], [], ls = '-', c = 'k', linewidth=3, label = '25 percent satellite galaxy assembly bias')
blue_line = mlines.Line2D([], [], ls = '-', c = 'b', linewidth=3, label = '75 percent satellite assembly bias')
red_line = mlines.Line2D([], [], ls = '-', c = 'r', linewidth=3, label = '100 percent satellite assembly bias')
first_legend = plt.legend(handles=[blue_line, red_line, black_line], frameon=False, loc='best', fontsize=17)

plt.xscale('log')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(xmin = 0.1, xmax=15)
plt.ylim(ymin = -0.1, ymax = 1)

fig3.savefig('/home/mj/public_html/assembly_plots/for_hogg/delta_wp_satsonly.png', bbox_extra_artists=[xlabel, ylabel], bbox_inches='tight')
