'''
Plotting modules
'''
from halotools.sim_manager import CachedHaloCatalog
import os.path as path

from Corrfunc import _countpairs
from Corrfunc.utils import read_catalog
from halotools.empirical_models.factories.mock_helpers import three_dim_pos_bundle
from halotools.mock_observables.catalog_analysis_helpers import return_xyz_formatted_array
from halotools.empirical_models import enforce_periodicity_of_box
from halotools.empirical_models import HodModelFactory
from halotools.empirical_models import TrivialPhaseSpace, AssembiasZheng07Cens
from halotools.empirical_models import AssembiasZheng07Sats , HeavisideAssembias
from halotools.empirical_models import PrebuiltHodModelFactory
from halotools.mock_observables import tpcf
from halotools.empirical_models.factories.mock_helpers import three_dim_pos_bundle
from halotools.mock_observables import FoFGroups
from halotools.mock_observables.pair_counters import npairs_3d
from halotools.mock_observables.catalog_analysis_helpers import return_xyz_formatted_array
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors
from Corrfunc.utils import read_catalog
from halotools.empirical_models import NFWPhaseSpace
from matplotlib import lines as mlines

halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')
import os
import os.path as path 

import h5py
import corner
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import colorConverter
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors
plt.switch_backend("Agg")

def composite_model(Mr):

    cens_occ_model =  AssembiasZheng07Cens(threshold = -1.*Mr)
    cens_prof_model = TrivialPhaseSpace()
    sats_occ_model =  AssembiasZheng07Sats(threshold = -1.*Mr)
    sats_prof_model = NFWPhaseSpace()
    
    return HodModelFactory(
               centrals_occupation = cens_occ_model,
               centrals_profile = cens_prof_model,
               satellites_occupation = sats_occ_model,
               satellites_profile = sats_prof_model)

class MCMC_model(object):

    def __init__(self, Mr):
        
        self.Mr = Mr
        self.model = composite_model(Mr)
        self.halocat = CachedHaloCatalog(simname = 'bolplanck', redshift = 0, halo_finder = 'rockstar')

        ###pair counter settings ###

        self.boxsize = self.halocat.Lbox
        self.nthreads = 1
        self.pimax = 40.0
        self.binfile = path.join(path.dirname(path.abspath(__file__)),
                        "../", "bin")
        self.autocorr = 1

    def __call__(self, theta, prior_range):
        return self._sum_stat(theta, prior_range=prior_range)

    def _sum_stat(self, theta, prior_range=None):
        
        self.model.param_dict['logM0'] = theta[0]
        self.model.param_dict['sigma_logM'] = theta[1]
        self.model.param_dict['logMmin'] = theta[2]
        self.model.param_dict['alpha'] = theta[3]
        self.model.param_dict['logM1'] = theta[4]
        self.model.param_dict['mean_occupation_centrals_assembias_param1']= theta[5]
        self.model.param_dict['mean_occupation_satellites_assembias_param1']= theta[6]

        self.model.populate_mock(self.halocat) 
        x = self.model.mock.galaxy_table['x']
        y = self.model.mock.galaxy_table['y']
        z = self.model.mock.galaxy_table['z']
        vz = self.model.mock.galaxy_table['vz']
        # applying RSD
        pos = return_xyz_formatted_array(x, y, z, velocity = vz, velocity_distortion_dimension = 'z')
        # enforcing PBC
        pos = enforce_periodicity_of_box(pos, self.boxsize)
        pos = pos.astype(np.float32)
        x, y, z = pos[:,0] , pos[:,1] , pos[:,2]
        results_wp = _countpairs.countpairs_wp(self.boxsize, self.pimax, 
                                           self.nthreads,
                                           self.binfile, x, y, z)
        wp = np.array(results_wp)[:,3]
        return wp

prior_min = [10., 0.05, 10., 0.85, 12., -1. , -1.]
prior_max = [14.5, 1.5, 14., 1.45, 15., 1. , 1.]  

def plot_time_mcmc(Nwalkers, Nchains, filename):

    sample = h5py.File(filename , "r")["mcmc"]
    npars = sample.shape[2]
    fig , axes = plt.subplots(npars , 1 , sharex=True, figsize=(10, 12))

    for i in xrange(npars):
        axes[i].plot(sample[:, :, i], color="b", alpha=.4 , lw = .5)
	axes[i].yaxis.set_major_locator(MaxNLocator(5))
        axes[i].set_ylim([prior_min[i], prior_max[i]])
        axes[i].set_xlim(0, Nchains)
        #axes[i].set_ylabel(labels[i], fontsize=25)

    axes[-1].set_xlabel("Step Number", fontsize=25)
    fig.tight_layout(h_pad=0.0)
    fig_file = "mcmc_time.pdf"
    plt.savefig(fig_file)
    plt.close()

def plot_overlay_corner(Nchains1 , Nchains2, Nburns1, Nburns2, Mr , style , filename1, filename2):

    from truth import truth

    

    sample1 = h5py.File(filename1 , "r")["mcmc"]
    sample1 = sample1[Nburns1:Nchains1, : , :]
    sample1 = sample1.reshape(sample1.shape[0]*sample1.shape[1] , sample1.shape[2])
   
    sample2 = h5py.File(filename2 , "r")["mcmc"]
    sample2 = sample2[Nburns2:Nchains2, : , :]
    sample2 = sample2.reshape(sample2.shape[0]*sample2.shape[1] , sample2.shape[2])

    sample3 = np.zeros((sample2.shape[0] , sample1.shape[1]))
    sample3[:,0:5] = sample2
    sample3[:,5:7] += 110 
    prior_range = np.zeros((len(prior_min),2))
    prior_range[:,0] = np.array(prior_min)
    prior_range[:,1] = np.array(prior_max) 
    prior_range2 = prior_range.copy()
    prior_range2[5,0] = 100
    prior_range2[5,1] = 120
    prior_range2[6,0] = 100
    prior_range2[6,1] = 120
    
    fig = corner.corner(
            sample3,
            labels=[
                r'$\log\;\mathcal{M}_{0}}$',
                r'$\sigma_\mathtt{log\;M}}$',
                r'$\log\;\mathcal{M}_\mathtt{min}}$',
                r'$\alpha$',
                r'$\log\;\mathcal{M}_{1}}$',
                r'$\mathcal{A}_{\rm cen}}$',
                r'$\mathcal{A}_{\rm sat}}$'
                ],
            label_kwargs={'fontsize': 25},
            range=prior_range2,
            quantiles=[0.5],
            title_args={"fontsize": 12},
            plot_datapoints=False,
            fill_contours=True,
            levels=[0.68, 0.95],
            color='#FF7F0E' ,
            scale_hist = False,
            bins=20,
            smooth = 1.)

    corner.corner(
            sample1,
            labels=[
                r'$\log\;\mathcal{M}_{0}}$',
                r'$\sigma_\mathtt{log\;M}}$',
                r'$\log\;\mathcal{M}_\mathtt{min}}$',
                r'$\alpha$',
                r'$\log\;\mathcal{M}_{1}}$',
                r'$\mathcal{A}_{\rm cen}}$',
                r'$\mathcal{A}_{\rm sat}}$'
                ],
            label_kwargs={'fontsize': 25},
            range=prior_range,
            quantiles=[0.5],
            title_args={"fontsize": 12},
            plot_datapoints=False,
            fill_contours=True,
            levels=[0.68, 0.95],
            color = '#1F77B4',
            scale_hist = False,
            bins=20,
            smooth=1. , fig = fig)
 
    import matplotlib.lines as mlines

    blue_line = mlines.Line2D([], [], color='#1F77B4' , label=r'$Heaviside \; \; AB$')
    red_line = mlines.Line2D([], [], color='#FF7F0E' , label=r'$Standard \; \; HOD$')

    #plt.legend(handles=[blue_line,red_line], bbox_to_anchor=(0., 2.0, 2., .0), loc= 5 )


    fig_name = ''.join(['post',
         str(Mr),str(style),  
        '.pdf'])
    fig.savefig(fig_name, bbox_inches='tight', dpi=150) 
    plt.close()
    return None
 
def plot_corner_mcmc(Nchains , Nburns, Mr , style , filename):

    sample = h5py.File(filename , "r")["mcmc"]
    npars = sample.shape[2]
    nwalkers = sample.shape[1]
    sample = sample[Nburns:Nchains, : , :]
    sample = sample.reshape(sample.shape[0]*sample.shape[1] , sample.shape[2])
    print sample.shape    
    #prior_min, prior_max = PriorRange('first_try' , Mr)
    prior_range = np.zeros((len(prior_min),2))
    prior_range[:,0] = np.array(prior_min)
    prior_range[:,1] = np.array(prior_max) 

    fig = corner.corner(
            sample,
            labels=[
                r'$\log\;\mathcal{M}_{0}}$',
                r'$\sigma_\mathtt{log\;M}}$',
                r'$\log\;\mathcal{M}_\mathtt{min}}$',
                r'$\alpha$',
                r'$\log\;\mathcal{M}_{1}}$',
                r'$\mathcal{A}_{\rm cen}}$',
                r'$\mathcal{A}_{\rm sat}}$',
                ],
            label_kwargs={'fontsize': 25},
            range=prior_range,
            quantiles=[0.16,0.5,0.84],
            show_titles=True,
            title_args={"fontsize": 12},
            plot_datapoints=False,
            fill_contours=True,
            levels=[0.68, 0.95],
            color='#ee6a50',
            bins=25,
            smooth=1.)
    fig_name = ''.join(['post',
         str(Mr),str(style),  
        '.pdf'])
    fig.savefig(fig_name, bbox_inches='tight', dpi=150) 
    plt.close()
    return None 

def plot_corner_mcmc_hod(Nchains , Nburns, Mr , style , filename):

    sample = h5py.File(filename , "r")["mcmc"]
    npars = sample.shape[2]
    nwalkers = sample.shape[1]
    sample = sample[Nburns:Nchains, : , :]
    sample = sample.reshape(sample.shape[0]*sample.shape[1] , sample.shape[2])
    print sample.shape    
    #prior_min, prior_max = PriorRange('first_try' , Mr)
    prior_range = np.zeros((len(prior_min)-2,2))
    prior_range[:,0] = np.array(prior_min)[:-2]
    prior_range[:,1] = np.array(prior_max)[:-2] 

    fig = corner.corner(
            sample,
            labels=[
                r'$\log\;\mathcal{M}_{0}}$',
                r'$\sigma_\mathtt{log\;M}}$',
                r'$\log\;\mathcal{M}_\mathtt{min}}$',
                r'$\alpha$',
                r'$\log\;\mathcal{M}_{1}}$'
                ],
            label_kwargs={'fontsize': 25},
            range=prior_range,
            quantiles=[0.16,0.5,0.84],
            show_titles=True,
            title_args={"fontsize": 12},
            plot_datapoints=False,
            fill_contours=True,
            levels=[0.68, 0.95],
            color='#ee6a50',
            bins=25,
            smooth=1.)
    fig_name = ''.join(['post',
         str(Mr),str(style),  
        '.pdf'])
    fig.savefig(fig_name, bbox_inches='tight', dpi=150) 
    plt.close()
    return None 


def plot_predictions(Mr, nburnins, nchains, assembly = True, clotter = False):

    if (assembly == True):
        filename = 'Mr'+str(Mr)+'.hdf5'
    else:
        filename = 'adhoc_Mr'+str(Mr)+'.hdf5'

    sample = h5py.File(filename , "r")["mcmc"]
    npars = sample.shape[2]
    nwalkers = sample.shape[1]
    sample = sample[nchains-nburnins:nchains, : , :]
    sample = sample.reshape(sample.shape[0]*sample.shape[1] , sample.shape[2])
    print np.percentile(sample , [16,50,84] , axis = 0)
    if (clotter==False):
        model = MCMC_model(Mr) 
    
        model_wp = []
        for i in xrange(len(sample) - nwalkers*5 , len(sample)-1):
            print i
            model_wp.append(model._sum_stat(sample[i] , prior_range = None))
            np.savetxt("model_wp_assembly_"+str(Mr)+".dat" , np.array(model_wp)) 
    else:

        theta_best = np.median(sample, axis = 0)
        print theta_best
        model = MCMC_model(Mr)
        wp_best = model._sum_stat(theta_best , prior_range = None)
        if Mr == 19.0:
           cov = np.loadtxt("../../data/wpxicov_dr72_bright0_mr19.0_z0.064_nj400")[:12,:12]
           data_wp = np.loadtxt("../../data/wpxi_dr72_bright0_mr19.0_z0.064_nj400")[:,1]
        if Mr == 20.0:
           cov = np.loadtxt("../../data/wpxicov_dr72_bright0_mr20.0_z0.106_nj400")[:12,:12]
           data_wp = np.loadtxt("../../data/wpxi_dr72_bright0_mr20.0_z0.106_nj400")[:,1]
        prettyplot()
        pretty_colors=prettycolors()
        fig = plt.figure(1, figsize=(16,12))
        ax = fig.add_subplot(111)
        wps = np.loadtxt("model_wp_assembly_"+str(Mr)+".dat")
        a, b, c, d, e = np.percentile(wps, [2.5, 16, 50, 84, 97.5], axis=0)
        rbin = np.loadtxt(path.join(path.dirname(path.abspath(__file__)),
                        "../", "bin"))
        rbin_center = np.mean(rbin , axis = 1)
        
        ax.fill_between(rbin_center, a, e, color=pretty_colors[3], alpha=0.3, edgecolor="none")
        ax.fill_between(rbin_center, b, d, color=pretty_colors[3], alpha=0.5, edgecolor="none")
        ax.errorbar(rbin_center, data_wp, np.diag(cov)**.5, markersize=0, lw=0, capsize=3, elinewidth=1.5) 
        #ax.plot(rbin_center , wp_best , color = "red", lw=1.5) 
        ax.set_xlabel(r'$r_{p}[\mathtt{Mpc}]$', fontsize=27)
        ax.set_ylabel(r'$w_{p}(\mathtt{r_{p}})$', fontsize=27)
        ax.set_title(r"posterior prediction of $w_p$ for $\mathrm{M_r}$<-"+str(Mr))
        ax.set_yscale('log') 
        ax.set_xscale('log')
        ax.set_xticklabels([])
        ax.set_xlim([0.05, 25.])
        #ax.set_ylim([0.09, 1000.])
        fig.savefig("posterior_prediction"+str(Mr)+".pdf", bbox_inches='tight')
        plt.close()
        return None

def plot_occupations(Mr, nburnins, nchains, assembly = True , clotter = False , style = 'sdss'):

    model = composite_model(Mr)
    npts = 1e3
    mass = np.logspace(11, 14, npts)
    prettyplot()
    pretty_colors=prettycolors()
    fig = plt.figure(1, figsize=(16,12))

    if (assembly == True):
        filename = 'mcmc_chain_Mr'+str(Mr)+'_style_'+str(style)+'.hdf5'
    else:
        filename = 'adhoc_Mr'+str(Mr)+'.hdf5'

    sample = h5py.File(filename , "r")["mcmc"]
    npars = sample.shape[2]
    nwalkers = sample.shape[1]
    sample = sample[nchains-2:nchains, : , :]
    print sample.shape 
    sample = sample.reshape(sample.shape[0]*sample.shape[1] , sample.shape[2])
    if (clotter == False):
    	nsat_old = []
   	nsat_young = []
    	ncen_old = []
    	ncen_young = []

    	for i in xrange(len(sample)):
        	print i     
        	model.param_dict['logM0'] =  sample[i][0]
		model.param_dict['sigma_logM'] =  sample[i][1]
		model.param_dict['logMmin'] =  sample[i][2]
		model.param_dict['alpha'] =  sample[i][3]
		model.param_dict['logM1'] =  sample[i][4]
        	model.param_dict['mean_occupation_centrals_assembias_param1'] = sample[i][5]
        	model.param_dict['mean_occupation_satellites_assembias_param1'] = sample[i][6]
        	ncen_old.append(model.mean_occupation_centrals(prim_haloprop = mass, sec_haloprop_percentile=0))
        	ncen_young.append(model.mean_occupation_centrals(prim_haloprop = mass, sec_haloprop_percentile=1))
        	nsat_old.append(model.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=0))
        	nsat_young.append(model.mean_occupation_satellites(prim_haloprop = mass, sec_haloprop_percentile=1))
    	nsat_old = np.array(nsat_old)
    	nsat_young = np.array(nsat_young)
    	ncen_old = np.array(ncen_old)
    	ncen_young = np.array(ncen_young)
    	np.savetxt("nsat_old"+str(Mr)+"_"+str(style)+".dat" , nsat_old)
    	np.savetxt("nsat_young"+str(Mr)+"_"+str(style)+".dat" , nsat_young)
    	np.savetxt("ncen_old"+str(Mr)+"_"+str(style)+".dat" , ncen_old)
    	np.savetxt("ncen_young"+str(Mr)+"_"+str(style)+".dat" , ncen_young)

    else:

       nsat_old = np.loadtxt("nsat_old"+str(Mr)+"_"+str(style)+".dat")
       nsat_young = np.loadtxt("nsat_young"+str(Mr)+"_"+str(style)+".dat")
       ncen_old = np.loadtxt("ncen_old"+str(Mr)+"_"+str(style)+".dat")
       ncen_young = np.loadtxt("ncen_young"+str(Mr)+"_"+str(style)+".dat")
       
       a1, b1, c1 = np.percentile(nsat_young, [16, 50, 84], axis=0)
       a2, b2, c2 = np.percentile(nsat_old, [16, 50, 84], axis=0)
       a3, b3, c3 = np.percentile(ncen_young, [16, 50, 84], axis=0)
       a4, b4, c4 = np.percentile(ncen_old, [16, 50, 84], axis=0)
        
       fig = plt.figure()
       ax = fig.add_subplot(111)

       xlabel = ax.set_xlabel(r'$M_{\rm vir} [M_{\odot}]$', fontsize=25)
       ylabel = ax.set_ylabel(r'$\langle N_{\rm s}\rangle$', fontsize=25)
       title = ax.set_title(r'$\langle N_{\rm s} \rangle$ for $\mathrm{M_{r}}$ < -'+str(Mr)+',catalog='+str(style) , fontsize=20)

       
       ax.plot(mass, b1, color='blue', linewidth=3.5)
       ax.fill_between(mass, a1, c1 , color='blue', alpha = 0.1)
       ax.plot(mass, b2, color='red', linewidth=3.5)
       ax.fill_between(mass, a2, c2 , color='red', alpha = 0.1)
       ax.plot(mass, 0.5*(b1 + b2), '--', color='k', linewidth=2.5)

       plt.loglog()
       plt.xlim(xmin=1e11, xmax=1e14)
       plt.ylim(ymin=5e-3, ymax=100)
       plt.xticks(fontsize=20)
       plt.yticks(fontsize=20)

       blue_line = mlines.Line2D([], [], ls = '-', c = 'b', linewidth=3, label = 'high-concentration halos')
       red_line = mlines.Line2D([], [], ls = '-', c = 'r', linewidth=3, label = 'low-concentration halos')
       black_line = mlines.Line2D([], [], ls = '--', c = 'k', linewidth=3, label = 'all halos')
       first_legend = plt.legend(handles=[blue_line, red_line, black_line], frameon=False, loc='best', fontsize=17)
       fig.savefig("nsats"+str(Mr)+"_"+str(style)+".pdf", bbox_extra_artists=[xlabel, ylabel], bbox_inches='tight')
       
       fig = plt.figure()
       ax = fig.add_subplot(111)

       xlabel = ax.set_xlabel(r'$M_{\rm vir} [M_{\odot}]$', fontsize=25)
       ylabel = ax.set_ylabel(r'$\langle N_{\rm c}\rangle$', fontsize=25)
       title = ax.set_title(r'$\langle N_{\rm c} \rangle$ for $\mathrm{M_{r}}$ < -'+str(Mr)+',catalog='+str(style), fontsize=20)

       
       ax.plot(mass, b3, color='blue', linewidth=3.5)
       ax.fill_between(mass, a3, c3 , color='blue', alpha = 0.1)
       ax.plot(mass, b4, color='red', linewidth=3.5)
       ax.fill_between(mass, a4, c4 , color='red', alpha = 0.1)
       ax.plot(mass, 0.5*(b3 + b4), '--', color='k', linewidth=2.5)

       plt.loglog()
       plt.xlim(xmin=1e11, xmax=1e14)
       plt.ylim(ymin=5e-3, ymax=100)
       plt.xticks(fontsize=20)
       plt.yticks(fontsize=20)

       blue_line = mlines.Line2D([], [], ls = '-', c = 'b', linewidth=3, label = 'high-concentration halos')
       red_line = mlines.Line2D([], [], ls = '-', c = 'r', linewidth=3, label = 'low-concentration halos')
       black_line = mlines.Line2D([], [], ls = '--', c = 'k', linewidth=3, label = 'all halos')
       first_legend = plt.legend(handles=[blue_line, red_line, black_line], frameon=False, loc='best', fontsize=17)
       fig.savefig("ncens"+str(Mr)+"_"+str(style)+".pdf", bbox_extra_artists=[xlabel, ylabel], bbox_inches='tight')
    return None 

if __name__=='__main__':
 
   dire = "/export/bbq2/mj/chains/"
   #filename = "mcmc_chain_Mr19.0_style_AM.hdf5"
   #filename = "Mr19.0-f.hdf5"
   #filename = "Mr20.0-group.hdf5"

   filename1 = dire+"group_mcmc_chain_Mr20.0.hdf5"
   filename2 = dire+"adhoc_group_mcmc_chain_Mr20.0.hdf5"

   Nchains1 , Nburns1 = 3100 , 1100
   Nchains2 , Nburns2 = 2300 , 300
   
   Mr = 20.0
   style = "gmf"
   #plot_time_mcmc(Nwalkers = 140, Nchains = 2200, filename=filename)
   #plot_predictions(19.0 , 8000 , 20000, True , True)
   #plot_occupations(20. , 2000 , 2100 , True , True , "SHAM")
   #plot_corner_mcmc_hod(Nchains = 3000 , Nburns = 500, Mr = 18.5, style = "wp-hod", filename=filename)
   plot_overlay_corner(Nchains1 , Nchains2, Nburns1, Nburns2, Mr , style , filename1, filename2) 
