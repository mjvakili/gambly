import numpy as np
import matplotlib.pyplot as plt
from matplotlib import lines as mlines

def load_wp(Mr): 
    ''' loads wp 
    
    Parameters
    ----------
    Mr : (int)
        Absolute r-band magnitude threshold. Default is M_r = -21
    '''
    if Mr == 21.:
        data_file = '../../dat/wpxi_dr72_bright0_mr21.0_z0.159_nj400' 
    if Mr == 20.5:
        data_file = '../../dat/wpxi_dr72_bright0_mr20.5_z0.132_nj400'
    if Mr == 20.:
        data_file = '../../dat/wpxi_dr72_bright0_mr20.0_z0.106_nj400'
    if Mr == 19.5:
        data_file = '../../dat/wpxi_dr72_bright0_mr19.5_z0.085_nj400'
    if Mr == 19.:
        data_file = '../../dat/wpxi_dr72_bright0_mr19.0_z0.064_nj400'
    wp = np.loadtxt(data_file)[:,1]

    return wp



def load_wp_covariance(Mr): 
    ''' loads the jackknife covariance matrix associated with  wp
    
    Parameters
    ----------
    Mr : (int)
        Absolute r-band magnitude threshold. Default is M_r = -21
    '''
    if Mr == 21.:
        data_file = '../../dat/wpxicov_dr72_bright0_mr21.0_z0.159_nj400' 
    if Mr == 20.5:
        data_file = '../../dat/wpxicov_dr72_bright0_mr20.5_z0.132_nj400'
    if Mr == 20.:
        data_file = '../../dat/wpxicov_dr72_bright0_mr20.0_z0.106_nj400'
    if Mr == 19.5:
        data_file = '../../dat/wpxicov_dr72_bright0_mr19.5_z0.085_nj400'
    if Mr == 19.:
        data_file = '../../dat/wpxicov_dr72_bright0_mr19.0_z0.064_nj400'
    wpcov = np.loadtxt(data_file)[:12 , :12]

    return wpcov

def plot_wp(Mr , style):

    if style == "SHAM":
       style = "SHAM"
    elif style == "AM":
       style = "AM"

    wp_data = load_wp(Mr)
    wpcov_data = load_wp_covariance(Mr)
    wperr_data  = np.diag(wpcov_data)**.5

    wp_model = np.loadtxt("../../dat/wp_"+str(style)+"_Mr"+str(Mr)+".dat")
    wpcov_model = np.loadtxt("../../dat/wpcov_"+str(style)+"_Mr"+str(Mr)+".dat")
    print wpcov_model.shape
    wperr_model = np.diag(wpcov_model)**.5
    bins = np.loadtxt("../diagnosis/bin")
    r  =  np.mean(bins , axis = 1)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)

    ax.set_title(str(style)+r" clustering vs SDSS clustering for M$_{r}$<-"+str(Mr))
    ax.set_xlabel(r"$r$" , fontsize = 20)
    ax.set_ylabel(r"$w_{p}$" , fontsize = 20)
    ax.errorbar(r , wp_data , wperr_data , fmt = ".b")         
    ax.errorbar(r , wp_model , wperr_model , fmt = ".r")
    line1 = mlines.Line2D([], [], ls = '-', c = 'blue', linewidth=1.5, 
                           label='SDSS')
    line2 = mlines.Line2D([], [], ls = '-', c = 'red', linewidth=1.5, 
                           label=str(style)) 
   
    plt.legend(handles=[line1 , line2], frameon=False, loc='best', fontsize=20)
     
    plt.loglog()
    plt.savefig("../figs/"+str(style)+"_data"+str(Mr)+".pdf")


if __name__ == "__main__":

   Mrs = [19 , 20 , 21]
   styles = ["SHAM","AM"]
   for Mr in Mrs:
       for style in styles:
           plot_wp(Mr , style)         
