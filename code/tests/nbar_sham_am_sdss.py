import numpy as np
import matplotlib.pyplot as plt
from matplotlib import lines as mlines


def load_nbar(Mr):
    '''load the number density of the data
    '''
    if Mr == 21.:
        nbar = 1.16 * 10**-3.
    if Mr == 20.5:
        nbar = 3.13 * 10**-3.
    if Mr == 20.:
        nbar = 6.37 * 10**-3.
    if Mr == 19.5:
        nbar = 11.64 * 10**-3.
    if Mr == 19.:
        nbar = 15.66 * 10**-3.

    return nbar


def load_nbar_variance(Mr):
    '''load the variance of the number density of the data
    '''
    if Mr == 21.:
        nbarerr= 0.12 * 10**-3.
    if Mr == 20.5:
        nbarerr = 0.3 * 10**-3.
    if Mr == 20.:
        nbarerr = 0.75 * 10**-3.
    if Mr == 19.5:
        nbarerr = 1.27 * 10**-3.
    if Mr == 19.:
        nbarerr = 2.06 * 10**-3.

    return nbarerr ** 2.


def plot_nbar(style):

    if style == "SHAM":
       style = "SHAM"
    elif style == "AM":
       style = "AM"
    Mrs  = [19 , 20 , 21]
    nbar_data , nbarerr_data = [] , []
    nbar_model , nbarerr_model = [] , []
    for Mr in Mrs:

       nbar_data.append(load_nbar(Mr))
       nbarerr_data.append(load_nbar_variance(Mr)**.5)

       nbar_model.append(np.loadtxt("../../dat/nbar_"+str(style)+"_Mr"+str(Mr)+".dat"))
       nbarerr_model.append(np.loadtxt("../../dat/nbarcov_"+str(style)+"_Mr"+str(Mr)+".dat")**.5)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)

    ax.set_title(str(style)+r" number density vs SDSS number density")
    ax.set_xlabel(r"$M_{r}$" , fontsize = 20)
    ax.set_ylabel(r"number density" , fontsize = 20)
    ax.errorbar(Mrs , nbar_data , nbarerr_data , fmt = ".b")         
    ax.errorbar(Mrs , nbar_model , nbarerr_model , fmt = ".r")
    line1 = mlines.Line2D([], [], ls = '-', c = 'blue', linewidth=1.5, 
                           label='SDSS')
    line2 = mlines.Line2D([], [], ls = '-', c = 'red', linewidth=1.5, 
                           label=str(style)) 
    ax.set_xlim([18. , 22.]) 
    plt.legend(handles=[line1 , line2], frameon=False, loc='best', fontsize=20)
     
    plt.savefig("../figs/"+str(style)+"_data_numberdensity.pdf")


if __name__ == "__main__":

   styles = ["SHAM","AM"]
   for style in styles:
       plot_nbar(style)         
