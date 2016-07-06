'''

Test that sham.py is delivering


'''
import h5py
import numpy as np
import os.path

# Local ---
import sham 
from sham import LFClass
from sham import SMFClass

# Plotting ----
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors 



def HaloMF(catalog='bolshoi', scatter=0.0, m_kind='mag_r', source='blanton'): 
    ''' 
    Compare the SHAMed LF/SMF to the analytic LF/SMFs

    Parameters
    ----------
    scatter : float
        Float that specifies the scatter in the SMHM relation, which 
        affects the SHAM masses
    source : string
        String that specifies which SMF is used for SHAM
    '''
    prettyplot()
    pretty_colors = prettycolors()

    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)
    
    # dictionary that specifies SHAM properites
    sham_dict = { 
            'm_kind': m_kind, 
            'scat': scatter, 
            'source': source, 
            'sham_prop': 'Vpeak'
            }
    haloh = sham.shamHalos(catalog=catalog, sham_dict=sham_dict)
    haloh.Read()
        
    dlogm = 0.1
    if m_kind == 'mag_r': 
        haloh.mag_r *= -1.
        m_arr = np.arange(-24.0, -16., dlogm)
    elif m_kind == 'mstar': 
        m_arr = np.arange(6.0, 12.1, dlogm)

    mag_r_bin, phi = MF_est(getattr(haloh, m_kind), 
            dlogm=dlogm, m_arr=m_arr,
            box=haloh.Lbox, h=haloh.little_h)
        
    analytic_m_kind, analytic_phi = [], [] 
    if m_kind == 'mag_r': 
        MF = LFClass(source=source, hubble=haloh.little_h)
        for mi in m_arr: 
            analytic_m_kind.append(mi - 0.5 * dlogm)
            analytic_phi.append(MF.numden(-mi, -mi+dlogm)/dlogm) 
    elif m_kind == 'mstar': 
        MF = SMFClass(source=source)
        for mi in m_arr: 
            analytic_m_kind.append(mi + 0.5 * dlogm)
            analytic_phi.append(MF.numden(mi, mi+dlogm)/dlogm) 

    sub.plot(mag_r_bin, phi, lw=2, ls='-', c=pretty_colors[3], label='SHAM-ed')

    sub.plot(analytic_m_kind, analytic_phi, lw=4, ls='--', c='k', label='Analytic')

    sub.set_yscale('log')
    sub.set_ylim([10**-7, 10**-1])
    if m_kind == 'mag_r': 
        sub.set_xlim([-17.8, -24.5])
    elif m_kind == 'mstar': 
        sub.set_xlim([8., 12.])
    sub.legend(loc='upper right')
    fig_file = ''.join(['figs/'
        'HaloMF', 
        '.', catalog, 
        '.', m_kind, 
        '.', source, 
        '.', str(round(scatter,2)),
        '.png']) 
    fig.savefig(fig_file, bbox_inches='tight')
    plt.close()
    return None


def MF_est(m_kind, dlogm=None, m_arr=None, box=None,  h=0.7):
    '''
    Calculate MF for m_kind
    '''
    if not dlogm: 
        dlogm = 0.1 # log M bin size
    if not box: 
        box = 250   # 250 Mpc/h Box length
    if m_arr is None: 
        m_arr = np.arange(0.0, 12.1, dlogm)

    vol = box ** 3  # box volume
    mass, phi = [], [] 
    for mi in m_arr: 
        mass_bin = np.where((m_kind > mi) & (m_kind <= mi+dlogm))
        ngal_bin = np.float(len(m_kind[mass_bin]))
        
        mass.append(mi + 0.5 * dlogm)
        phi.append(np.float(ngal_bin)/vol/dlogm * h**3) 

    return np.array(mass), np.array(phi)



if __name__=='__main__': 
    HaloMF(catalog='smdpl', scatter=0.0, m_kind='mag_r', source='blanton')
    HaloMF(catalog='smdpl', scatter=0.2, m_kind='mag_r', source='blanton')
