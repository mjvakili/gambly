'''

Code for generating the sham catalogs.

Author : ChangHoon Hahn

'''
import os 
import h5py
import numpy as np
from numpy import log10, Inf
from scipy import integrate, interpolate, ndimage


def DownloadedCatalog(catalog='bolshoi'):  
    ''' Take .list format downloaded catalog and load it into an .hdf5 file format
    so that the halo catalog can be easily read 
    ''' 
    current_dir = os.path.dirname(os.path.realpath(__file__))
    if current_dir != '/home/mj/assembly/code':
        raise ValueError("This function is only meant to be run on Chang's user directory on Sirocco!") 
    
    if catalog == 'bolshoi':        # Bolshoi Box
        list_file = ''.join([
            '/export/bbq2/mj/', 
            'hlist_1.00231.list'
            ])
        a_scale = 1.00231
        redshift = 0.0  # close enough
    elif catalog == 'smdpl':        # Small MultiDark Planck
        list_file = ''.join([
            '/export/bbq2/mj/', 
            'smdpl_hlist_1.00000.list'
            ])
        a_scale = 1.00
        redshift = 0.0  # close enough

    # read in first few lines of the file to get the columns and column indicies
    f = open(list_file, 'r') 
    first_line = f.readline() 
    col_list = first_line[1:].replace(')\n', '').replace(') ',',').split(',')   # don't touch! stirng parsing magic

    columns, column_indices = [], [] 
    for ccc in col_list: 
        columns.append('('.join(ccc.rsplit('(')[:-1]))
        column_indices.append(ccc.rsplit('(')[-1])
    
    list_data = np.loadtxt(list_file)   # this takes forever

    # save all columns to hdf5
    if catalog == 'bolshoi': 
        hdf5_file = ''.join([
            '/export/bbq2/mj/', 
            catalog, '_a1.00231.hdf5'])
    elif catalog == 'smdpl': 
        hdf5_file = ''.join([
            '/export/bbq2/mj/', 
            catalog, '_a1.00000.hdf5'])

    f = h5py.File(hdf5_file, 'w') 
    grp = f.create_group('data') 
    grp.attrs['a_scale'] = a_scale 
    grp.attrs['redshift'] = redshift
    for i_col, col in enumerate(columns): 
        if ('id' in col) or ('ID' in col): 
            data_column = list_data[:,i_col].astype('int') 
        elif col == '': 
            continue 
        else: 
            data_column = list_data[:,i_col]

        grp.create_dataset(col, data=data_column)
    
    f.close()
    return None 


class Halos(object): 
    ''' 
    Class to deal with halo catalogs downloaded from Peter Behroozi's  
    Rockstar website. 
    > http://hipacc.ucsc.edu/Bolshoi/MergerTrees.html
    '''
    
    def __init__(self, catalog='bolshoi'):
        '''
        '''
        self.column_list = ['id', 'upid', 'x', 'y', 'z', 'vx', 'vy', 'vz', 
                'Vpeak', 'Mpeak', 'vrms', 'mvir', 'rvir', 'Vmax@Mpeak', 'Mpeak_Scale']
        if catalog not in ['bolshoi', 'smdpl']:
            raise NotImplementedError("Catalog not included yet") 
        self.catalog = catalog 
   
    def File(self):
        '''
        '''
        file_dir = '/export/bbq2/mj/'
        if self.catalog == 'bolshoi': 
            self.file_name = ''.join([file_dir, 'bolshoi_a1.00231.hdf5']) 
        elif self.catalog == 'smdpl':
            self.file_name = ''.join([file_dir, 'smdpl_a1.00000.hdf5']) 
        return self.file_name
    
    def Read(self): 
        ''' Read in hdf5 file 
        '''
        file = self.File()

        f = h5py.File(file, 'r') 
        grp = f['data'] 
        for attr in grp.attrs.keys():
            setattr(self, attr, grp.attrs[attr]) 

        # import *select* columns
        cols = self.Columns()
        for col in cols:
            if col == 'Mpeak': 
                setattr(self, col, np.log10(grp[col][:]))
            elif col == 'Vmax@Mpeak':
                setattr(self, "VmaxMpeak", grp[col][:])
            else: 
                setattr(self, col, grp[col][:])
        return None 
        
    def Columns(self): 
        return self.column_list 


class shamHalos(object): 
    ''' 
    Class to deal with halo catalogs with SHAMed M_r or stellar mass 
    '''
    def __init__(self, catalog='bolshoi', sham_dict=None):  
        
        self.catalog = catalog
        self.sham_dict = sham_dict  # dictionary specifying the SHAM choices 
        self.column_list = None

    def ReadHaloCatalog(self):
        ''' Read in halo catalog using the Halos object class and import 
        the columns into object.
        '''
        halocat = Halos(catalog=self.catalog) 
        halocat.Read() 
        for col in halocat.Columns():
            if col == "Vmax@Mpeak":
                setattr(self, "VmaxMpeak", getattr(halocat, "VmaxMpeak"))
            else: 
                setattr(self, col, getattr(halocat, col)) 
        self.a_scale = halocat.a_scale
        self.redshift = halocat.redshift

        self.column_list = halocat.column_list
        return None 

    def SHAMassign(self, m_kind='mstar', scat=0, source='', sham_prop='Mpeak'):
        '''
        Assign Mag_r or M_star via abundance matching.
        '''
        self.ReadHaloCatalog()

        #self.sham_dict = {'m_kind': m_kind, 'scat': scat, 'source': source, 'sham_prop': sham_prop}
        if m_kind not in self.column_list: 
            self.column_list.append(m_kind)

        vol = (self._CatalogBox() / self._little_h()) ** 3
        print 'Box Length', self._CatalogBox()
        print 'Box Hubble', self._little_h()

        if m_kind == 'mstar':
            if not source:
                source = 'li-drory-march'
            redshift = self.z_scale 
            if redshift < 0.1:
                redshift = 0.1
            MF = SMFClass(source=source, redshift=redshift, scat=scat, hubble=self._little_h())
        elif m_kind == 'mag_r':
            if source == 'cool_ages':
                redshift = self.z_scale 
                if redshift < 0.1:
                    redshift = 0.1
                MF = LFClass(source=source, scat=scat, hubble=self._little_h(), redshift=redshift)
            else:
                if not source:
                    source = 'blanton'
                MF = LFClass(source, scat, self._little_h())
        else:
            raise ValueError('not recognize m_kind = %s' % m_kind)
    
        if sham_prop == 'tailored': 
            # special request SHAM : 
            # v_vir * (v_max / v_vir)^0.57
            sham_attr = np.zeros(len(self.vrms))
            delta  = 360
            omegam = 0.3
            omegal = 0.7 
            H_scale = omegam / (self.Mpeak_Scale)**3.  + omegal
            v_vir = (0.5 * delta * H_scale ** 2. * 4.302**2. * 10.**14. )**(1./6) * (10. ** self.Mpeak)**1./3
            #v_vir = ((100 ** 0.5) * (4.302 * 10 ** -7.) * (10. ** self.Mpeak))**1./3 
            vvir_notzero = np.where(v_vir != 0.) 
            sham_attr[vvir_notzero] = v_vir[vvir_notzero] * (
                    self.VmaxMpeak[vvir_notzero] / v_vir[vvir_notzero])**0.57
        else: 
            sham_attr = getattr(self, sham_prop)
        print 'SHAM attribute = ', sham_attr.min(), sham_attr.max()
        m_kind_attr = np.zeros(len(sham_attr), np.float32)
        if m_kind == 'mstar':
            MF.initialize_redshift(redshift)
        elif m_kind == 'mag_r':
            if source == 'cool_ages':
                MF.initialize_redshift(redshift)

        # maximum number of objects in volume to assign given SMF/LF threshold
        num_max = int(round(MF.numden(MF.mmin) * vol))
        sis = elements(sham_attr, [0.001, Inf])
        siis_sort = np.argsort(sham_attr[sis]).astype(sis.dtype)[::-1][:num_max]
        num_sums = arange_length(num_max) + 1
        if scat:
            if m_kind == 'mstar': 
                scats = np.random.normal(np.zeros(num_max), MF.scat).astype(np.float32)
            elif m_kind == 'mag_r': 
                scats = np.random.normal(np.zeros(num_max), 2.5 * MF.scat).astype(np.float32)
            #print MF.m_scat(num_sums / vol) + scats
            m_kind_attr[sis[siis_sort]] = MF.m_scat(num_sums / vol) + scats
        else:
            m_kind_attr[sis[siis_sort]] = MF.m(num_sums / vol)

        setattr(self, m_kind, m_kind_attr)
        return None 

    def File(self): 
        ''' File name of SHAMed catalog
        '''
        if self.sham_dict is None: 
            raise ValueError 
        halocat = Halos(catalog=self.catalog) 
        halo_file = halocat.File()
        halo_file = halo_file.rsplit('.hdf5')[0]

        sham_file = ''.join([
            halo_file, 
            '.', self.sham_dict['m_kind'], 
            '.source_', self.sham_dict['source'], 
            '.scatter', str(round(self.sham_dict['scat'], 2)), 
            '.', self.sham_dict['sham_prop'], 
            '.hdf5']) 
        return sham_file

    def Read(self): 
        sham_file = self.File() 
        f = h5py.File(sham_file, 'r') 
        grp = f['data'] 
        for col in grp.keys(): 
            setattr(self, col, grp[col][:]) 

        for attr in grp.attrs.keys():
            setattr(self, attr, grp.attrs[attr])
        return None

    def Write(self):
        ''' Write SHAMed halo catalog to file  
        '''
        if self.sham_dict is None: 
            raise ValueError 

        if self.column_list is None or self.sham_dict['m_kind'] not in self.column_list: 
            self.SHAMassign(
                    m_kind=sham_dict['m_kind'], 
                    scat=sham_dict['scat'], 
                    source=sham_dict['source'], 
                    sham_prop=sham_dict['sham_prop']
                    )

        sham_file = self.File() 
        f = h5py.File(sham_file, 'w') 
        grp = f.create_group('data') 

        columns = self.Columns()
        for i_col, col in enumerate(columns): 

            if col == 'M200b': 
                col_data = getattr(self, col)

                id = getattr(self, 'id') 
                upid = getattr(self, 'upid') 
                nonneg = np.where(upid != -1)[0]

                sub_id, sub_upid = intersection_index(id, upid[nonneg]) 
                col_data[nonneg[sub_upid]] = col_data[sub_id]
            elif col == 'Mpeak':
                col_data = np.log10(getattr(self, col))
            elif col == 'Vmax@Mpeak':
                col_data = getattr(self, 'VmaxMpeak')
            else: 
                col_data = getattr(self, col)

            grp.create_dataset(col, data=col_data)

        grp.attrs['a_scale'] = self.a_scale 
        grp.attrs['redshift'] = self.redshift
        grp.attrs['Lbox'] = self._CatalogBox()
        grp.attrs['little_h'] = self._little_h() 
        
        f.close()

    def Columns(self):
        return self.column_list 

    def _CatalogBox(self): 
        if self.catalog == 'multidark':
            L_box = 1000. # Mpc/h
        elif self.catalog == 'smdpl': 
            L_box = 400. # Mpc/h
        elif self.catalog == 'bolshoi': 
            L_box = 250. # Mpc/h

        return L_box

    def _little_h(self): 
        if self.catalog == 'bolshoi': 
            h_little = 0.68
        elif self.catalog == 'smdpl': 
            h_little = 0.6777
        else: 
            h_little = 0.7

        return h_little


class SMFClass:
    '''
    Relate number density [dnumden / dlog(M_star/M_sun)] <-> stellar mass [log10(M_star/M_sun)]
    using fits to observed stellar mass functions.
    All SMFs assume input Hubble constant.
    '''
    def __init__(self, source='li-march', redshift=0.1, scat=0, hubble=0.7):
        '''
        Import SMF source, redshift, log scatter in M_star at fixed Msub.
        '''
        self.source = source
        self.scat = scat
        self.hubble = hubble
        if source == 'li':
            '''
            Li & White 2009. z = 0.1 from SDSS. Chabrier IMF. Complete to 1e8 M_sun/h^2.
            '''
            self.redshifts = np.array([0.1])
            self.mchars = np.array([10.525]) - 2 * log10(hubble)    # {M_sun}
            self.amplitudes = np.array([0.0083]) * hubble ** 3    # {Mpc ^ -3 / log(M/M_sun)}
            self.slopes = np.array([-1.155])
            self.initialize_redshift(redshift)
        elif source == 'baldry':
            '''
            Baldry et al 2008. z = 0.1 from SDSS. diet Salpeter IMF = 0.7 Salpeter.
            Complete to 1e8 M_sun.
            '''
            h_them = 0.7    # their assumed hubble constant
            self.redshifts = np.array([0.1])
            # covert to Chabrier
            self.mchars = (np.array([10.525]) + 2 * log10(h_them / hubble) + log10(1 / 1.6 / 0.7))
            self.amplitudes = np.array([0.00426]) * (hubble / h_them) ** 3
            self.amplitudes2 = np.array([0.00058]) * (hubble / h_them) ** 3
            self.slopes = np.array([-0.46])
            self.slopes2 = np.array([-1.58])
            self.initialize_redshift(redshift)
        elif source == 'cole-march':
            '''
            Marchesini et al 2009. 1.3 < z < 4.0. Kroupa IMF.
            z = 0.1 from Cole et al 2001 (2dF), converting their Salpeter to Kroupa.
            *** In order to use out to z ~ 4, made evolution flat from z = 3.5 to 4.
            '''
            self.redshifts = np.array([0.1, 1.6, 2.5, 3.56, 4.03])
            self.mchars = np.array([10.65, 10.60, 10.65, 11.07, 11.07]) - 2 * log10(hubble)
            # converted to {Mpc ^ -3 dex ^ -1}
            self.amplitudes = np.array([90.00, 29.65, 11.52, 1.55, 1.55]) * 1e-4 * hubble ** 3
            self.slopes = np.array([-1.18, -1.00, -1.01, -1.39, -1.39])
            self.make_splines()
            self.initialize_redshift(redshift)
        elif source == 'li-march':
            '''
            Marchesini et al 2009, using Li & White at z = 0.1.
            '''
            self.redshifts = np.array([0.1, 1.6, 2.5, 3.56, 4.03])
            self.mchars = np.array([10.525, 10.60, 10.65, 11.07, 11.07]) - 2 * log10(hubble)
            self.amplitudes = (np.array([0.0083, 0.002965, 0.00115, 0.000155, 0.000155]) *
                               hubble ** 3)
            self.slopes = np.array([-1.155, -1.00, -1.01, -1.39, -1.39])
            self.make_splines()
            self.initialize_redshift(redshift)
        elif source == 'li-march-extreme': 
            '''
            More extreme version of Marchesini et al 2009, using Li & White at z = 0.1.
            '''
            self.redshifts = np.array([0.1, 1.6, 2.5, 3.56, 4.03])
            self.mchars = np.array([10.525, 10.60, 10.65, 11.07, 11.07]) - 2 * log10(hubble)
            self.amplitudes = (np.array([0.0083, 0.00001, 0.00001, 0.00001, 0.000001]) *
                               hubble ** 3)
            self.slopes = np.array([-1.155, -1.00, -1.01, -1.39, -1.39])
            self.make_splines()
            self.initialize_redshift(redshift)
        elif source == 'constant-li': 
            '''
            Li & White at all redshifts 
            '''
            self.redshifts = np.arange(0.1, 4.03, 0.1) 
            self.mchars = np.repeat(10.525, len(self.redshifts)) - 2 * log10(hubble)
            self.amplitudes = (np.repeat(0.0083, len(self.redshifts))* hubble ** 3)
            self.slopes = np.repeat(-1.155, len(self.redshifts))
            self.make_splines()
            self.initialize_redshift(redshift)

        elif source == 'fontana':
            '''
            Fontana et al 2006. 0.4 < z < 4 from GOODS-MUSIC. Salpeter IMF.
            z = 0.1 from Cole et al 2001.
            '''
            h_them = 0.7    # their assumed hubble constant
            self.redshifts = np.array([0.1, 4.0])    # store redshift range of validity
            self.amplitude0 = 0.0035 * (hubble / h_them) ** 3    # to {Mpc ^ -3 / log10(M/M_sun)}
            self.amplitude1 = -2.2
            self.slope0 = -1.18
            self.slope1 = -0.082
            self.mchar0 = 11.16    # log10(M/M_sun)
            self.mchar1 = 0.17    # log10(M/M_sun)
            self.mchar2 = -0.07    # log10(M/M_sun)
            # convert to my hubble & Chabrier IMF
            self.mchar0 += 2 * log10(h_them / hubble) - log10(1.6)
            self.initialize_redshift(redshift)
        elif source == 'li-drory-march':
            '''
            Drory et al 2009. 0.3 < z < 1.0 from COSMOS.
            Chabrier IMF limited to 0.1 - 100 M_sun.
            Complete to (8.0, 8.6, 8.9, 9.1) M_sun/h^2 at z = (0.3, 0.5, 0.7, 0.9).
            Anchor to Li & White at z = 0.1, Marchesini et al at higher redshift.
            See Ilbert et al 2010 for alternate COSMOS version.
            '''
            h_them = 0.72    # their assumed hubble constant
            self.redshifts = np.array([0.3, 0.5, 0.7, 0.9])
            self.mchars = np.array([10.90, 10.91, 10.95, 10.92]) + 2 * log10(h_them / hubble)
            # convert to [Mpc ^ -3 dex^-1]
            self.amplitudes = (np.array([0.00289, 0.00174, 0.00216, 0.00294]) *
                               (hubble / h_them) ** 3)
            self.slopes = np.array([-1.06, -1.05, -0.93, -0.91])
            self.mchars2 = np.array([9.63, 9.70, 9.75, 9.85]) + 2 * log10(h_them / hubble)
            self.amplitudes2 = (np.array([0.00180, 0.00143, 0.00289, 0.00212]) *
                                (hubble / h_them) ** 3)
            self.slopes2 = np.array([-1.73, -1.76, -1.65, -1.65])
            # add li & white
            self.redshifts = np.append(0.1, self.redshifts)
            self.mchars = np.append(10.525 - 2 * log10(hubble), self.mchars)
            self.amplitudes = np.append(0.0083 * hubble ** 3, self.amplitudes)
            self.slopes = np.append(-1.155, self.slopes)
            self.mchars2 = np.append(self.mchars2[0], self.mchars2)
            self.amplitudes2 = np.append(0, self.amplitudes2)
            self.slopes2 = np.append(self.slopes2[0], self.slopes2)
            # add marchesini et al
            h_them = 0.7    # their assumed hubble constant
            self.redshifts = np.append(self.redshifts, [1.6, 2.5, 3.56, 4.03])
            self.mchars = np.append(self.mchars,
                np.array([10.60, 10.65, 11.07, 11.07]) - 2 * log10(hubble))
            self.amplitudes = np.append(self.amplitudes,
                                        np.array([0.002965, 0.00115, 0.000155, 0.000155]) *
                                        hubble ** 3)
            self.slopes = np.append(self.slopes, [-1.00, -1.01, -1.39, -1.39])
            self.mchars2 = np.append(self.mchars2, np.zeros(4) + self.mchars2[0])
            self.amplitudes2 = np.append(self.amplitudes2, np.zeros(4))
            self.slopes2 = np.append(self.slopes2, np.zeros(4) + self.slopes2[0])
            self.make_splines()
            self.initialize_redshift(redshift)
        elif source == 'li-drory-march_sameslope':
            '''
            Apply low-mass slope from Drory et al 2009 to Li & White, Marchesini et al.
            '''
            self.redshifts = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.6, 2.5, 3.56, 4.03])
            self.mchars = np.array([10.525, 10.61, 10.62, 10.66, 10.63, 10.60, 10.65, 11.07,
                                    11.07] - 2 * log10(hubble))
            self.amplitudes = np.array([0.0083, 0.00774, 0.00466, 0.00579, 0.00787, 0.00297,
                                        0.00115, 0.000155, 0.000155]) * hubble ** 3
            self.slopes = np.array([-1.155, -1.06, -1.05, -0.93, -0.91, -1.00, -1.01, -1.39, -1.39])
            self.mchars2 = (np.array([9.35, 9.34, 9.41, 9.46, 9.56, 9.41, 9.46, 9.83, 9.83]) -
                            2 * log10(hubble))
            self.amplitudes2 = np.array([0.00269, 0.00482, 0.00383, 0.00774, 0.00568, 0.000962,
                                         0.000375, 0.0000503, 0.0000503]) * hubble ** 3
            self.slopes2 = np.array([-1.70, -1.73, -1.76, -1.65, -1.65, -1.72, -1.74, -2.39, -2.39])
            self.make_splines()
            self.initialize_redshift(redshift)
        elif source == 'perez':
            '''
            Perez-Gonzalez et al 2008. 0.1 < z < 4.0 from Spitzer, Hubble, Chandra.
            Salpeter IMF.
            Complete to (8, 9.5, 10, 11) M_star at z = (0, 1, 2, 3).
            '''
            h_them = 0.7    # their assumed hubble constant
            self.redshifts = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.15, 1.45, 1.8, 2.25, 2.75, 3.25,
                                       3.75])
            self.mchars = np.array([11.16, 11.20, 11.26, 11.25, 11.27, 11.31, 11.34, 11.40, 11.46,
                                    11.34, 11.33, 11.36]) + 2 * log10(h_them / hubble)
            # convert to Chabrier IMF
            self.mchars -= log10(1.6)
            # convert to [Mpc ^ -3 dex ^ -1]
            self.amplitudes = (10 ** np.array([-2.47, -2.65, -2.76, -2.82, -2.91, -3.06, -3.27,
                                              - 3.49, -3.69, -3.64, -3.74, -3.94]) *
                               (hubble / h_them) ** 3)
            self.slopes = np.array([-1.18, -1.19, -1.22, -1.26, -1.23, -1.26, -1.29, -1.27, -1.26,
                                    - 1.20, -1.14, -1.23])
            self.make_splines()
            self.initialize_redshift(redshift)
        else:
            raise ValueError('not recognize source = %s' % source)

    def make_splines(self):
        '''
        Make spline fits to SMF fit parameters v redshift.
        Use 1st order spline (k) to avoid ringing.
        '''
        self.mchar_z_spl = interpolate.splrep(self.redshifts, self.mchars, k=1)
        self.slope_z_spl = interpolate.splrep(self.redshifts, self.slopes, k=1)
        self.amplitude_z_spl = interpolate.splrep(self.redshifts, self.amplitudes, k=1)
        if self.source in ('li-drory-march', 'li-drory-march_sameslope'):
            self.mchar2_z_spl = interpolate.splrep(self.redshifts, self.mchars2, k=1)
            self.slope2_z_spl = interpolate.splrep(self.redshifts, self.slopes2, k=1)
            self.amplitude2_z_spl = interpolate.splrep(self.redshifts, self.amplitudes2, k=1)

    def initialize_redshift(self, redshift=0.1):
        '''
        Make spline to get mass from number density.

        Import redshift.
        Find SMF fit parameters at redshift, correcting amplitude by * log(10) & slope
        by + 1 to make dndm call faster.
        '''
        if redshift < self.redshifts.min() - 1e-5 or redshift > self.redshifts.max() + 1e-5:
            raise ValueError('z = %.2f out of range for %s' % (redshift, self.source))
        self.redshift = redshift
        if self.source in ('li'):
            self.m_char = self.mchars[0]
            self.amplitude = self.amplitudes[0] * np.log(10)
            self.slope = self.slopes[0] + 1
        elif self.source in ('baldry'):
            self.m_char = self.mchars[0]
            self.mchar2 = self.mchars[0]
            self.amplitude = self.amplitudes[0] * np.log(10)
            self.amplitude2 = self.amplitudes2[0] * np.log(10)
            self.slope = self.slopes[0] + 1
            self.slope2 = self.slopes2[0] + 1
        elif self.source in ('cole-march', 'li-march', 'perez', 'constant-li', 'li-march-extreme'):
            self.m_char = interpolate.splev(redshift, self.mchar_z_spl)
            self.amplitude = interpolate.splev(redshift, self.amplitude_z_spl) * np.log(10)
            self.slope = interpolate.splev(redshift, self.slope_z_spl) + 1
        elif self.source == 'fontana':
            self.m_char = self.mchar0 + self.mchar1 * redshift + self.mchar2 * redshift ** 2
            self.amplitude = (self.amplitude0 * (1 + redshift) ** self.amplitude1) * np.log(10)
            self.slope = (self.slope0 + self.slope1 * redshift) + 1
        elif self.source in ('li-drory-march', 'li-drory-march_sameslope'):
            self.m_char = interpolate.splev(redshift, self.mchar_z_spl)
            self.amplitude = interpolate.splev(redshift, self.amplitude_z_spl) * np.log(10)
            self.slope = interpolate.splev(redshift, self.slope_z_spl) + 1
            self.mchar2 = interpolate.splev(redshift, self.mchar2_z_spl)
            self.amplitude2 = interpolate.splev(redshift, self.amplitude2_z_spl) * np.log(10)
            self.slope2 = interpolate.splev(redshift, self.slope2_z_spl) + 1
        self.make_numden_m_spline(self.redshift, self.scat)

    def dndm(self, m_star):
        '''
        Compute d(num-den) / d(log m) = ln(10) * amplitude * (10^(m_star - m_char)) ** (1 + slope) *
        exp(-10^(m_star - m_char)).

        Import stellar mass.
        '''
        m_rats = 10 ** (m_star - self.m_char)
        if 'drory' in self.source or self.source == 'baldry':
            dm2s = 10 ** (m_star - self.mchar2)
            return (self.amplitude * m_rats ** self.slope * np.exp(-m_rats) +
                    self.amplitude2 * dm2s ** self.slope2 * np.exp(-dm2s))
        else:
            return self.amplitude * m_rats ** self.slope * np.exp(-m_rats)

    def numden(self, m_min, m_max=14):
        '''
        Compute number density within range.

        Import stellar mass range.
        '''
        return integrate.quad(self.dndm, m_min, m_max)[0]

    def make_numden_m_spline(self, redshift=0.1, scat=0):
        '''
        Make splines to relate d(num-den) / d[log]m & num-den(> m) to m.

        Import redshift (if want to change), mass scatter [dex].
        '''
        iter_num = 30

        if redshift != self.redshift:
            self.initialize_redshift(redshift)
        if scat != self.scat:
            self.scat = scat
        dm = 0.01
        dm_scat_lo = 3 * scat    # extend fit for deconvolute b.c.'s
        dm_scat_hi = 0.5 * scat    # extend fit for deconvolute b.c.'s
        self.mmin = 7.3
        self.mmax = 12.3
        m_stars = np.arange(self.mmin - dm_scat_lo, self.mmax + dm_scat_hi, dm, np.float32)
        numdens = np.zeros(m_stars.size)
        dndms = np.zeros(m_stars.size)
        for mi in xrange(m_stars.size):
            # make sure numdens are monotonically decreasing even if = -infinity
            numdens[mi] = self.numden(m_stars[mi]) + 1e-9 * (1 - mi * 0.001)
            dndms[mi] = self.dndm(m_stars[mi]) + 1e-9 * (1 - mi * 0.001)
        # make no scatter splines
        self.log_numden_m_spl = interpolate.splrep(m_stars, log10(numdens))
        self.m_log_numden_spl = interpolate.splrep(log10(numdens)[::-1], m_stars[::-1])
        # at high z, smf not monotonically decreasing, so spline not work on below
        # self.m_log_dndm_spl = interpolate.splrep(log10(dndms)[::-1], m_stars[::-1])
        # make scatter splines
        if scat:
            # deconvolve osbserved smf assuming scatter to find unscattered one
            dndms_scat = deconvolute(dndms, scat, dm, iter_num)
            # chop off lower boundaries, unreliable
            m_stars = m_stars[dm_scat_lo / dm:]
            dndms_scat = dndms_scat[dm_scat_lo / dm:]
            # find spline to integrate over
            self.dndm_m_scat_spl = interpolate.splrep(m_stars, dndms_scat)
            numdens_scat = np.zeros(m_stars.size)
            for mi in xrange(m_stars.size):
                numdens_scat[mi] = interpolate.splint(m_stars[mi], m_stars.max(),
                                                      self.dndm_m_scat_spl)
                numdens_scat[mi] += 1e-9 * (1 - mi * 0.001)
            self.log_numden_m_scat_spl = interpolate.splrep(m_stars, log10(numdens_scat))
            self.m_log_numden_scat_spl = interpolate.splrep(log10(numdens_scat)[::-1],
                                                            m_stars[::-1])

    def m(self, num_den):
        '''
        Get mass at threshold.

        Import threshold number density.
        '''
        return interpolate.splev(log10(num_den), self.m_log_numden_spl).astype(np.float32)

    def m_scat(self, num_den):
        '''
        Get mass at threshold, using de-scattered source.

        Import threshold number density.
        '''
        return interpolate.splev(log10(num_den), self.m_log_numden_scat_spl).astype(np.float32)

    def m_dndm(self, dn_dm):
        '''
        Get mass at d(num-den)/d[log]m.

        Import d(num-den) / d[log]m.
        '''
        return interpolate.splev(log10(dn_dm), self.m_log_dndm_spl)

    def dndm_scat(self, m):
        '''
        Get d(num-den) / d[log]m at m, using de-scattered source.

        Import mass.
        '''
        return interpolate.splev(m, self.dndm_m_scat_spl)

    def numden_scat(self, m_min, m_max=14):
        '''
        Get num-den(>[log]m) at m, using de-scattered source.

        Import mass.
        '''
        return integrate.quad(self.dndm_scat, m_min, m_max)[0]
        #return 10 ** (interpolate.splev(m, self.log_numden_m_scat_spl))


class LFClass(SMFClass):
    '''
    Relate number density [Mpc ^ -3] <-> magnitude/luminosity using spline fit to luminosity
    functions.

    Import spline querying functions from SMFClass.
    '''
    def __init__(self, source='blanton', scat=0, hubble=0.7, redshift=0.1):
        '''
        Import source, log-normal scatter.
        '''
        self.source = source
        self.scat = scat
        self.hubble = hubble
        if source == 'norberg':
            # Norberg et al 2002: 2dF r-band at z ~ 0.1.
            self.m_char = -19.66
            self.amplitude = 1.61e-2 * hubble ** 3    # Mpc ^ -3
            self.slope = -1.21
        elif source == 'blanton':
            # Blanton et al 03: SDSS r-band z ~ 0.1.
            self.m_char = -20.44
            self.amplitude = 1.49e-2 * hubble ** 3    # Mpc ^ -3
            self.slope = -1.05
        elif source == 'sheldon':
            # Sheldon et al 07: SDSS i-band z = 0.25. Valid for Mag < -19.08 (0.19L*).
            self.m_char = -20.9    # Hansen et al 09 catalog has -20.8
            self.amplitude = 1.02e-2 * hubble ** 3    # Mpc ^ -3
            self.slope = -1.21
        elif source == 'cool_ages': 
            # Cool et al 2012: AGES.
            self.redshifts = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.65])
            self.mchars = np.array([-20.58, -20.81, -20.81, -20.99, -21.29, -21.38]) 
            self.amplitudes = (np.array([1.59e-2, 1.52e-2, 1.24e-2, 1.44e-2, 1.08e-2, 1.05e-2]) * hubble ** 3)    # Mpc ^ -3
            self.slopes = np.repeat(-1.05, len(self.redshifts)) 
            self.make_splines()
            self.initialize_redshift(redshift)
        else:
            raise ValueError('not recognize source = %s in LFClass' % source)

        if source != 'cool_ages': 
            self.make_numden_m_spline(scat, redshift=None)

    def dndm(self, mag):
        '''
        Get d(num-den) / d(mag).

        Import (positive) magnitude.
        '''
        mag *= -1.
        return (np.log(10) / 2.5 * self.amplitude *
                10 ** ((self.slope + 1) / 2.5 * (self.m_char - mag)) *
                np.exp(-10 ** ((self.m_char - mag) / 2.5)))

    def numden(self, m_min, m_max=25):
        '''
        Get number density within range.

        Import (positive) magnitude range.
        '''
        return integrate.quad(self.dndm, m_min, m_max)[0]
    
    def initialize_redshift(self, redshift=0.1):
        '''
        Make spline to get mass from number density.

        Import redshift.
        Find SMF fit parameters at redshift, correcting amplitude by * log(10) & slope
        by + 1 to make dndm call faster.
        '''
        if redshift < self.redshifts.min() - 1e-5:# or redshift > self.redshifts.max() + 1e-5:
            raise ValueError('z = %.2f out of range for %s' % (redshift, self.source))
        self.redshift = redshift
        self.m_char = interpolate.splev(redshift, self.mchar_z_spl, ext=0)
        self.amplitude = interpolate.splev(redshift, self.amplitude_z_spl, ext=0) 
        self.slope = interpolate.splev(redshift, self.slope_z_spl, ext=0)
        self.make_numden_m_spline(scat = self.scat, redshift = self.redshift)

    def make_numden_m_spline(self, scat=0, redshift=0.1):
        '''
        Make splines to relate d(num-den)/d(mag) & num-den(> mag) to mag.

        Import scatter [dex].
        '''
        try: 
            if redshift != self.redshift:
                self.initialize_redshift(redshift)
        except AttributeError:
            pass 
        if scat != self.scat:
            self.scat = scat    # convert scatter in log(lum) to scatter in magnitude
        mag_scat = 2.5 * self.scat
        deconvol_iter_num = 30
        dmag = 0.01
        dmag_scat_lo = 2 * mag_scat    # extend fit for b.c.'s of deconvolute
        dmag_scat_hi = 1 * mag_scat
        self.mmin = 17.0
        #self.mmax = 23.3
        self.mmax=24.
        mags = np.arange(self.mmin - dmag_scat_lo, self.mmax + dmag_scat_hi, dmag, np.float32)
        numdens = np.zeros(mags.size)
        dndms = np.zeros(mags.size)
        for mi in xrange(len(mags)):
            numdens[mi] = np.abs(self.numden(mags[mi]))
            dndms[mi] = self.dndm(mags[mi])
        #print 'numden ', numdens[:10]
        #print mags[:10]
        # make no scatter splines
        self.log_numden_m_spl = interpolate.splrep(mags, log10(numdens))
        self.dndm_m_spl = interpolate.splrep(mags, dndms)
        self.m_log_numden_spl = interpolate.splrep(log10(numdens)[::-1], mags[::-1])
        # make scatter splines
        if self.scat:
            # deconvolve observed lf assuming scatter to find unscattered one
            dndms_scat = deconvolute(dndms, mag_scat, dmag, deconvol_iter_num)
            # chop off boundaries, unreliable
            #print mags.min(), mags.max() 
            #mags = mags[dmag_scat_lo / dmag:-dmag_scat_hi / dmag]
            #dndms_scat = dndms_scat[dmag_scat_lo / dmag:-dmag_scat_hi / dmag]
            #print mags.min(), mags.max() 
            # find spline to integrate over
            self.dndm_m_scat_spl = interpolate.splrep(mags, dndms_scat)
            numdens_scat = np.zeros(mags.size)
            for mi in xrange(mags.size):
                numdens_scat[mi] = np.abs(interpolate.splint(mags[mi], mags.max(), self.dndm_m_scat_spl))
                numdens_scat[mi] += 1e-9 * (1 - mi * 0.001)
            self.log_numden_m_scat_spl = interpolate.splrep(mags, log10(numdens_scat))
            self.m_log_numden_scat_spl = interpolate.splrep(log10(numdens_scat)[::-1], mags[::-1])



# Utility functions 
def deconvolute(y_conv, scatter, x_wid, iter_num=10):
    '''
    Get deconvolved version via Lucy routine.

    Import gaussian convoluted function, scatter, bin width, number of iterations.
    '''
    yit = y_conv
    for _ in xrange(iter_num):
        ratio = y_conv / ndimage.filters.gaussian_filter1d(yit, scatter / x_wid)
        yit = yit * ndimage.filters.gaussian_filter1d(ratio, scatter / x_wid)
        # this is part of lucy's routine, but seems less stable
        #yit = yit * ratio
    return yit


def elements(vals, lim=[-Inf, Inf], vis=None, vis_2=None, get_indices=False, dtype=np.int32):
    '''
    Get the indices of the input values that are within the input limit, that also are in input vis
    index array (if defined).
    Either of limits can have same range as vals.

    Import array, range to keep, prior indices of vals array to keep,
    other array to sub-sample in same way, whether to return selection indices of input vis array.
    '''
    if not isinstance(vals, np.ndarray):
        vals = np.array(vals)
    # check if input array
    if vis is None:
        vis = np.arange(vals.size, dtype=dtype)
    else:
        vals = vals[vis]
    vis_keep = vis
    # check if limit is just one value
    if np.isscalar(lim):
        keeps = (vals == lim)
    else:
        # sanity check - can delete this eventually
        if isinstance(lim[0], int) and isinstance(lim[1], int):
            if lim[0] == lim[1]:
                raise ValueError('input limit = %s, has same value' % lim)
            if lim[0] != lim[1] and 'int' in vals.dtype.name:
                print '! elements will not keep objects at lim[1] = %d' % lim[1]
        if not np.isscalar(lim[0]) or lim[0] > -Inf:
            keeps = (vals >= lim[0])
        else:
            keeps = None
        if not np.isscalar(lim[1]) or lim[1] < Inf:
            if keeps is None:
                keeps = (vals < lim[1])
            else:
                keeps *= (vals < lim[1])
        elif keeps is None:
            keeps = np.arange(vals.size, dtype=dtype)
    if get_indices:
        if vis_2 is not None:
            return vis_keep[keeps], vis_2[keeps], np.arange(vis.size, dtype=dtype)[keeps]
        else:
            return vis_keep[keeps], np.arange(vis.size, dtype=dtype)[keeps]
    else:
        if vis_2 is not None:
            return vis_keep[keeps], vis_2[keeps]
        else:
            return vis_keep[keeps]


def arange_length(array_or_length_or_imin=None, imax=None, dtype=np.int32):
    '''
    Get arange corresponding to input limits or input array size.

    Import array or array length or starting value (if latter, also need ending value).
    '''
    if imax is None:
        if np.isscalar(array_or_length_or_imin):
            num = array_or_length_or_imin
        else:
            num = len(array_or_length_or_imin)
        return np.arange(num, dtype=dtype)
    else:
        return np.arange(array_or_length_or_imin, imax, dtype=dtype)


def intersection_index(arr1, arr2):  
    """ 
    Find the indicies of the intersecting elements of arr1 and arr2.
    Takes approximately < 1 second
    """
    sort_arr1_indices = np.argsort(arr1)
    sort_arr2_indices = np.argsort(arr2)

    sorted_arr1 = arr1[sort_arr1_indices]
    sorted_arr2 = arr2[sort_arr2_indices]

    arr1_in1d = np.in1d(sorted_arr1, sorted_arr2)
    arr2_in1d = np.in1d(sorted_arr2, sorted_arr1)

    arr1_intersect_indices = sort_arr1_indices[arr1_in1d]
    arr2_intersect_indices = sort_arr2_indices[arr2_in1d]

    return arr1_intersect_indices, arr2_intersect_indices 



if __name__=='__main__': 
    #DownloadedCatalog(catalog='smdpl')
    sham_dict = { 
            'm_kind': 'mag_r', 
            'scat': 0.17, 
            'source': 'blanton', 
            'sham_prop': 'tailored'
            }
    shame = shamHalos(catalog='smdpl', sham_dict=sham_dict)
    shame.ReadHaloCatalog()
    shame.Write()

