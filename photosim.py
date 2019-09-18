import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from glob import glob
from astropy.cosmology import Planck15 as cosmo
from os import system
from os.path import isfile
from os.path import isdir
from os import makedirs
from tqdm import tqdm
from pandas import read_csv

# dirlist = ['./cache/neb_em', './cache/neb_noem', './cache/noneb_noem']
# if not all(list(map(isdir, dirlist))):
#     for directory in dirlist:
#         makedirs(directory)


with open('/home/adam/Research/fsps/src/sps_vars.f90', 'r') as readfile:
    lines = readfile.readlines()
    if ('MILES 0' in lines[8]) or ('MIST 0' in lines[14]):
        system('fsps_mm')

import fsps


font = {'family':'Roboto', 'weight':'light'}

Lsun = 3.848e33 # erg/s
c_angstroms = 3e8 * 10**10 # in AA/s


class photosim:

    def __init__(self, inputfolder = './adam_synth_spectra/'):

        self.starpop = None


    def redshift(self, redshift, wavelength, l_nu):

        if redshift == 0:
            return wavelength, l_nu
        else:

            new_wavelength = wavelength * (1.+redshift)
            flux = l_nu * (1+redshift) / (4*np.pi*np.power(cosmo.luminosity_distance(redshift).to('cm').value,2))

            return new_wavelength, flux


    def find_spectrum(self, tage, metallicity = 1., imf_type = 1, sfh_type = 3, sfh_params = {'sfr':np.array([1.,1.]), 't':np.array([0,13])}, dust_type = 2, Av = 0, emline = True, nebcont = True, increase_ssp = True, delay_csf = False, peraa = False):

        # Cache and retreive spectra

    
        # fname_params = (Av, sfh_params['tau'])

        # fname = '%.2f_%.5f.spec' % fname_params

        # fstub = './cache/'

        # if emline and nebcont:
        #     fstub = fstub + 'neb_em/'
        # elif nebcont:
        #     fstub = fstub + 'neb_noem/'
        # elif not emline and not nebcont:
        #     fstub = fstub + 'noneb_noem/'

        # if isfile(fstub + fname):

        #     waves, lum = np.loadtxt(fstub + fname, unpack = True)

        # else:


        if sfh_type == 6:

            time = np.linspace(0, 10, 250)
            sfr = np.exp(time/sfh_params['tau'])
            if not 't' in sfh_params.keys():
                sfh_params['t'] = time
            if not 'sfr' in sfh_params.keys():
                sfh_params['sfr'] = sfr

            return self.find_spectrum(tage, metallicity = metallicity, imf_type = imf_type, sfh_type = 3, sfh_params = sfh_params, dust_type = dust_type, emline = emline, nebcont = nebcont, peraa = peraa, Av = Av)

            
        if self.starpop == None:
            self.starpop = fsps.StellarPopulation(zcontinuous=1, add_neb_emission = True, nebemlineinspec = emline, add_neb_continuum = nebcont,
                imf_type = imf_type, dust_type = dust_type, sfh = sfh_type, logzsol = np.log10(metallicity), tage = tage, gas_logz = np.log10(metallicity), dust2 = Av/1.08574)
        else:
            self.starpop.params['imf_type'] = imf_type
            self.starpop.params['logzsol'] = np.log10(metallicity)
            self.starpop.params['gas_logz'] = np.log10(metallicity)
            self.starpop.params['dust_type'] = dust_type
            self.starpop.params['sfh'] = sfh_type
            self.starpop.params['nebemlineinspec'] = emline
            self.starpop.params['add_neb_continuum'] = nebcont
            self.starpop.params['dust2'] = Av/1.08574

        if sfh_type == 3:
            # Form stars at 1Msun/yr for 1Gyr, then spike to 200Msun/yr
            # starpop.set_tabular_sfh(np.array([0,0.999,1,1.1]), np.array([1,1,10,10]))
            self.starpop.set_tabular_sfh(sfh_params['t'], sfh_params['sfr'])
            if delay_csf:
                tage = tage + 1.

        elif sfh_type == 1:

            if not 'tau' in sfh_params.keys():
                sfh_params['tau'] = 1. 

            #Set everything to zeros if they don't already exist

            for thiskey in ['const', 'sf_start', 'sf_trunc', 'tburst', 'fburst']:
                if thiskey not in sfh_params.keys():
                    sfh_params[thiskey] = 0.

            for thiskey in list(sfh_params.keys()):
                self.starpop.params[thiskey] = sfh_params[thiskey]
       
        waves, lum = self.starpop.get_spectrum(tage = tage, peraa = peraa)

        # peraa doesn't work for python FSPS
        # if peraa:
        #     lum = lum * c_angstroms / (waves**2)

        lum = lum * Lsun

        # np.savetxt(fstub + fname, np.vstack((waves, lum)).T, fmt = '%e   %.10e')

        if increase_ssp and sfh_type == 0:
            lum = lum * 10.**7

        return waves, lum

        
