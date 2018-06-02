import pandas as pd
import numpy as np
import os
from astropy.io import fits
import apogee_tools as ap

FULL_PATH  = os.path.realpath(__file__)
BASE = os.path.split(os.path.split(os.path.split(FULL_PATH)[0])[0])[0]
AP_PATH = os.environ['APOGEE_DATA'] 

import time

def getTelluric(**kwargs):

    """
    Get a telluric spectrum

    Parameters
    ----------
    wavelow:  int
              lower bound of the wavelength range

    wavehigh: int
              upper bound of the wavelength range

    airmass:  str
              airmass of the telluric model
    alpha:    float
              the power alpha parameter of the telluric model

    Returns
    -------
    telluric: arrays
              telluric wavelength and flux within specified region
    """

    alpha   = kwargs.get('alpha', 1)
    airmass = kwargs.get('airmass', ap.fix_param["airmass"])
    cut_rng = kwargs.get('cut_rng', [ap.data['orders'][0][0], ap.data['orders'][-1][-1]])

    am_key = {'1.0':'10', '1.5':'15'}

    if 'cut_rng' not in kwargs:
        print("Reading in telluric model between default wavelength {} and {} angstrom. \
            Use key word 'cut_rng'=[wl_min, wl_max] to specify otherwise.".format(cut_rng[0], cut_rng[1]))

    tfile     = 'pwv_R300k_airmass{}/LBL_A{}_s0_w005_R0300000_T.fits'.format(str(airmass), am_key[str(airmass)])
    tellurics = fits.open(BASE + '/libraries/TELLURIC/' + tfile)

    tell_wave = np.array(tellurics[1].data['lam'] * 10000)
    tell_flux = np.array(tellurics[1].data['trans'])**(alpha)

    #cut model wl grid wl array to bounds of telluric spectrum
    cut = np.where((tell_wave > cut_rng[0]) & (tell_wave < cut_rng[1]))[0]

    tell_sp = ap.Spectrum(wave=tell_wave[cut], flux=tell_flux[cut])

    return tell_sp


def applyTelluric(mdl, tell_sp, **kwargs):

    """
    Apply telluric model to PHOENIX model (or whatever grid you're using)
    @Elizabeth Moreno, Jessica Birky
    """

    alpha = kwargs.get('alpha', 1)
    method = kwargs.get('method', 'fast')
    
    mdl_wave  = mdl.wave
    mdl_flux  = mdl.flux
    tell_wave = tell_sp.wave
    tell_flux = tell_sp.flux

    min_wave = min(min(tell_wave), min(mdl_wave))
    max_wave = max(max(tell_wave), max(mdl_wave))
    cut1 = np.where((mdl_wave > min_wave) & (mdl_wave < max_wave))[0]
    cut2 = np.where((tell_wave > min_wave) & (tell_wave < max_wave))[0]

    mdl_wave = mdl_wave[cut1]
    mdl_flux = mdl_flux[cut1]
    tell_wave = tell_wave[cut2]
    tell_flux = tell_flux[cut2]

    if len(mdl_wave) > len(tell_wave): 

        cut3 = np.where((tell_wave > mdl_wave[0]) & (tell_wave < mdl_wave[-1]))[0]

        tell_wave = tell_wave[cut3]
        tell_flux = tell_flux[cut3]

        mdl_rs  = ap.integralResample(mdl_wave, mdl_flux, tell_wave, method=method)

        #Convolve model and telluric spectrum
        conv_flux = mdl_rs * tell_flux

        telluric_spec = ap.Spectrum(wave=tell_wave, flux=conv_flux)

    else:

        cut3 = np.where((mdl_wave > tell_wave[0]) & (mdl_wave < tell_wave[-1]))[0]

        mdl_wave = mdl_wave[cut3]
        mdl_flux = mdl_flux[cut3]

        tell_rs = ap.integralResample(tell_wave, tell_flux, mdl_wave, method=method)

        #Convolve model and telluric spectrum
        conv_flux = tell_rs * mdl_flux

        telluric_spec = ap.Spectrum(wave=mdl_wave, flux=conv_flux)

    
    return telluric_spec

