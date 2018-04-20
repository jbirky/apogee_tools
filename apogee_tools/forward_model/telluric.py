import pandas as pd
import numpy as np
import os
from astropy.io import fits
import apogee_tools as ap

FULL_PATH  = os.path.realpath(__file__)
BASE = os.path.split(os.path.split(os.path.split(FULL_PATH)[0])[0])[0]
AP_PATH = os.environ['APOGEE_DATA'] 



def getTelluric(wavelow=21600, wavehigh=21900, **kwargs):

    """
    Get a telluric spectrum.

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
    airmass = kwargs.get('airmass', '1.0')

    am_key = {'1.0':'10', '1.5':'15'}

    tfile     = 'pwv_R300k_airmass{}/LBL_A{}_s0_w005_R0300000_T.fits'.format(airmass, am_key[airmass])
    tellurics = fits.open(BASE + '/libraries/TELLURIC/' + tfile)

    tell_wave = np.array(tellurics[1].data['lam'] * 10000)
    tell_flux = np.array(tellurics[1].data['trans'])**(alpha)

    #cut model wl grid wl array to bounds of telluric spectrum
    cut = np.where( (tell_wave > wavelow) & (tell_wave < wavehigh) )

    return tell_wave[cut], tell_flux[cut]


def applyTelluric(mdl, **kwargs):

    """
    Apply telluric model to PHOENIX model (or whatever grid you're using)
    @Elizabeth Moreno, Jessica Birky
    """

    alpha = kwargs.get('alpha', 1)
    airmass = kwargs.get('airmass', '1.0')

    am_key = {'1.0':'10', '1.5':'15'}

    mdl_wave  = mdl.wave
    mdl_flux  = mdl.flux

    # tellurics = pd.read_csv(BASE + '/libraries/lw_solartrans_apogee.csv')
    # tell_wave = np.asarray(tellurics['wave'] * 10000.0)
    # tell_flux = np.asarray(tellurics['trans'])**(alpha)

    tfile = 'pwv_R300k_airmass{}/LBL_A{}_s0_w005_R0300000_T.fits'.format(airmass, am_key[airmass])
    tellurics = fits.open(BASE + '/libraries/TELLURIC/' + tfile)

    tell_wave = np.array(tellurics[1].data['lam'] * 10000)
    tell_flux = np.array(tellurics[1].data['trans'])**(alpha)

    #Resample higher res spectrum to lower res spectrum
    try: #if res tell > res mdl, resample tell to mdl

        #cut telluric spectrum wl array to bounds of model wl grid
        cut = np.where((mdl_wave > tell_wave[0]) & (mdl_wave < tell_wave[-1]))[0][5:-5]
        mdl_wave = mdl_wave[cut]
        mdl_flux = mdl_flux[cut]

        tell_rs = ap.integralResample(tell_wave, tell_flux, mdl_wave)
        print('Downsampled telluric spectrum resolution to model resolution.')

        #Convolve model and telluric spectrum
        conv_flux = tell_rs * mdl_flux

        telluric_spec = ap.Spectrum(wave=mdl_wave, flux=conv_flux)

    except: #if res mdl > res tell, resample mdl to tell 

        #cut model wl grid wl array to bounds of telluric spectrum
        cut = np.where((tell_wave > mdl_wave[0]) & (tell_wave < mdl_wave[-1]))[0][5:-5]
        tell_wave = tell_wave[cut]
        tell_flux = tell_flux[cut]

        mdl_rs  = ap.integralResample(mdl_wave, mdl_flux, tell_wave)
        print('Downsampled model resolution to telluric spectrum resolution.')

        #Convolve model and telluric spectrum
        conv_flux = mdl_rs * tell_flux

        telluric_spec = ap.Spectrum(wave=tell_wave, flux=conv_flux)
    
    return telluric_spec

