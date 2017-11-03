import pandas as pd
import numpy as np
import os
import apogee_tools as ap
import splat

FULL_PATH  = os.path.realpath(__file__)
BASE = os.path.split(os.path.split(FULL_PATH)[0])[0]
AP_PATH = os.environ['APOGEE_DATA'] 


def applyTelluric(wave, flux, **kwargs):

    """
    Apply telluric model to PHOENIX model (or whatever grid you're using)
    @ Elizabeth Moreno
    """

    tellurics = pd.read_csv(BASE + '/libraries/lw_solartrans_apogee.csv')

    WaveHigh  = np.asarray(tellurics['wave'] * 10000.0)
    TransHigh = np.asarray(tellurics['trans'])

    #Resampling
    TransLow  = splat.integralResample(WaveHigh, TransHigh, wave)

    #Getting the flux with the transmission 
    FluxWithTrans = TransLow * flux

    telluric_flux = FluxWithTrans
    telluric_spec = ap.Spectrum(wave=mdl.wave, flux=telluric_flux)
    
    return telluric_spec