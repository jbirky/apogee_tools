import pandas as pd
import numpy as np
import os
import apogee_tools as ap
import splat

FULL_PATH  = os.path.realpath(__file__)
BASE = os.path.split(os.path.split(FULL_PATH)[0])[0]
AP_PATH = os.environ['APOGEE_DATA'] 


def applyTelluric(mdl, **kwargs):

    """
    Apply telluric model to PHOENIX model (or whatever grid you're using)
    @ Elizabeth Moreno
    """

    alpha = kwargs.get('alpha', 1)

    tellurics = pd.read_csv(BASE + '/libraries/lw_solartrans_apogee.csv')

    WaveLow  = mdl.wave
    FluxLow  = mdl.flux

    WaveHigh = np.asarray(tellurics['wave'] * 10000.0)
    TransHigh = np.asarray(tellurics['trans'])**(alpha)

    #Resampling
    TransLow = splat.integralResample(WaveHigh, TransHigh, WaveLow)

    #Getting the flux with the transmission 
    FluxWithTrans = TransLow * FluxLow

    telluric_flux = FluxWithTrans
    telluric_spec = ap.Spectrum(wave=mdl.wave, flux=telluric_flux)
    # telluric_spec.flux = telluric_spec.flux**(alpha)
    
    return telluric_spec