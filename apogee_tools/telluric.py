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

    mdl_wave  = mdl.wave
    mdl_flux  = mdl.flux

    tell_wave = np.asarray(tellurics['wave'] * 10000.0)
    tell_flux = np.asarray(tellurics['trans'])**(alpha)

    #Resample higher res spectrum to lower res spectrum
    try: 
        #if res tell > res mdl, resample tell to mdl
        tell_rs = splat.integralResample(tell_wave, tell_flux, mdl_wave)

        #Convolve model and telluric spectrum
        conv_flux = tell_rs * mdl_flux

    except: 
        #if res mdl > res tell, resample mdl to tell 
        mdl_rs  = splat.integralResample(mdl_wave, mfl_flux, tell_wave)

        #Convolve model and telluric spectrum
        conv_flux = mdl_rs * tell_flux
        
    telluric_spec = ap.Spectrum(wave=mdl_wave, flux=conv_flux)
    
    return telluric_spec