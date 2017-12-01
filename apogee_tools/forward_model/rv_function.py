import apogee_tools as ap

def rvShift(wave, **kwargs): 

    """
    Shift wavelenth of spectrum object by specific rv.
    @Jessica Birky

    Input: 'rv' : radial velocity (km/s)

    """
    
    rv = kwargs.get('rv', 0) 

    shift   = 1. + rv/299792.
    rv_wave = wave*shift

    return rv_wave


def rvShiftSpec(sp, **kwargs): 

    """
    Rv shift, given spectrum object, return spectrum object
    @Jessica Birky

    Input: 'rv' : radial velocity (km/s)

    """
    
    rv = kwargs.get('rv', 0) 

    shift   = 1. + rv/299792.
    rv_wave = sp.wave*shift
    rv_spec = ap.Spectrum(wave=rv_wave, flux=sp.flux)

    return rv_spec