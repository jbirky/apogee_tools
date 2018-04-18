import numpy as np
import matplotlib.pyplot as plt
import copy
import apogee_tools as ap


def continuum(data, mdl, **kwargs):
    
    """
    This function returns a continuum multiplied model.
    @Jessica Birky, Dino Chih-Chun Hsu

    Input:  'data'  : the data used in the fitting as a polynomial
            'mdl'   : the model being corrected
            'deg'   : the degree of the fitting polynomial. The default vaule is 5.
            'bands' : list of wavelenth ranges in over which separate polynomials will be fit;
                      any wavelengths outside of these ranges will be replaced with nan in the returned model
    
    Output: 'mdl_poly' : spectrum object of model times polynomial
    """
    
    deg = kwargs.get('deg', 5) 
    bands = kwargs.get('bands', [[data.wave[0], data.wave[-1]]])

    # Cut data to region of bands
    dcut_rng = np.where((data.wave >= bands[0][0]) & (data.wave <= bands[-1][1]))
    data_wave = data.wave[dcut_rng]
    data_flux = data.flux[dcut_rng]
    
    # Cut model to wavelength region of the data
    data_rng = np.where((mdl.wave >= data_wave[0]) & (mdl.wave <= data_wave[-1]))[0]
    mdl_wave = mdl.wave[data_rng]
    mdl_flux = mdl.flux[data_rng]

    # Model flux resampled down to data resolution
    mdl_int = np.interp(data_wave, mdl_wave, mdl_flux)
    mdl_res = ap.Spectrum(wave=data_wave, flux=mdl_int)

    # Data divided by model
    mdldiv  = data_flux/mdl_res.flux
    
    # Fit a polynomial to each band separately
    band_indices = []
    band_flux = []
    
    for bnd in bands:
        band_rng = np.where((mdl_res.wave >= bnd[0]) & (mdl_res.wave <= bnd[-1]))[0]
        band_indices.append(band_rng)
        
        # Temporarily remove nans from mdldiv before fitting continuum polynomial
        rm_nan = np.isfinite(mdldiv[band_rng])
        pcont = np.polyfit(mdl_res.wave[band_rng][rm_nan], mdldiv[band_rng][rm_nan], deg)
        
        poly_mult = mdl_res.flux[band_rng] * np.polyval(pcont, mdl_res.wave[band_rng])
        band_flux.append(poly_mult)
        
    # Flux at the gaps between the bands should be nans
    gap_flux = []
    wave_gaps = [[bands[i][1], bands[i+1][0]] for i in range(len(bands)-1)]
    
    for gap in wave_gaps:
        gap_rng = np.where((mdl_res.wave >= gap[0]) & (mdl_res.wave <= gap[-1]))[0]      
        gap_flux.append([np.nan for g in gap_rng])
    
    # Make sure that output spectrum object has continuous wavelength array
    poly_flux = [band_flux[0]]
    
    for i in range(len(wave_gaps)):
        poly_flux.append(gap_flux[i])
        poly_flux.append(band_flux[i+1])
    
    # Return model flux times polynomial
    mdl_poly = ap.Spectrum(wave=mdl_res.wave, flux=np.concatenate(poly_flux))

    return mdl_poly