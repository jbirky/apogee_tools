import numpy as np
import matplotlib.pyplot as plt
import copy


def continuum(data, mdl, **kwargs):

	"""
	This function returns a continuum corrected model.
	@Dino Chih-Chun Hsu

	data: the data used in the fitting as a polynomial
	mdl: the model being corrected. 
	contfitdeg: the degree of the fitting polynomial. The default vaule is 5.
	"""

	deg = kwargs.get('deg', 1) 

	flux_in_rng = np.where((mdl.wave > data.wave[0]) & (mdl.wave < data.wave[-1]))[0]
	mdl_wave = mdl.wave[flux_in_rng]
	mdl_flux = mdl.flux[flux_in_rng]
	   
	mdl_flux_vals = [x for x in mdl_flux if str(x) != 'nan']
	mdl_flux = mdl_flux/max(mdl_flux_vals)
	mdl_flux = [1 if str(x) == 'nan' else x for x in mdl_flux]

	mdl_int = np.interp(data.wave, mdl_wave, mdl_flux)
	mdlcont = ap.Spectrum(wave=mdl_wave, flux=mdl_flux)
	mdldiv  = data.flux/mdl_int

	mdldiv_vals = [x for x in mdldiv if str(x) != 'nan']
	mdldiv = [np.mean(mdldiv_vals) if str(x) == 'nan' else x for x in mdldiv]
	pcont  = np.polyfit(data.wave, mdldiv, deg)
	    
	mdlcont.flux *= np.polyval(pcont, mdlcont.wave)

	return mdlcont

