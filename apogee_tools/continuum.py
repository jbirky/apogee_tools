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
	
	deg = kwargs.get('deg', 5)

	# interpolate the shape of the data to be the same as that of the model
	data_int = np.interp(mdl.wave, data.wave, data.flux)

	mdlcont = copy.deepcopy(mdl)
	mdldiv  = data_int/mdl.flux
	pcont   = np.polyfit(mdl.wave, mdldiv, deg)
	mdlcont.flux *= np.polyval(pcont, mdlcont.wave)

	return mdlcont
