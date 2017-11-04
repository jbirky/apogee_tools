import numpy
import apogee_tools as ap
import matplotlib.pyplot as plt
import copy

def continuum(data, mdl, contfitdeg=5):

	"""
	This function returns a continuum corrected model.
	@Dino Chih-Chun Hsu

	data: the data used in the fitting as a polynomial
	mdl: the model being corrected. 
	contfitdeg: the degree of the fitting polynomial. The default vaule is 5.
	"""
	
	if contfitdeg is None:
		contfitdeg = 5

	# interpolate the shape of the data to be the same as that of the model
	data_int = numpy.interp(mdl.wave, data.wave, data.flux)

	mdlcont = copy.deepcopy(mdl)
	mdldiv = data_int/mdl.flux
	pcont = numpy.polyfit(mdl.wave, mdldiv, contfitdeg)
	mdlcont.flux *= numpy.polyval(pcont, mdlcont.wave)

	return mdlcont
