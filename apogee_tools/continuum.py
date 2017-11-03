#!/usr/bin/env python
#
# This is a function for continuum correction
#
# Oct. 27 2017
# @Dino Hsu

import numpy
import apogee_tools as ap
import matplotlib.pyplot as plt
import copy

def continuum(data, mdl, contfitdeg=5):
	"""
	This function returns a continuum corrected model.
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


# Below is a test plot
def test():
	ap.download('2M03290406+3117075', type='aspcap')
	mdl      = ap.getModel(params=[3200, 5.0, 0.0], grid='BTSETTLb', xrange=[15200,16940])
	data = ap.Spectrum(id='2M03290406+3117075', type='aspcap')
	mdlcont = continuum(data, mdl)
	data_int = numpy.interp(mdl.wave, data.wave, data.flux)

	plt.figure(figsize=(12,4))
	plt.plot(mdl.wave, mdl.flux, color='k', alpha=.8, label='the original model')
	plt.plot(mdl.wave, mdlcont.flux, color='c', alpha=.8, label='the continuum corrected model')
	plt.plot(mdl.wave, data_int, color='r', alpha=.8, label='the data used in polynomial fitting')
	plt.xlabel('Wavelength (Angstrom)')
	plt.ylabel('Flux')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
	plt.show()
	plt.close()






