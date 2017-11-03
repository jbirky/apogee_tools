import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='serif')
import apogee_tools as ap
import os


def makeModel(**kwargs):

	"""
	Input:  'params' : dictionary of parameters specified like {'teff': 3051, 'logg': 5.2, 'z': -0.25, 'vsini': 10., 'rv': -12, 'alpha': 0.8}

	Output: 
	"""

	params = kwargs.get('params')
	plot   = kwargs.get('plot', False)

	labels = [params['teff'], params['logg'], params['z']]

	#Interpolate model grids at give teff, logg, fe/h
	interp_sp = ap.interpolateGrid(labels=labels)

	#Apply radial velocity
	rv_sp   = ap.spec_tools.rvShiftSpec(interp_sp, rv=params['rv'])

	#Apply rotational velocity broadening
	rot_sp  = ap.applyVsini(rv_sp, vsini=params['vsini'])

	#Apply telluric spectrum
	tell_sp = ap.applyTelluric(rot_sp)

	if plot == True:
		plt.figure(1, figsize=(10,6))  
		plt.plot(interp_sp.wave, interp_sp.flux, label=r'Teff = %s, logg = %s, Fe/H = %s'%(params['teff'], params['logg'], params['z']))
		plt.plot(rv_sp.wave, rv_sp.flux, label=r'RV (%s km/s)'%(params['rv']))
		plt.plot(rot_sp.wave, rot_sp.flux, label=r'RV + rot (%s km/s)'%(params['vsini']))
		plt.plot(tell_sp.wave, tell_sp.flux, label=r'RV + rot + telluric')
		plt.xlim(interp_sp.wave[0], interp_sp.wave[-1])
		plt.ylim(0.7, 1.1)
		plt.legend(frameon=False)
		plt.show()

	return tell_sp


def returnModelFit(data, synth_mdl, **kwargs):

	"""
	Function to be fed into the MCMC.
	
	Input:  'data'   : spectrum obj of data
			'params' : parameters to synthesize a model with 'makeModel'

	Output: 'chi' : chi-squared fit between data and synthesized model
	"""

	params = kwargs.get('params')

	synth_mdl = makeModel(params=params)
	chi = ap.compareSpectra(data, synth_mdl)

	return chi

