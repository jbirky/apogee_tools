import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='serif')
import apogee_tools as ap
import os


def initialize(**kwargs):

	# Read in configuration parameters

	instrument = ap.data["instrument"]

	# Read initilization and step parameters
	init_par = {key:ap.init[key] for key in ap.grid['parname']}
	step_par = {key:ap.step[key] for key in ap.grid['parname']}

	if instrument == 'APOGEE':

		# Read in data that model will be generated for
		data = ap.Spectrum(id=ap.data['ID'], type='apvisit', visit=ap.data['visit'])

		# Get APOGEE lsf fiber number
		fiber = ap.searchVisits(id_name=ap.data['ID'])[3][ap.data['visit']-1]

		return init_par, step_par, fiber

	elif instrument == 'NIRSPEC':

		return init_par, step_par


def makeModel(**kwargs):

	"""
	Returns the synthesized model of a given set of parameters, and plots (optional).

	Input:  'params' : dictionary of parameters specified like {'teff': 3051, 'logg': 5.2, 'fe_h': -0.25, 'vsini': 10., 'rv': -12, 'alpha': 0.8}
			'fiber'  : APOGEE fiber for particular ap1D visit spectrum you are fitting model to

	Output: 'synth_model' : spectrum object synthesized to input paramters
	"""

	params = kwargs.get('params')
	fiber  = kwargs.get('fiber', 40)
	plot   = kwargs.get('plot', False)
	res    = kwargs.get('res', '300k')
	grid   = kwargs.get('grid', 'phoenix').lower()

	mdl_name = r'Teff = {}, logg = {}, Fe/H = {}, vsini = {}, rv = {}, $\alpha$ = {}'.format(params['teff'], params['logg'], params['z'], params['vsini'], params['rv'], params['alpha'])
	labels = [params['teff'], params['logg'], params['z']]

	#Interpolate model grids at give teff, logg, fe/h
	interp_sp = ap.interpolateGrid(labels=labels, res=res, grid=grid)
	interp_sp.flux = interp_sp.flux/max(interp_sp.flux)

	#Apply radial velocity
	rv_sp   = ap.spec_tools.rvShiftSpec(interp_sp, rv=params['rv'])

	#Apply rotational velocity broadening
	rot_sp  = ap.applyVsini(rv_sp, vsini=params['vsini'])

	#Apply telluric spectrum
	tell_sp = ap.applyTelluric(rot_sp, alpha=params['alpha'])

	#Apply APOGEE LSF function
	lsf_sp  = ap.convolveLsf(tell_sp, fiber=fiber)

	synth_model = ap.Spectrum(wave=lsf_sp.wave, flux=lsf_sp.flux, name=mdl_name)

	if plot == True:
		xrange = kwargs.get('xrange', [interp_sp.wave[0], interp_sp.wave[-1]])
		yrange = kwargs.get('yrange', [-.4, 1.1])

		plt.figure(1, figsize=(16,6))  
		plt.plot(interp_sp.wave, interp_sp.flux, alpha=.7, linewidth=1, label=r'Teff = %s, logg = %s, Fe/H = %s'%(params['teff'], params['logg'], params['z']))
		plt.plot(rv_sp.wave, rv_sp.flux-.15, label=r'RV (%s km/s)'%(params['rv']), alpha=.7, linewidth=1)
		plt.plot(rot_sp.wave, rot_sp.flux-.3, label=r'RV + rot (%s km/s)'%(params['vsini']), alpha=.7, linewidth=1)
		plt.plot(tell_sp.wave, tell_sp.flux-.45, label=r'RV + rot + telluric ($\alpha$ = %s)'%(params['alpha']), alpha=.7, linewidth=1)
		plt.plot(lsf_sp.wave, lsf_sp.flux-.8, label=r'RV + rot + telluric + lsf', alpha=.7, linewidth=1)
		plt.xlim(xrange)
		plt.ylim(yrange)
		plt.legend(loc='lower left', frameon=False, fontsize=12)
		plt.ylabel(r'$F_{\lambda}$ + offset', fontsize=15)
		plt.xlabel(r'$\lambda$', fontsize=15)
		plt.show()

	return synth_model


def returnModelFit(data, synth_mdl, **kwargs):

	"""
	Function to be fed into the MCMC.

	Input:  'data'   : spectrum obj of data
			'params' : parameters to synthesize a model with 'makeModel'

	Output: 'chi' : chi-squared fit between data and synthesized model
	"""

	params = kwargs.get('params')
	fiber  = kwargs.get('fiber', 40)

	synth_mdl = makeModel(params=params, fiber=fiber)
	chi = ap.compareSpectra(data, synth_mdl)

	return chi

