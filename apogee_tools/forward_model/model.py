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
	init_param = {key:ap.init[key] for key in ap.init.keys()}
	step_param = {key:ap.step[key] for key in ap.init.keys()}

	init_theta = {key:ap.init[key] for key in ap.model["theta"]}
	step_theta = {key:ap.step[key] for key in ap.model["theta"]}

	tell_sp = ap.getTelluric(airmass=str(ap.fix_param["airmass"]), cut_rng=[ap.data['orders'][0][0], ap.data['orders'][-1][-1]])

	if instrument == 'APOGEE':

		# Read in data that model will be generated for
		data = ap.Apogee(id=ap.data['ID'], type='apvisit', visit=ap.data['visit'])

		# Get APOGEE lsf fiber number
		fiber = ap.searchVisits(id_name=ap.data['ID'])[3][ap.data['visit']-1]

		return init_param, step_param, init_theta, step_theta, fiber

	elif instrument == 'NIRSPEC':


		### DINO'S NIRSPEC READING CODE ###


		return init_par, step_par, init_theta, step_theta, fiber, tell_sp


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
	bands  = kwargs.get('bands', [[15200,15800],[15860,16425],[16475,16935]])

	# mdl_name = r'Teff = {}, logg = {}, Fe/H = {}, vsini = {}, rv = {}, $\alpha$ = {}'.format(params['teff'], params['logg'], params['fe_h'], params['vsini'], params['rv'], params['alpha'])
	labels = [params['teff'], params['logg'], params['fe_h']]

	#Interpolate model grids at give teff, logg, fe/h
	interp_sp = ap.interpolateGrid(labels=labels, res=res, grid=grid)
	interp_sp.flux = interp_sp.flux/max(interp_sp.flux)

	#Apply radial velocity
	rv_sp   = ap.rvShiftSpec(interp_sp, rv=params['rv'])

	#Apply rotational velocity broadening
	rot_sp  = ap.applyVsini(rv_sp, vsini=params['vsini'])

	#Apply telluric spectrum
	tell_sp = ap.applyTelluric(rot_sp, alpha=params['alpha'])

	#Apply APOGEE LSF function
	lsf_sp  = ap.convolveLsf(tell_sp, fiber=fiber)

	# synth_model = ap.Spectrum(wave=lsf_sp.wave, flux=lsf_sp.flux, name=mdl_name)
	synth_model = ap.Spectrum(wave=lsf_sp.wave, flux=lsf_sp.flux)

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


def returnModelFit(data, theta, **kwargs):

	"""
	Function to be fed into the MCMC.

	Input:  'data'   : spectrum obj of data
			'params' : parameters to synthesize a model with 'makeModel'

	Output: 'chi' : chi-squared fit between data and synthesized model
	"""

	params = kwargs.get('params')
	fiber  = kwargs.get('fiber', 40)

	# normalize and apply sigma clip to data flux
	data.mask(sigma=ap.data["sigma_clip"], pixel_buffer=[0,0])
	data.flux = data.flux/max(data.flux[np.isfinite(data.flux)])

	# synthesize model
	synth_mdl = ap.makeModel(params=theta, fiber=fiber)

	# multiply model by continuum polynomial
	cont_sp = ap.continuum(data, synth_mdl, bands=ap.data["orders"], deg=ap.fix_param["cont_deg"], norm=True)

	# plt.figure(figsize=[12,5])
	# plt.plot(data.wave, data.flux)
	# plt.plot(cont_sp.wave, cont_sp.flux)
	# plt.show()
	# plt.close()

	# return chi-squared fit
	chisq = ap.compareSpectra(data, cont_sp, fit_scale=False)[0]

	return chisq


def fitMCMC(**kwargs):

	return None

