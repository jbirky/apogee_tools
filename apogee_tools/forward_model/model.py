import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='serif')
import time
import os

import apogee_tools as ap
import apogee_tools.apogee_hack.spec.lsf as lsf


def initialize(**kwargs):

	t0 = time.time()

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

		# Evaluate lsf array
		xlsf = np.linspace(-7., 7., 43)
		lsf_array = lsf.eval(xlsf, fiber=[fiber])

		if ap.out['print_report'] == True:
			print('\n[{}s] MCMC initialization step complete.'.format(time.time() - t0))

		return init_param, step_param, init_theta, step_theta, fiber, tell_sp, lsf_array

	elif instrument == 'NIRSPEC':


		### DINO'S NIRSPEC READING CODE ###


		return init_param, step_param, init_theta, step_theta, fiber, tell_sp


def makeModel(**kwargs):

	"""
	Returns the synthesized model of a given set of parameters, and plots (optional).

	Input:  'params' : dictionary of parameters specified like {'teff': 3051, 'logg': 5.2, 'fe_h': -0.25, 'vsini': 10., 'rv': -12, 'alpha': 0.8}
			'fiber'  : APOGEE fiber for particular ap1D visit spectrum you are fitting model to

	Output: 'synth_model' : spectrum object synthesized to input paramters
	"""

	params = kwargs.get('params')
	plot   = kwargs.get('plot', False)
	res    = kwargs.get('res', '300k')
	grid   = kwargs.get('grid', 'phoenix').lower()
	bands  = kwargs.get('bands', [[15200,15800],[15860,16425],[16475,16935]])
	method = kwargs.get('method', 'fast')

	# mdl_name = r'Teff = {}, logg = {}, Fe/H = {}, vsini = {}, rv = {}, $\alpha$ = {}'.format(params['teff'], params['logg'], params['fe_h'], params['vsini'], params['rv'], params['alpha'])
	labels = [params['teff'], params['logg'], params['fe_h']]

	out_str = " ".join(["{}={}".format(key, params[key]) for key in params.keys()])
	if ap.out['print_report'] == True:
		print('\n##################################################')
		print('Making model:', out_str)


	# =============================================================
	# 1. Interpolate model grids at give teff, logg, fe/h
	# =============================================================
	
	t0 = time.time()

	if ap.fix_param["interp_method"] == 'cannon':
		interp_sp = ap.interpolateGrid(labels=labels, res=res, grid=grid)
	else:
		interp_sp = ap.interpolateGrid(labels=labels)
	interp_sp.flux = interp_sp.flux/max(interp_sp.flux)

	t1 = time.time()

	if ap.out['print_report'] == True:
		print("[{}s] Interpolated model".format(str(t1-t0)))

	# =============================================================
	# 2. Apply radial velocity
	# =============================================================

	rv_sp = ap.rvShiftSpec(interp_sp, rv=params['rv'])

	t2 = time.time()

	if ap.out['print_report'] == True:
		print("[{}s] Shifted radial velocity".format(str(t2-t1)))

	# =============================================================
	# 3. Apply rotational velocity broadening
	# =============================================================

	rot_sp = ap.applyVsini(rv_sp, vsini=params['vsini'])

	t3 = time.time()

	if ap.out['print_report'] == True:
		print("[{}s] Applied vsini broadening".format(str(t3-t2)))

	# =============================================================
	# 4. Apply telluric spectrum
	# =============================================================

	telluric_model = kwargs.get('telluric', ap.getTelluric(cut_rng=[ap.data['orders'][0][0], ap.data['orders'][-1][-1]]))
	tell_sp = ap.applyTelluric(rot_sp, telluric_model, alpha=params['alpha'], method=method)

	t4 = time.time()

	if ap.out['print_report'] == True:
		print("[{}s] Convolved telluric model".format(str(t4-t3)))

	# =============================================================
	# 5. Apply APOGEE LSF function
	# =============================================================

	# Specify either lsf array or fiber number. 
	# Note for multiple iterations of the function specifying the lsf array will be faster.

	if 'lsf' in kwargs:
		lsf_array  = kwargs.get('lsf')
		lsf_sp     = ap.convolveLsf(tell_sp, lsf=lsf_array)

	elif ('lsf' not in kwargs) and ('fiber' in kwargs):
		import apogee_tools.apogee_hack.spec.lsf as lsf

		xlsf 	  = np.linspace(-7., 7., 43)
		fiber     = kwargs.get('fiber')
		lsf_array = lsf.eval(xlsf, fiber=fiber)
		lsf_sp    = ap.convolveLsf(tell_sp, lsf=lsf_array)

	else:
		print("Error convolving LSF: ap.convolveLsf(tell_sp) requires key word input 'lsf' or 'fiber'.")

	t5 = time.time()

	if ap.out['print_report'] == True:
		print("[{}s] Applied LSF broadening \n".format(str(t5-t4)))

	# =============================================================

	cut = np.where((lsf_sp.wave > ap.data["orders"][0][0]) & (lsf_sp.wave < ap.data["orders"][-1][1]))[0]
	synth_model = ap.Spectrum(wave=lsf_sp.wave[cut], flux=lsf_sp.flux[cut])

	if plot == True:
		xrange = kwargs.get('xrange', [interp_sp.wave[0], interp_sp.wave[-1]])
		yrange = kwargs.get('yrange', [-.5, 1.1])

		plt.figure(1, figsize=(16,6))  
		plt.plot(interp_sp.wave, interp_sp.flux, alpha=.7, linewidth=1, label=r'$\lambda_{grid}$: %s'%(ap.model["grid_name"]))
		plt.plot(rv_sp.wave, rv_sp.flux-.15, label=r'$\lambda_{grid}$ * (1+$\frac{rv}{c}$)', alpha=.7, linewidth=1)
		plt.plot(rot_sp.wave, rot_sp.flux-.3, label=r'$\lambda_{grid}$ * (1+$\frac{rv}{c}$) * vsini', alpha=.7, linewidth=1)
		plt.plot(tell_sp.wave, tell_sp.flux-.45, label=r'$\lambda_{grid}$ * (1+$\frac{rv}{c}$) * vsini * telluric', alpha=.7, linewidth=1)
		plt.plot(lsf_sp.wave, lsf_sp.flux-.8, label=r'$\lambda_{grid}$ * (1+$\frac{rv}{c}$) * vsini * telluric * lsf', alpha=.7, linewidth=1)
		txt = r"Teff=%s, logg=%s, [Fe/H]=%s, rv=%s, vsini=%s, $\alpha$=%s"%(params['teff'], params['logg'], params['fe_h'], params['rv'], params['vsini'], params['alpha'])
		plt.text(xrange[-1]-20, yrange[0]+.03, txt, ha='right', fontsize=14)
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
	lsf    = kwargs.get('lsf')
	plot   = kwargs.get('plot', False)
	telluric_model = kwargs.get('telluric', ap.getTelluric(cut_rng=[ap.data['orders'][0][0], ap.data['orders'][-1][-1]]))

	# normalize and apply sigma clip to data flux
	data.mask(sigma=ap.data["sigma_clip"], pixel_buffer=ap.data["pixel_buffer"])
	data.flux = data.flux/np.nanmax(data.flux)

	# synthesize model
	synth_mdl = ap.makeModel(params=theta, lsf=lsf, telluric=telluric_model)

	# multiply model by continuum polynomial
	cont_sp = ap.continuum(data, synth_mdl, bands=ap.data["orders"], deg=ap.fix_param["cont_deg"])

	print('data median flux', np.median(data.flux), np.nanmedian(data.flux))
	print('model median flux', np.median(cont_sp.flux), np.nanmedian(cont_sp.flux))
	# return chi-squared fit
	chisq = ap.compareSpectra(data, cont_sp, fit_scale=False)[0]

	if plot == True:
		plt.figure(figsize=[12,5])
		plt.plot(data.wave, data.flux, label=str(data.name), color='k', alpha=.7, linewidth=1)
		plt.plot(cont_sp.wave, cont_sp.flux, label=r'$\chi^2=%s$'%(str(chisq)), color='r', alpha=.7, linewidth=1)
		plt.legend(loc='upper right')
		plt.ylabel(r'$F_{\lambda}$ + offset', fontsize=15)
		plt.xlabel(r'$\lambda$', fontsize=15)
		plt.show()
		plt.close()

	return chisq

