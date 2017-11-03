import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='serif')
import apogee_tools as ap
import os


def prepData(apogee_id, **kwargs):

	"""
	Download and load ap1D spectrum, get lsf fiber number

	Input:  'apogee_id' : ex '2M00452143+1634446'
			'visit'     : visit number that will be read in

	Output: 'data'  : spectrum obj of given object and given visit
			'fiber' : fiber number of ap1d spectrum
	"""

	#frameid default to 1, where does this come from?--should investigate
	frameid = kwargs.get('frameid', 1) 

	ap_id, plates, mjds, fibers = ap.searchVisits(id_name=apogee_id)
	fiberid = fibers[visit]

	# url   = get_1dspec_urls(apogee_id, nvisit=visit, frame_id=frameid)
	# fiber = url['visit%s'%(visit)]['fiberid']-1

	data = []

	return data, fiberid


def makeModel(**kwargs):

	"""
	Input:  'params' : dictionary of parameters specified like {'teff': 3051, 'logg': 5.2, 'z': -0.25, 'vsini': 10., 'rv': -12, 'alpha': 0.8}

	Output: 'tell_sp' : spectrum object synthesized to input paramters
	"""

	params = kwargs.get('params')
	fiber  = kwargs.get('fiber', 40)
	plot   = kwargs.get('plot', False)

	labels = [params['teff'], params['logg'], params['z']]

	#Interpolate model grids at give teff, logg, fe/h
	interp_sp = ap.interpolateGrid(labels=labels)

	#Apply radial velocity
	rv_sp   = ap.spec_tools.rvShiftSpec(interp_sp, rv=params['rv'])

	#Apply rotational velocity broadening
	rot_sp  = ap.applyVsini(rv_sp, vsini=params['vsini'])

	#Apply APOGEE LSF function
	lsf_sp  = ap.convolveLsf(rot_sp, fiber=fiber)

	#Apply telluric spectrum
	tell_sp = ap.applyTelluric(lsf_sp)

	if plot == True:
		xrange = kwargs.get('xrange', [interp_sp.wave[0], interp_sp.wave[-1]])
		yrange = kwargs.get('yrange', [0.7, 1.1])

		plt.figure(1, figsize=(10,6))  
		plt.plot(interp_sp.wave, interp_sp.flux, label=r'Teff = %s, logg = %s, Fe/H = %s'%(params['teff'], params['logg'], params['z']))
		plt.plot(rv_sp.wave, rv_sp.flux, label=r'RV (%s km/s)'%(params['rv']))
		plt.plot(rot_sp.wave, rot_sp.flux, label=r'RV + rot (%s km/s)'%(params['vsini']))
		plt.plot(lsf_sp.wave, lsf_sp.flux, label=r'RV + rot + lsf')
		plt.plot(tell_sp.wave, tell_sp.flux, label=r'RV + rot + lsf + telluric')
		plt.xlim(xrange)
		plt.ylim(yrange)
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
	fiber  = kwargs.get('fiber', 40)

	synth_mdl = makeModel(params=params, fiber=fiber)
	chi = ap.compareSpectra(data, synth_mdl)

	return chi

