import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

#Get the path of apogee_tools file
FULL_PATH  = os.path.realpath(__file__)
BASE, NAME = os.path.split(FULL_PATH)

AP_PATH = os.environ['APOGEE_DATA']

"""
Code modularization:

	0. Process data
		- read in, continuum normalize, mask, etc.

	1. interpolateGrid(teff, logg, fe_h)
		- Interpolate teff, logg, fe_h parameters from model grid
		- Use Cannon model coefficients; already continuum normalized

	2. _rvShift(wave, **kwargs)
		- Apply radial velocity shift to the model

	3. smoothVSINI(mspec, **kwargs)
		- Apply vsini

	4. applyTelluric()
		- Apply telluric correction to model

	5. compareSpectra(data, model)
		- calculate chi-squared fit

-------------------------------------------------
	fitModel(data, params)
		- Master function that combines 1-5 +MCMC
		- Returns chi-squared value between data, and model at given parameters


User input for code:

	md = Model(data=data, params=params)
	chi = md.fitModel()

"""

class Model():

	def __init__(self, **kwargs):

		# Required inputs:
		self.data   = kwargs.get('data')
		self.params = kwargs.get('params')

		par_keys = ['teff', 'logg', 'fe_h', 'rv', 'vsini', 'alpha']

		self.teff  = params['teff']
		self.logg  = params['logg']
		self.fe_h  = params['fe_h']
		self.rv    = params['rv']
		self.vsini = params['vsini']
		self.alpha = params['alpha']


	def labelToSpec(labels):

		"""
		Input set of labels, dot product with Cannon model coefficients, return set of fluxes.
		"""

		labels  = np.array(labels)
		nlabels = labels.shape[1]

		pivots = np.load(BASE+'/libraries/cannon_phoenix/phn_cannon_pivots.npy')
		scales = np.load(BASE+'/libraries/cannon_phoenix/phn_cannon_scales.npy')
		coeffs = np.load(BASE+'/libraries/cannon_phoenix/phn_cannon_coeffs.npy')

		scaled_labels = [(lbl - pivots) / scales for lbl in labels]

		label_vecs = [list(_get_lvec(lbl)) for lbl in scaled_labels]
		label_vecs = np.column_stack(([1 for l in label_vecs], label_vecs))

		interp_fluxes = np.dot(coeffs, np.transpose(label_vecs))
		interp_fluxes = np.transpose(synth_fluxes)
		
		return interp_fluxes


	def interplateGrid(self):

		"""
		Interpolate grids over teff, logg, fe_h using The Cannon.
		"""

		labels = [self.teff, self.logg, self.fe_h]

		self.interp_flux = labelToSpec(labels)


	def applyTelluric(self):

		"""
		Apply telluric model to PHOENIX model (or whatever grid you're using)
		"""


	def synthesizeModel(self):

		"""
		Input:  params : dictionary of parameters, stored with the following keys:
						  ex: {'teff': 1000, 'logg': 4.5, 'z':0.0, 'vsini': 30., 'rv': -12, 'alpha': 0.8}

		Output: synth_model : apogee_tools Spectrum object; 
							  model spectrum synthesized to given parameters
		"""



	def fitModel(data, params):

		"""
		Wrap together all of the functions of this class.
		Returns chi-squared value between data, and model at given parameters.
		"""

		interpolateGrid()
		applyTelluric()
		synthesizeModel()

		synth_model = Spectrum(wave=self.mdl_wave, flux=self.synth_flux, name=self.params)

		chi = compareSpectra(data, synth_model)[0]

		return chi, synth_model

