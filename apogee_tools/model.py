import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# from apogee import Spectrum

#Get the path of apogee_tools file
FULL_PATH  = os.path.realpath(__file__)
BASE = os.path.split(os.path.split(FULL_PATH)[0])[0]

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

	4. addTelluric()
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
		# self.data   = kwargs.get('data')
		self.params = kwargs.get('params')

		par_keys = ['teff', 'logg', 'fe_h', 'rv', 'vsini']

		self.teff  = self.params['teff']
		self.logg  = self.params['logg']
		self.fe_h  = self.params['fe_h']
		self.rv    = self.params['rv']
		self.vsini = self.params['vsini']
		# self.alpha = self.params['alpha']

		self.wave = np.load(BASE+'/libraries/cannon_phoenix/phn_cannon_wl.npy')


	def interpolateGrid(self):

		"""
		Input set of labels, dot product with Cannon model coefficients, return set of fluxes.
		"""

		labels  = np.array([[self.teff, self.logg, self.fe_h]])
		nlabels = labels.shape[1]

		pivots = np.load(BASE+'/libraries/cannon_phoenix/phn_cannon_pivots.npy')
		scales = np.load(BASE+'/libraries/cannon_phoenix/phn_cannon_scales.npy')
		coeffs = np.load(BASE+'/libraries/cannon_phoenix/phn_cannon_coeffs.npy')

		scaled_labels = [(lbl - pivots) / scales for lbl in labels]

		label_vecs = [list(_get_lvec(lbl)) for lbl in scaled_labels]
		label_vecs = np.column_stack(([1 for l in label_vecs], label_vecs))

		self.iflux = np.array(np.dot(coeffs, np.transpose(label_vecs)))
		self.iflux = np.transpose(self.iflux)[0]
		
		return self.iflux


	# def applyLSF(self):
		

	def addTelluric(self):

		"""
		Apply telluric model to PHOENIX model (or whatever grid you're using)
		"""

		mdl_obj = Spectrum(wave=self.wave, flux=self.flux)
		self.tflux = applyTelluric(mdl_obj)

		return self.tflux


	def fitModel(self, data):

		"""
		Wrap together all of the functions of this class.
		Returns chi-squared value between data, and model at given parameters.
		
		Input:  params : dictionary of parameters, stored with the following keys:
						  ex: {'teff': 1000, 'logg': 4.5, 'z':0.0, 'vsini': 30., 'rv': -12, 'alpha': 0.8}
		"""

		self.interpolateGrid()
		self.addTelluric()

		synth_model = Spectrum(wave=self.wave, flux=self.tflux, name=self.params)

		self.chi = compareSpectra(data, synth_model)[0]

		return self.chi, self.params


def _get_lvec(labels):

    """
    Constructs a label vector for an arbitrary number of labels
    Assumes that our model is quadratic in the labels
    @Anna Ho

    Parameters
    ----------
    labels: numpy ndarray
        pivoted label values for one star

    Returns
    -------
    lvec: numpy ndarray
        label vector
    """

    nlabels = len(labels)

    # specialized to second-order model
    linear_terms = labels 
    quadratic_terms = np.outer(linear_terms, linear_terms)[np.triu_indices(nlabels)]
    lvec = np.hstack((linear_terms, quadratic_terms))

    return lvec
