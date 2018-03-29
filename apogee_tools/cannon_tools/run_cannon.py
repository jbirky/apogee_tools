import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
from astropy.io import ascii
from pathlib import Path
import apogee_tools as ap

try:
	from TheCannon import apogee, dataset, model
except:
	print("apogee_tools: Optional dependency TheCannon is not installed. To use cannon_tools, run 'pip install TheCannon' \n")

import os
os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2016/bin/x86_64-darwin'


def loadLabels(filename, **kwargs):
	
    data = ascii.read(filename)
    ids  = data['ID']
    inds = ids.argsort()
    ids  = ids[inds]

    lbl_names = kwargs.get('lbl_names', ['SPT'])

    nlbl = len(lbl_names)
    params = []

    for i in range(nlbl):
    	p = data[lbl_names[i]]
    	p = p[inds]
    	params.append(p)
    
    tr_label = np.vstack(params).T
    
    return tr_label


def initializeTrainingSet(**kwargs):

	"""
	Read in training data set for The Cannon, and normalize the spectra

	Input:  'data' : director of folder containing aspcap data
		    'ref'  : csv file of reference labels

	Output: ids, wl, tr_flux, tr_ivar, tr_label = training_set_info

			'ids'      : spectrum ids
			'wl'       : wavelength array 
			'tr_flux'  : array of training fluxes
			'tr_ivar'  : array of inverse variances
			'tr_label' : array of training set labels

			ds : dataset object for training set spectra
	"""

	# required
	data_dir  = kwargs.get('data') 
	lbl_file  = kwargs.get('ref', None)

	# optional
	e_vers = kwargs.get('eiler', False) # if True, will use Eiler's Cannon settings
	uncert = kwargs.get('uncert', False) 
	save_base = kwargs.get('save_base') #specifies name normalized flux will be save as
	ds_ranges = kwargs.get('ds_ranges', [[371,3192], [3697,5997], [6461,8255]])

	# main directory of cannon run
	main_dir = Path(data_dir).parent

	# optional:
	lbl_names = kwargs.get('lbl_names', ['SPT'])

	tr_ID, wl, tr_flux, tr_ivar = apogee.load_spectra(data_dir)

	#Labels can be none if this is a test dataset, instead of training dataset
	if lbl_file == None:
		tr_label = None
	else:
		tr_label = loadLabels(lbl_file, lbl_names=lbl_names)

	test_ID = tr_ID
	test_flux = tr_flux
	test_ivar = tr_ivar

	if uncert == False:
		tr_delta    = None
		coeff_old   = None
		scatter_old = None

	ids = ['2M'+ID.split('2M')[1].split('.fits')[0] for ID in tr_ID]

	if e_vers == True: # If using Eiler's Cannon
		ds = dataset.Dataset(wl, tr_ID, tr_flux, tr_ivar, tr_label, tr_delta, test_ID, test_flux, test_ivar, coeff_old=coeff_old, scatter_old=scatter_old)
	else: 
		ds = dataset.Dataset(wl, tr_ID, tr_flux, tr_ivar, tr_label, test_ID, test_flux, test_ivar)
	ds.set_label_names(lbl_names)
	ds.ranges = ds_ranges

	if not os.path.exists(str(main_dir) + '/norm_fluxes/'):
		os.makedirs(str(main_dir) + '/norm_fluxes/')

	try:
		save_flux = str(main_dir) + '/norm_fluxes/' + save_base + '_norm_tr_flux.npy'
		save_ivar = str(main_dir) + '/norm_fluxes/' + save_base + '_norm_tr_ivar.npy'
	except:
		save_flux = str(main_dir) + '/norm_fluxes/' + '_norm_tr_flux.npy'
		save_ivar = str(main_dir) + '/norm_fluxes/' + '_norm_tr_ivar.npy'

	# Continuum normalization
	try:

		norm_tr_flux = np.load(save_flux)
		norm_tr_ivar = np.load(save_ivar)
		norm_test_flux = norm_tr_flux
		norm_test_ivar = norm_tr_ivar

		print('Loaded', save_flux)
		print('Loaded', save_ivar)

	except:

		pseudo_tr_flux, pseudo_tr_ivar = ds.continuum_normalize_training_q(\
			    q=0.90, delta_lambda=50)

		contmask = ds.make_contmask(pseudo_tr_flux, pseudo_tr_ivar, frac=0.07)
		ds.set_continuum(contmask)

		cont = ds.fit_continuum(3, "sinusoid")

		norm_tr_flux, norm_tr_ivar, norm_test_flux, norm_test_ivar = \
				ds.continuum_normalize(cont)

		np.save(save_flux, norm_tr_flux)
		np.save(save_ivar, norm_tr_ivar)

		print('Saved', save_flux)
		print('Saved', save_ivar)

	ds.tr_ID = ids
	ds.tr_flux = norm_tr_flux
	ds.tr_ivar = norm_tr_ivar
	ds.test_flux = norm_test_flux
	ds.test_ivar = norm_test_ivar

	return ds


def runCannon(ds, **kwargs):

	"""
	Run The Cannon for Teff and [Fe/H]. See https://annayqho.github.io/TheCannon/apogee_tutorial.html#apogee-tutorial

	Input:  'ds' : dataset object output from initializeTrainingSet()

	Output: wavelength, IDs, training flux, training labels,
			synthesized flux, synthesized labels
	"""

	# Fit model
	md = model.CannonModel(2, None)
	md.fit(ds)

	# Infer labels
	label_errs  = md.infer_labels(ds)
	test_labels = ds.test_label_vals

	coeffs = md.coeffs
	scatter = md.scatters
	test_labels = test_labels

	if ds.tr_label.shape[1] == 1:
		par1 = test_labels.T[0]
		
		pivot1, scale1 = _getPivotsAndScales(par1)

		lbl1 = [(t - pivot1)/scale1 for t in par1]
		lbl1_sq = [t**2 for t in lbl1]

		# Create label vector
		label_n = []
		for i in range(test_labels.shape[0]):
		    l_i = [1, lbl1[i], lbl1_sq[i]]
		    label_n.append(l_i)
		label_n = np.array(label_n)

	elif ds.tr_label.shape[1] == 2:
		par1 = list(map(itemgetter(0), test_labels))
		par2 = list(map(itemgetter(1), test_labels))

		pivot1, scale1 = _getPivotsAndScales(par1)
		pivot2, scale2 = _getPivotsAndScales(par2)

		lbl1 = [(t - pivot1)/scale1 for t in par1]
		lbl2 = [(m - pivot2)/scale2 for m in par2]
		lbl1_sq = [t**2 for t in lbl1]
		lbl2_sq = [f**2 for f in lbl2]
		lbl1_lbl2 = [lbl2[i]*lbl1[i] for i in range(len(par1))]

		# Create label vector
		label_n = []
		for i in range(test_labels.shape[0]):
		    l_i = [1, lbl1[i], lbl2[i], lbl1_sq[i], lbl2_sq[i], lbl1_lbl2[i]]
		    label_n.append(l_i)
		label_n = np.array(label_n)

	synth_fluxes = np.dot(coeffs, np.transpose(label_n))
	synth_fluxes = np.transpose(synth_fluxes)

	return md, synth_fluxes, test_labels


def crossValidate(ds, **kwargs):

    """
    Cross-validation test of Cannon test spectra
    Input:  'ds' : dataset object output from initializeTrainingSet()
    Output: plot training label vs. inferred label left out of set
            return training label and inferred label
    """

    # optional
    label_names = kwargs.get('lbl_names', ['Teff', 'Fe/H'])
    save_dir = kwargs.get('save_dir', 'cross_validation/')

    # ds: Data set of all objects (including n); continuum normalized in initializeTrainingSet()
    wl = ds.wl
    N = len(ds.tr_ID)

    # Training labels and cross-validated labels
    trn_labels, crv_labels = [], []
    
    ds.set_label_names(label_names) 
    md, synth_fluxes, test_labels = ap.runCannon(ds)
    np.save(save_dir+'model_coeff_fullset', md.coeffs)
    np.save(save_dir+'synth_flux_fullset', synth_fluxes)
    np.save(save_dir+'test_labels_fullset', test_labels)

    # Remove nth spectrum from the training set and run Cannon model on the rest of spectra
    for n in range(N): 

        tr_ids   = list(ds.tr_ID)
        tr_flux  = list(ds.tr_flux)
        tr_ivar  = list(ds.tr_ivar)
        tr_label = list(ds.tr_label)

        tr_ids.pop(n)
        tr_flux.pop(n)
        tr_ivar.pop(n)
        tr_label.pop(n)

        id_minus_n = np.array(tr_ids)
        fl_minus_n = np.array(tr_flux)
        vr_minus_n = np.array(tr_ivar)
        tr_label_minus_n = np.array(tr_label)

        ds_minus_n = dataset.Dataset(wl, id_minus_n, fl_minus_n, vr_minus_n, tr_label_minus_n, \
            id_minus_n, fl_minus_n, vr_minus_n)

        ds_minus_n.set_label_names(label_names)       
        md_minus_n, synth_fl_minus_n, te_label_minus_n = runCannon(ds_minus_n)

        # Find cross-validation label
        crv_label_n = test_labels[n]
        crv_labels.append(crv_label_n)
        
        np.save(save_dir+'model_coeff_'+str(n), md_minus_n.coeffs)
        np.save(save_dir+'synth_flux_'+str(n), synth_fl_minus_n)
        np.save(save_dir+'test_labels_'+str(n), te_label_minus_n)
        np.save(save_dir+'crv_label_'+str(n), crv_label_n)

        print('Labeled [%s/%s] sources.\n'%(n+1, N))

    trn_labels = ds.tr_label

    return trn_labels, crv_labels


def labelToSpec(labels, coeffs):

	nlabels = labels.shape[1]

	scaled_labels = []
	for i in range(nlabels):
		p, s = _getPivotsAndScales(labels.T[i])
		slbl = [(t - p)/s for t in labels.T[i]]
		scaled_labels.append(slbl)

	label_vec = np.array([_get_lvec(lbl) for lbl in np.array(scaled_labels).T])
	label_vec = np.insert(label_vec, 0, 1, axis=1)

	synth_fluxes = np.dot(coeffs, label_vec.T).T

	return synth_fluxes


def synthesizeFlux(ds, **kwargs):

	order = kwargs.get('order', 2)

	md = model.CannonModel(order, None)
	md.fit(ds)

	md.infer_labels(ds)
	test_labels = ds.test_label_vals

	synth_fluxes = labelToSpec(test_labels, md.coeffs)

	return ds, synth_fluxes


def fitCannonModel(ds, **kwargs):
    
    order = kwargs.get('order', 2)
    
    md = model.CannonModel(order, None)
    md.fit(ds)
    
    coeffs = md.coeffs
    md.infer_labels(ds)
    test_labels = ds.test_label_vals
    
    nlabels = test_labels.shape[1]
    
    scaled_labels = []
    for i in range(nlabels):
        p, s = _getPivotsAndScales(test_labels.T[i])
        slbl = [(t - p)/s for t in test_labels.T[i]]
        scaled_labels.append(slbl)
    
    label_vec = np.array([_get_lvec(lbl) for lbl in np.array(scaled_labels).T])
    label_vec = np.insert(label_vec, 0, 1, axis=1)
    
    synth_fluxes = np.dot(coeffs, label_vec.T).T
    
    return md, ds, synth_fluxes


def _getPivotsAndScales(label_vals):

	"""
	Function scales the labels see https://github.com/annayqho/TheCannon
	"""

	qs = np.percentile(label_vals, (2.5, 50, 97.5), axis=0)
	pivots = qs[1]
	scales = (qs[2] - qs[0])/4.

	return pivots, scales


def _get_lvec(labels):

	"""
	Return quadratic label vector, see https://github.com/annayqho/TheCannon
	"""

	nlabels = len(labels)
	# specialized to second-order model
	linear_terms = labels 
	quadratic_terms = np.outer(linear_terms, linear_terms)[np.triu_indices(nlabels)]
	lvec = np.hstack((linear_terms, quadratic_terms))
	
	return lvec


def scaleLabels(labels):
    
    scaled_labels = []
    for i in range(labels.shape[1]):
        p, s = _getPivotsAndScales(labels.T[i])
        slbl = [(t - p)/s for t in labels.T[i]]
        scaled_labels.append(slbl)
        
    return np.array(scaled_labels).T


def interpolateGrids(**kwargs):

	"""
	Interpolate a set of grids using The Cannon; return spectrum objects
	"""

	# required
	prange = kwargs.get('prange', [[3000,4000], [4.0,5.5], [-0.5,0.5]])
	irange = kwargs.get('irange', prange) #interpolation range--must be subset of prange
	grid   = kwargs.get('grid', 'PHOENIX')

	# optional
	step   = kwargs.get('step', [100, .5, .5]) #interpolation step size
	deg    = kwargs.get('deg', 2) #degree of cannon polynomial fit
	sp_dir = kwargs.get('sp_dir', 'normalized_specs/')
	xrange = kwargs.get('xrange', [15200,16940])
	lbl_names = kwargs.get('lbl_names', ['T_{eff}', '\log g', '[Fe/H]'])

	save_name = '%s%s_t%s_%s_l%s_%s_z%s_%s_' %(sp_dir, grid, prange[0][0], prange[0][1], prange[1][0], prange[1][1], prange[2][0], prange[2][1])
	save_flux = save_name + 'flux'
	save_ivar = save_name + 'ivar'

	teffs = np.arange(prange[0][0], prange[0][1]+10, 100)
	loggs = np.arange(prange[1][0], prange[1][1]+.1, 0.5)
	fe_hs = np.arange(prange[2][0], prange[2][1]+.1, 0.5)

	param_list = [teffs, loggs, fe_hs]
	params = np.array([list(x) for x in np.array(np.meshgrid(*param_list)).T.reshape(-1,len(param_list))])

	try:

		tr_flux = np.load(save_flux + '.npy')
		tr_ivar = np.load(save_ivar + '.npy')
		wl = ap.getModel(params=params[0], grid=grid, xrange=xrange).wave
		tr_label = params
		tr_ID = tr_label

		print('Loaded', save_flux + '.npy')
		print('Loaded', save_ivar + '.npy')

		ds = dataset.Dataset(wl, tr_ID, tr_flux, tr_ivar, tr_label, tr_ID, tr_flux, tr_ivar)
		ds.set_label_names(lbl_names)

	except:

		tr_ID, tr_flux, tr_ivar = [], [], []
		for par in params:
			mdl = ap.getModel(params=par, grid=grid, xrange=xrange)
			tr_flux.append(np.array(mdl.flux))
			tr_ID.append(par)
			
			ivar = [100 for i in mdl.wave]
			tr_ivar.append(ivar)
		wl = mdl.wave
		tr_label = np.array(tr_ID)
		tr_flux = np.array(tr_flux)
		tr_ivar = np.array(tr_ivar)

		np.save('normalized_specs/phn_wl.npy', wl)

		ds = dataset.Dataset(wl, tr_ID, tr_flux, tr_ivar, tr_label, tr_ID, tr_flux, tr_ivar)
		ds.set_label_names(lbl_names)

		pseudo_tr_flux, pseudo_tr_ivar = ds.continuum_normalize_training_q(q=0.90, delta_lambda=50)
		contmask = ds.make_contmask(pseudo_tr_flux, pseudo_tr_ivar, frac=0.07)
		ds.set_continuum(contmask)
		cont = ds.fit_continuum(3, "sinusoid")
		norm_tr_flux, norm_tr_ivar, norm_test_flux, norm_test_ivar = ds.continuum_normalize(cont)
		ds.tr_flux = norm_tr_flux
		ds.tr_ivar = norm_tr_ivar
		ds.test_flux = norm_test_flux
		ds.test_ivar = norm_test_ivar

		np.save(save_flux, norm_tr_flux)
		np.save(save_ivar, norm_tr_ivar)

		print('Saved ', save_flux)
		print('Saved ', save_ivar)

	md = model.CannonModel(deg, None)
	md.fit(ds)

	coeffs = md.coeffs

	#Interpolation labels
	i_teffs = np.arange(irange[0][0], irange[0][1]+ 10, step[0])
	i_loggs = np.arange(irange[1][0], irange[1][1]+.01, step[1])
	i_fe_hs = np.arange(irange[2][0], irange[2][1]+.01, step[2])

	iparam_list = [i_teffs, i_loggs, i_fe_hs]
	iparams = np.array([list(x) for x in np.array(np.meshgrid(*iparam_list)).T.reshape(-1,len(iparam_list))])

	synth_fluxes = labelToSpec(iparams, coeffs)

	synth_specs = []
	for i in range(len(synth_fluxes)):
		sp = ap.Spectrum(wave=ds.wl, flux=synth_fluxes[i], name=iparams[i])
		synth_specs.append(sp)

	if kwargs.get('save_specs', False) == True:
		np.save('normalized_specs/PHOENIX_synth_fluxes', synth_fluxes)
		np.save('normalized_specs/PHOENIX_synth_wave', ds.wl)
		np.save('normalized_specs/PHOENIX_synth_coeffs', coeffs)

	return synth_specs


