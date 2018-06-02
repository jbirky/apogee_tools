import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import emcee
import apogee_tools as ap
import corner
import argparse
import os


# =============================================================
# MCMC prior and likelihood functions
# =============================================================

def lnlike(theta, lsf, tell_sp):

	"""
	Log-likelihood, computed from chi-squared
	"""

	# If theta is entered as a list, make it into a dictionary
	theta_keys = [key for key in ap.init.keys()]
	if type(theta) == np.ndarray:
		theta = dict(zip(theta_keys, theta))

	# Choose the appropriate Spectrum class to read the data
	if ap.data['instrument'] == 'APOGEE':
		data = ap.Apogee(id=ap.data['ID'], type=ap.data["dtype"], visit=ap.data['visit'])
	else:
		print('No Spectrum class to read data for instrument', ap.data['instrument'])

	chisq = ap.returnModelFit(data, theta, lsf=lsf, telluric=tell_sp)

	print('\n chisq', chisq, '\n')

	return -0.5 * chisq


def lnprior(theta):

	"""
	Specifies a flat prior
	"""

	# keys = theta.keys()
	theta_keys = [key for key in ap.init.keys()]
	if type(theta) == np.ndarray:
		theta = dict(zip(theta_keys, theta))
	keys = theta.keys()

	for k in keys:
		if (ap.prior[k][0] < theta[k] < ap.prior[k][1]):
			pass
		else:
			return -np.inf
			break
	return 0.0


def lnprob(theta, lsf, tell_sp):

	lnp = lnprior(theta)
	if not np.isfinite(lnp):
	    return -np.inf

	return lnp + lnlike(theta, lsf, tell_sp)


#########################################################################################

"""
Conventions:
param - full list of parameters
theta - only parameters that are being sampled
"""

if __name__ == "__main__":

	if 'config.yaml' not in os.listdir():
		print('\nError: config.yaml not found in the current working directory. \
			Using default file found inside apogee_tools.\n')

	# =============================================================
	# Command line input
	# =============================================================

	parser = argparse.ArgumentParser(description='Specify plotting directory.')
	parser.add_argument("plot", action="store", type=str)
	args = parser.parse_args()


	# =============================================================
	# Testing...
	# =============================================================

	if 'make_model' in args.plot:

		init_param, step_param, init_theta, step_theta, fiber, tell_sp, lsf = ap.initialize()
		mdl = ap.makeModel(params=init_param, lsf=lsf, telluric=tell_sp, plot=True)

	if 'test_fit' in args.plot:

		data = ap.Apogee(id=ap.data['ID'], type=ap.data["dtype"], visit=ap.data['visit'])
		chi_sq = ap.returnModelFit(data, init_param, lsf=lsf, plot=True)

		print('chi^2', chi_sq)

	if 'test_telluric' in args.plot:

		mdl = ap.getModel(params=[3200, 5.0, 0.0], grid='BTSETTL', xrange=[15200,16940])
		tell_sp = ap.applyTelluric(mdl, ap.getTelluric(), cut_rng=[min(mdl.wave), max(mdl.wave)])

		tell_sp.plot()


	# =============================================================
	# Run MCMC!
	# =============================================================

	if (ap.out["mcmc_sampler"] == True) or ('mcmc' in args.plot):

		init_param, step_param, init_theta, step_theta, fiber, tell_sp, lsf = ap.initialize()

		theta_keys = list(init_theta.keys())
		theta_vals = list(init_theta.values())

		ndim = len(init_theta)
		nsteps = ap.mcmc["nsteps"]
		nwalkers = ap.mcmc["nwalkers"]

		pos = [list(init_theta.values()) + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(lsf, tell_sp))
		sampler.run_mcmc(pos, nsteps)

		np.save('sampler_chain', sampler.chain[:, :, :])

		samples = sampler.chain[:, :, :].reshape((-1, ndim))

		np.save('samples', samples)


	# =============================================================
	# Output corner/walker plots
	# =============================================================

	lbl = ['Teff', 'logg', '[Fe/H]', 'rv', 'vsini', r'$\alpha$']

	if ap.out["corner"] == True:

		try:
			fig = corner.corner(samples, labels=lbl, truths=theta_vals)
			fig.savefig("triangle.png")

		except:
			print('Traingle plot failed.')


	if 'corner' in args.plot:

		try:
			samples = np.load('samples.npy')
			fig = corner.corner(samples, labels=lbl, truths=theta_vals)
			fig.savefig("triangle.png")

		except:
			print('Traingle plot failed.')


	if (ap.out["walkers"] == True) or ('walkers' in args.plot):

		try:
			sampler_chain = np.load('sampler_chain.npy')

			ndim = len(sampler.chain.T)

			fig, ax = plt.subplots(ndim, sharex=True, figsize=[8,12])
			for i in range(ndim):
				ax[i].plot(sampler_chain.T[i], '-k', alpha=0.2);
				ax[i].set_ylabel(str(lbl[i]))
				if i == ndim:
					ax[i].set_xlabel(step)
			plt.tight_layout()
			plt.savefig('Walkers.png', dpi=300, bbox_inches='tight')
			plt.show()
			plt.close()

		except:
			print('Corner plot failed.')

