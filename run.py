import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import emcee
import apogee_tools as ap
import corner


def lnlike(theta):

	theta_keys = [key for key in ap.init.keys()]
	if type(theta) == np.ndarray:
		theta = dict(zip(theta_keys, theta))

	# mdl  = ap.makeModel(params=theta, fiber=124)
	data = ap.Spectrum(id=ap.data['ID'], type=ap.data["dtype"], visit=ap.data['visit'])

	chisq = ap.returnModelFit(data, theta, fiber=124)

	print('\n chisq', chisq, '\n')

	return -0.5 * chisq


def lnprior(theta):

	keys = theta.keys()

	for k in keys:
		if (ap.prior[k][0] < theta[k] < ap.prior[k][1]):
			pass
		else:
			return -np.inf
			break
	return 0.0


def lnprob(theta):

	lnp = lnprior(theta)
	if not np.isfinite(lp):
	    return -np.inf

	return lnp + lnlike(theta, x, y, yerr)


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
	# Pre-MCMC initialization
	# =============================================================

	init_param, step_param, init_theta, step_theta, fiber, tell_sp = ap.initialize()

	theta_keys = list(init_theta.keys())
	theta_vals = list(init_theta.values())

	ndim = len(init_theta)
	nsteps = ap.mcmc["nsteps"]
	nwalkers = ap.mcmc["nwalkers"]


	# =============================================================
	# Run MCMC!
	# =============================================================

	if ap.out["mcmc_sampler"] == True:

		pos = [list(init_theta.values()) + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike)
		sampler.run_mcmc(pos, nsteps)

		np.save('sampler_chain', sampler.chain[:, :, :])

		samples = sampler.chain[:, :, :].reshape((-1, ndim))

		np.save('samples', samples)


	# =============================================================
	# Output corner/walker plots
	# =============================================================

	if ap.out["corner"] == True:

		try:
			fig = corner.corner(samples, labels=theta_keys, truths=theta_vals)
			fig.savefig("triangle.png")

		except:
			print('Traingle plot failed.')


	if ap.out["walkers"] == True:

		try:
			sampler_chain = np.load('sampler_chain.npy')
			lbl = ['Teff', 'logg', '[Fe/H]', 'rv', 'vsini', r'$\alpha$']

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
