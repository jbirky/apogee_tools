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

	print('\n theta', theta, '\n')

	# mdl  = ap.makeModel(params=theta, fiber=124)
	data = ap.Spectrum(id=ap.data['ID'], type='ap1d', visit=ap.data['visit'])

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

	init_param, step_param, init_theta, step_theta, fiber = ap.initialize()

	theta_keys = list(init_theta.keys())
	theta_vals = list(init_theta.values())

	ndim = len(init_theta)
	nsteps = ap.mcmc["nsteps"]
	nwalkers = ap.mcmc["nwalkers"]

	print(lnlike(init_theta))

	# mdl = makeModel(params=init_theta, fiber=fiber)

	# pos = [list(init_theta.values()) + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

	# sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike)
	# sampler.run_mcmc(pos, nsteps)

	# samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

	# fig = corner.corner(samples, labels=theta_keys, truths=theta_vals)
	# fig.savefig("triangle.png")

	# fit = ap.fitMCMC(init_par, step_par, fiber)