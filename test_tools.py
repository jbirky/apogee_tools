import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='serif')
import apogee_tools as ap
import os

if __name__ == '__main__':

#--------------------------------------------------------------------
# Downloading/Reading

	# ap.download('2M03425325+2326495', type='ap1d', visit=1, frame=1)
	# data = ap.Spectrum(id='2M03425325+2326495', type='ap1d', visit=1)
	# data.plot(items=['spec', 'noise'], yrange=[0, 2e-13])
	# data.mask(sigma=[1,.5], pixel_buffer=[0,3])
	# data.plot(items=['spec'], yrange=[0, 4e-14])

	# params = {'teff': 3051, 'logg': 5.2, 'z': -0.25, 'vsini': 10., 'rv': -12, 'alpha': .1}
	# mdl = ap.makeModel(params=params, fiber=40, xrange=[data.wave[0],data.wave[-1]])
	# mdl.plot(yrange=[.7,1.2])

#--------------------------------------------------------------------
# Searching

	# params = ['TEFF', 'LOGG', 'M_H']
	# ranges = [[-10000,4000], [0,5], [-2,2]]
	# source_table = ap.multiParamSearch(par=params, select=ranges)

#--------------------------------------------------------------------
# Test Model functions

	params = {'teff': 3100, 'logg': 5.5, 'z': -0.3, 'vsini': 0.1, 'rv': -41, 'alpha': .5}

# Replicate Chris T's test script

	interp_sp = ap.interpolateGrid(labels=[params['teff'], params['logg'], params['z']], res='500k')
	interp_sp.flux = interp_sp.flux/max(interp_sp.flux)
	rv_sp   = ap.spec_tools.rvShiftSpec(interp_sp, rv=params['rv'])
	rot_sp  = ap.applyVsini(rv_sp, vsini=params['vsini'])
	tell_sp = ap.applyTelluric(rot_sp, alpha=params['alpha'], airmass='1.0')
	lsf_sp  = ap.convolveLsf(tell_sp, fiber=124)

	plt.figure(1, figsize=(16,6))  
	plt.plot(interp_sp.wave, interp_sp.flux, alpha=.7, linewidth=1, label=r'Teff = %s, logg = %s, Fe/H = %s'%(params['teff'], params['logg'], params['z']))
	plt.plot(rv_sp.wave, rv_sp.flux-.15, label=r'RV (%s km/s)'%(params['rv']), alpha=.7, linewidth=1)
	plt.plot(rot_sp.wave, rot_sp.flux-.3, label=r'RV + rot (%s km/s)'%(params['vsini']), alpha=.7, linewidth=1)
	plt.plot(tell_sp.wave, tell_sp.flux-.45, label=r'RV + rot + telluric ($\alpha$ = %s)'%(params['alpha']), alpha=.7, linewidth=1)
	plt.plot(lsf_sp.wave, lsf_sp.flux-.8, label=r'RV + rot + telluric + lsf', alpha=.7, linewidth=1)
	plt.xlim([15900,16200]) #[15190,16950]
	plt.ylim(-.4, 1.1)
	plt.legend(loc='lower left', frameon=False, fontsize=12)
	plt.ylabel(r'$F_{\lambda}$ + offset', fontsize=15)
	plt.xlabel(r'$\lambda$', fontsize=15)
	plt.title('Teff = %s, logg = %s, Fe/H = %s'%(params['teff'], params['logg'], params['z']))
	plt.show()

# Test makeModel function (found in model.py), which is the module containing the test code above

	# ap.makeModel(params=params, fiber=40, plot=True, xrange=[15678,15694])

