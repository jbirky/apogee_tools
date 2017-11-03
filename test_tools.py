import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='serif')
import apogee_tools as ap
import os

if __name__ == '__main__':

#--------------------------------------------------------------------
# Downloading

	# ap.download('2M15141711+0044474', type='aspcap')
	# data = ap.Spectrum(id='2M15141711+0044474', type='aspcap')
	# data.plot(items=['spectrum', 'apModel'], save=True)

#--------------------------------------------------------------------
# Searching

	# params = ['TEFF', 'LOGG', 'M_H']
	# ranges = [[-10000,4000], [0,5], [-2,2]]
	# source_table = ap.multiParamSearch(par=params, select=ranges)

#--------------------------------------------------------------------
# Test Model functions

	params = {'teff': 3051, 'logg': 5.2, 'z': -0.25, 'vsini': 10., 'rv': -12, 'alpha': 0.8}

# Replicate Chris T's test script

	# labels = [params['teff'], params['logg'], params['z']]

	# interp_sp = ap.interpolateGrid(labels=labels)

	# rv_sp   = ap.spec_tools.rvShiftSpec(interp_sp, rv=params['rv'])

	# rot_sp  = ap.applyVsini(rv_sp, vsini=params['vsini'])

	# lsf_sp  = ap.convolveLsf(rot_sp, fiber=40)

	# tell_sp = ap.applyTelluric(lsf_sp)

	# plt.figure(1, figsize=(10,6))  
	# plt.plot(interp_sp.wave, interp_sp.flux, label=r'Teff = %s, logg = %s, Fe/H = %s'%(params['teff'], params['logg'], params['z']))
	# plt.plot(rv_sp.wave, rv_sp.flux, label=r'RV (%s km/s)'%(params['rv']))
	# plt.plot(rot_sp.wave, rot_sp.flux, label=r'RV + rot (%s km/s)'%(params['vsini']))
	# plt.plot(lsf_sp.wave, lsf_sp.flux, label=r'RV + rot + lsf')
	# plt.plot(tell_sp.wave, tell_sp.flux, label=r'RV + rot + lsf + telluric')
	# plt.xlim(15678, 15694)
	# plt.ylim(0.7, 1.1)
	# plt.legend(frameon=False)
	# plt.show()

# Test makeModel function (found in model.py), which is the module containing the test code above

	ap.makeModel(params=params, fiber=40, plot=True, xrange=[15678,15694])

