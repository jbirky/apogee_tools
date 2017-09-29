import numpy
from functools import wraps
import apogee_tools as ap
import apogee_tools.apogee_hack.spec.lsf as lsf
import apogee_tools.rotbroaden as rot
from apogee_tools.apogee_hack.spec.plot import apStarWavegrid
import matplotlib.pyplot as plt

PLOT = True # Set this to True if you want to plot the results

#### Need something here to pull a 1-d spectrum

#### Need something here to pull a 1-d spectrum


### Work on the model
# Get a model
mdl      = ap.getModel(params=[3200, 5.0, 0.0], grid='BTSETTLb', xrange=[15200,16940])

# shift the RV
rv_wave  = ap.spec_tools._rvShift(mdl.wave, rv=10)

# Rotation broaden
newspec1 = rot.broaden(rv_wave, mdl.flux, vbroad=10)

# Apply the LSF (This is borrowed from Jo Bovy's code)
fiber    = 40 # APOGEE Fiber. I just picked an APOGEE fiber. In principle this will come from the 1-d spectrum
xlsf     = numpy.linspace(-7.,7.,43)
lsf1     = lsf.eval(xlsf, fiber=fiber)
newspec2 = lsf.convolve(rv_wave, newspec1, xlsf=xlsf, lsf=lsf1, vmacro=None)

if PLOT:
	plt.plot(mdl.wave, mdl.flux, label='original', alpha=0.5)
	plt.plot(rv_wave, mdl.flux, label='RV (10 km/s)', alpha=0.5)
	plt.plot(rv_wave, newspec, label='RV+rot (10 km/s)', alpha=0.5)
	plt.plot(apStarWavegrid(), newspec2[0], label='RV+rot+lsf', alpha=0.5)
	plt.xlim(15678, 15694)
	plt.ylim(0.48, 0.60)
	plt.legend(frameon=False)
	#plt.savefig('Test.png', dpi=600)
	plt.show()

