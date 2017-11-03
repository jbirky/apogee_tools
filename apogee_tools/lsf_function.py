import numpy as np
import apogee_tools as ap
import apogee_tools.apogee_hack.spec.lsf as lsf
import apogee_tools.rotbroaden as rot
from apogee_tools.apogee_hack.spec.plot import apStarWavegrid


def convolveLsf(spec, **kwargs):

	"""
	Input:  'spec'  : spectrum object 
			'fiber' : APOGEE fiber. Find this from 1D spectrum

	Output: 'lsf_spec' : lsf convolved spectrum object
	"""

	fiber = kwargs.get('fiber', 40)

	xlsf = np.linspace(-7., 7., 43)
	lsf  = lsf.eval(xlsf, fiber=fiber)

	lsf_flux = lsf.convolve(spec.wave, spec.flux, xlsf=xlsf, lsf=lsf, vmacro=None)
	lsf_spec = ap.Spectrum(wave=spec.wave, flux=lsf_flux, name=spec.name)

	return lsf_spec