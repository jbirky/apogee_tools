from __future__ import absolute_import, division, print_function, unicode_literals

from .core import Spectrum

from .cannon_tools.run_cannon import labelToSpec, loadLabels, initializeTrainingSet, synthesizeFlux, \
	runCannon, fitCannonModel, crossValidate, _getPivotsAndScales, _get_lvec
from .cannon_tools.plot_cannon import plotCrossValidation, plotCannonModels, plotSpectralSequence, plotBands

from .forward_model.lsf_function import convolveLsf
from .forward_model.model import makeModel, returnModelFit
from .forward_model.rotation_broaden import broaden, applyVsini
from .forward_model.rv_function import rvShift, rvShiftSpec
from .forward_model.synthesize_grid import interpolateGrid
from .forward_model.telluric import applyTelluric

from .info.features import lines

from .utils.ap1d import get_1dspec_urls, coadd_epoch, coadd_spectra
from .utils.continuum import continuum
from .utils.read import HDF5Convert, HDF5Interface, loadGrid, readModels, getModel
from .utils.search import searchStars, searchVisits, download, multiParamSearch, returnSimbadParams, returnAspcapTable
from .utils.spec_tools import calcScale, compareSpectra, subtractContinuum, integralResample


# import .apogee_hack.spec.lsf as lsf
# from .apogee_hack.spec.plot import apStarWavegrid