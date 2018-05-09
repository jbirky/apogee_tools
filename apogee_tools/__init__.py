from __future__ import absolute_import, division, print_function, unicode_literals

from .core import Spectrum, Apogee, Nirspec, ModelGrid

from .cannon_tools.run_cannon import labelToSpec, loadLabels, initializeTrainingSet, synthesizeFlux, \
	runCannon, fitCannonModel, crossValidate, _getPivotsAndScales, _get_lvec, scaleLabels
from .cannon_tools.plot_cannon import plotCrossValidation, plotCannonModels, plotSpectralSequence, plotBands

from .forward_model.lsf_function import convolveLsf
from .forward_model.model import initialize, makeModel, returnModelFit
from .forward_model.rotation_broaden import broaden, applyVsini
from .forward_model.rv_function import rvShift, rvShiftSpec
from .forward_model.synthesize_grid import interpolateGrid
from .forward_model.telluric import applyTelluric, getTelluric
# from .forward_model.read_model import *

from .info.features import lines

from .utils.ap1d import get_1dspec_urls, coadd_epoch, coadd_spectra
from .utils.continuum import continuum
from .utils.read import HDF5Convert, HDF5Interface, loadGrid, readModels, getModel
from .utils.read_lines import listLibraries, listSpecies, searchLines
from .utils.search import searchStars, searchVisits, download, multiParamSearch, returnSimbadParams, returnAspcapTable
from .utils.spec_tools import calcScale, compareSpectra, subtractContinuum, integralResample, splineInterpolate

from .instrument_tools.instrument import formatDesignation, getShortname, apogeeFile

import os
import yaml 


# Get directories of apogee_tools 
PATH 	  = os.path.realpath(__file__)
AP_TOOLS  = os.path.split(os.path.split(PATH)[0])[0]
LIBRARIES = AP_TOOLS + '/libraries'


# Read configuration file
try: # open from the current working directory
	f = open("config.yaml")
	config = yaml.load(f)
	f.close()

except: # otherwise read file from inside apogee_tools
	conf_path = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0] + '/config.yaml'

	f = open(conf_path)
	config = yaml.load(f)
	f.close()

data = config["data"]
workdir = config["workdir"]
fix_param = config["fix_param"]
out = config["out"]

model = config["model"]
mcmc = config["mcmc"]
init = config["init"]
step = config["step"]
prior = config["prior"]

instrument = data["instrument"]


# import .apogee_hack.spec.lsf as lsf
# from .apogee_hack.spec.plot import apStarWavegrid