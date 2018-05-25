import apogee_tools as ap
import os, sys
from astropy.utils.data import download_file
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
from astropy import units as u
import numpy
import glob
import pandas
import copy

# structure to store models that have been read in
MODELS_READIN = {}

SPECTRAL_MODELS = {\
#    'gaia': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/gaia/', 'name': 'AMES GAIA', 'citation': 'Hauschildt et al. (1999)', 'bibcode': '1999ApJ...525..871H', 'altname': ['nextgen,hauschildt,hauschildt99,hauschildt1999'], 'rawfolder': HOME_FOLDER+'/models/phoenix/nextgen/fullres/', 'default': {'teff': 2000., 'logg': 5.0, 'z': 0.0}}, \
    'btnextgen': {'instruments': {}, 'name': 'BT NextGen', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['nextgen-bt','btnextgen'], 'default': {'teff': 3000., 'logg': 5.0, 'z': 0.0, 'enrich': 0.}}, \
    'btsettl08': {'instruments': {'APOGEE-RAW'}, 'name': 'BTSettl08', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['allard','allard12','allard2012','btsettl','btsettled','btsettl08','btsettl2008','BTSettl2008'], 'default': {'teff': 2000., 'logg': 5.0, 'z': 0., 'enrich': 0.}}, \
    'btsettl15': {'instruments': {}, 'name': 'BTSettl15', 'citation': 'Allard et al. (2015)', 'bibcode': '2015A&A...577A..42B', 'altname': ['allard15','allard2015','btsettl015','btsettl2015','BTSettl2015'],  'default': {'teff': 1500., 'logg': 5.0, 'z': 0.}}, \
    'burrows06': {'instruments': {}, 'name': 'Burrows 2006', 'citation': 'Burrows et al. (2006)', 'bibcode': '2006ApJ...640.1063B', 'altname': ['burrows','burrows2006'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'cld': 'nc'}}, \
    'cond01': {'instruments': {}, 'name': 'AMES Cond', 'citation': 'Allard et al. (2001)', 'bibcode': '2001ApJ...556..357A', 'altname': ['cond','cond-ames','amescond'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
    'drift': {'instruments': {}, 'name': 'Drift', 'citation': 'Witte et al. (2011)', 'bibcode': '2011A&A...529A..44W', 'altname': ['witte','witte11','witte2011','helling'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.}}, \
    'dusty01': {'instruments': {}, 'name': 'AMES Dusty', 'citation': 'Allard et al. (2001)', 'bibcode': '2001ApJ...556..357A', 'altname': ['dusty','dusty-ames','amesdusty'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
    'madhusudhan11': {'instruments': {}, 'name': 'Madhusudhan 2011', 'citation': 'Madhusudhan et al. (2011)', 'bibcode': '2011ApJ...737...34M', 'altname': ['madhu','madhusudhan','madhu11','madhu2011','madhusudhan2011'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.,'cld': 'ae60', 'kzz': 'eq','fsed': 'eq'}}, \
    'morley12': {'instruments': {}, 'name': 'Morley 2012', 'citation': 'Morley et al. (2012)', 'bibcode': '2012ApJ...756..172M', 'altname': ['morley','morley2012'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'fsed': 'f5'}}, \
    'morley14': {'instruments': {}, 'name': 'Morley 2014', 'citation': 'Morley et al. (2014)', 'bibcode': '2014ApJ...787...78M', 'altname': ['morley2014'], 'default': {'teff': 300., 'logg': 5.0, 'z': 0., 'fsed': 'f5', 'cld': 'h50'}}, \
    'nextgen99': {'instruments': {}, 'name': 'Phoenix NextGen', 'citation': 'Hauschildt et al. (1999)', 'bibcode': '1999ApJ...525..871H', 'altname': ['nextgen,hauschildt,hauschildt99,hauschildt1999'], 'default': {'teff': 2000., 'logg': 5.0, 'z': 0.0}}, \
    'saumon12': {'instruments': {}, 'name': 'Saumon 2012', 'citation': 'Saumon et al. (2012)', 'bibcode': '2012ApJ...750...74S', 'altname': ['saumon','saumon2012'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.}}, \
#    'btcond': {'instruments': {}, 'name': 'BT Cond', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['dusty-cond','bt-cond'], 'rawfolder': '/Volumes/splat/models/btcond/ORIGINAL/', 'default': {'teff': 1500., 'logg': 5.0, 'z': 0.0, 'enrich': 0.0}}, \
#    'btdusty': {'instruments': {}, 'name': 'BT Dusty', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['dusty-bt','bt-dusty'], 'rawfolder': '/Volumes/splat/models/btdusty/ORIGINAL/', 'default': {'teff': 1500., 'logg': 5.0, 'z': 0.0}}, \
}

SPECTRAL_MODEL_PARAMETERS_INORDER = ['teff','logg','z','fsed','cld','kzz','ad','enrich','carbon','oxygen','broad','logpmin','logpmax']

SPECTRAL_MODEL_PARAMETERS = {\
    'teff': {'name': 'temperature', 'prefix': 't', 'unit': u.K, 'default': 1000.0, 'title': '$T_{eff}$ (K)', 'type': 'continuous'}, \
    'logg': {'name': 'gravity', 'prefix': 'g', 'unit': u.dex, 'default': 5.0, 'title': '$\log{g}$ (cgs)', 'type': 'continuous'}, \
    'z': {'name': 'metallicity', 'prefix': 'z', 'unit': u.dex, 'default': 0., 'title': '$[M/H]$', 'type': 'continuous'}, \
    'fsed': {'name': 'rainout', 'prefix': 'f', 'unit': u.m/u.m, 'default': 'nc', 'title': '$f_{sed}$', 'type': 'discrete'}, \
    'cld': {'name': 'cloud', 'prefix': 'c', 'unit': u.m/u.m, 'default': 'nc', 'title': 'Cloud or Condensation Treatment', 'type': 'discrete'}, \
    'kzz': {'name': 'mixing', 'prefix': 'k', 'unit': u.m/u.m, 'default': 'eq', 'title': '$log\ \kappa_{zz}$ (cgs)', 'type': 'discrete'},\
    'ad': {'name': 'adiabat', 'prefix': 'ad', 'unit': u.m/u.m, 'default': 1., 'title': 'Adiabatic Index', 'type': 'continuous'},\
    'enrich': {'name': 'alpha enrichment', 'prefix': 'en', 'unit': u.dex, 'default': 0., 'title': 'Alpha Element Enrichment', 'type': 'continuous'},\
    'carbon': {'name': 'carbon enrichment', 'prefix': 'ca', 'unit': u.dex, 'default': 0., 'title': 'Carbon Enrichment', 'type': 'continuous'},\
    'oxygen': {'name': 'oxygen enrichment', 'prefix': 'ox', 'unit': u.dex, 'default': 0., 'title': 'Oxygen Enrichment', 'type': 'continuous'},\
    'broad': {'name': 'broadening', 'prefix': 'br', 'unit': u.m/u.m, 'default': 'A', 'title': 'Alkali Line Broadening Prescription', 'type': 'discrete'},\
    'logpmin': {'name': 'log pressure top', 'prefix': 'pt', 'unit': u.dex, 'default': -8., 'title': 'log Minimum Pressure (bar)', 'type': 'continuous'},\
    'logpmax': {'name': 'log pressure bottom', 'prefix': 'pb', 'unit': u.dex, 'default': 4., 'title': 'log Maximum Pressure (bar)', 'type': 'continuous'},\
    'radius': {'name': 'radius', 'prefix': 'r', 'unit': u.Rsun, 'default': 0., 'title': 'Radius (R$_{\odot}$)', 'type': 'continuous'},\
    }



def checkLocal(inputfile):
    '''
    :Purpose: Checks if a file is present locally or within the SPLAT
                code directory
    :Example:
       >>> import splat
       >>> spl.checkLocal('spl.py')
       True  # found the code
       >>> spl.checkLocal('parameters.txt')
       False  # can't find this file
       >>> spl.checkLocal('SpectralModels/BTSettl08/parameters.txt')
       True  # found it
    '''
    #print(inputfile)
    #sys.exit()
    if not os.path.exists(os.path.normpath(inputfile)):
        return ''
    else:
        return inputfile



def isUnit(s):
    '''
    :Purpose: 
        Checks if something is an astropy unit quantity; written in response to the 
        many ways that astropy now codes unit quantities

    :Required Inputs: 
        :param s: quantity to be checked

    :Optional Inputs: 
        None

    :Output: 
        True or False

    :Example:
    >>> import splat
    >>> import astropy.units as u
    >>> print splat.isUnit(3)
        False
    >>> print splat.isUnit(3.*u.s)
        True
    >>> print splat.isUnit(3.*u.s/u.s)
        True
    >>> print splat.isUnit((3.*u.s/u.s).value)
        False
    '''
    return isinstance(s,u.quantity.Quantity) or \
        isinstance(s,u.core.Unit) or \
        isinstance(s,u.core.CompositeUnit) or \
        isinstance(s,u.core.IrreducibleUnit) or \
        isinstance(s,u.core.NamedUnit) or \
        isinstance(s,u.core.PrefixUnit)





def loadModelParameters(*args,**kwargs):
    '''
    Purpose: 
        Assistant routine for `loadModel()`_ that loads in the spectral model grid points.
    .. _`loadModel()` : api.html#splat_model.loadModel
    Required Inputs:
        :param: model: set of models to use; see options in `loadModel()`_ (default = 'BTSettl2008')
    
    Optional Inputs:
        **instrument = 'RAW'**: search specifically for an instrument-designed model 
        **pandas = False**: return a pandas Dataframe
    The parameters for `loadModel()`_ can also be used here.
    Output:
        A dictionary containing the individual parameter values for the grid points in the given model set (not all of these grid points will be filled); this dictionary includs a list of dictionaries containing the individual parameter sets.
    '''

# model set
    modelset = False
    if len(args) > 0: modelset = args[0]
    modelset = kwargs.get('modelset',modelset)
    modelset = kwargs.get('model',modelset)
    modelset = kwargs.get('set',modelset)
    #mset = checkSpectralModelName(modelset)
    #if mset == False:
    #    raise NameError('\nInput model set {} not in defined set of models:\n{}\n'.format(modelset,list(SPECTRAL_MODELS.keys())))
    mset = 'btsettl08'
    """
# instrument
    instrument = ''
    if len(args) > 1: instrument = args[1]
    instrument = kwargs.get('instrument',instrument)
    instrument = kwargs.get('instr',instrument)
    instr = checkInstrument(instrument)
    if instr != False: instrument = instr
    if instrument not in list(SPECTRAL_MODELS[mset]['instruments'].keys()):
        ins = list(SPECTRAL_MODELS[mset]['instruments'].keys())
        if 'ORIGINAL' in ins: ins.remove('ORIGINAL')
        if len(ins) == 0: raise ValueError('\nNo SPLAT-processed models for {}'.format(mset))
        instrument = ins[0]
    """
    instrument = 'APOGEE-RAW'
    instr      = ''
# folder for models
    #Get the path of apogee_tools file
    FULL_PATH  = os.path.realpath(__file__)
    BASE, NAME = os.path.split(FULL_PATH)
    ModelPath = ap.LIBRARIES + '/%s/%s/'%(ap.model["grid_name"], ap.data["instrument"])
    mfolder   = os.path.normpath(ModelPath)
    #print(mfolder)
    #mfolder = os.path.normpath(SPECTRAL_MODELS[mset]['instruments'][instrument])
    if not os.access(mfolder, os.R_OK):
#        raise NameError('\nInstrument setting {} is not defined for model set {}\n'.format(instrument,mset))
#        mfolder = os.path.normpath(SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/'+mset)
#        if not os.access(mfolder, os.R_OK):
        raise OSError('\nCould not find model folder {}\n'.format(mfolder))

    parameters = {'model': mset, 'instrument': instrument, 'parameter_sets': []}
    #for ms in list(SPECTRAL_MODELS[mset]['default'].keys()):
    #    parameters[ms] = []
    parameters['teff'] = [3000.]
    parameters['logg'] = [5.0]
    parameters['z']    = [0.]
#    print(parameters.keys())

# establish parameters from list of filenames
#    if kwargs.get('old',False) == False:
    mfiles = glob.glob(mfolder+'/*.txt')
    if instr == 'RAW' or len(mfiles) == 0: mfiles = glob.glob(mfolder+'/*.gz')
    if len(mfiles) == 0:
        raise ValueError('\nCould not find any model files in {}'.format(mfolder))
    for mf in mfiles:
        p = {'model': mset, 'instrument': instrument}
        sp = numpy.array(os.path.basename(mf).replace('.txt','').replace('.gz','').replace(mset+'_','').replace('_'+instrument,'').split('_'))
        if '' in sp: 
            sp = list(sp)
            sp.remove('')
            sp = numpy.array(sp)
        for ms in SPECTRAL_MODEL_PARAMETERS_INORDER:
#
            if ms in list(parameters.keys()):
#
                val = sp[[SPECTRAL_MODEL_PARAMETERS[ms]['prefix'] in l for l in sp]][0][len(SPECTRAL_MODEL_PARAMETERS[ms]['prefix']):]

                if SPECTRAL_MODEL_PARAMETERS[ms]['type'] == 'continuous': 
                    val = float(val)
                parameters[ms].append(val)
                p[ms] = val
        parameters['parameter_sets'].append(p)

#        if len(sp) >= len(SPECTRAL_MODEL_PARAMETERS_INORDER):
#            p = {'model': mset}
#            for i,ms in enumerate(SPECTRAL_MODEL_PARAMETERS_INORDER):
#                if sp[i] not in parameters[ms]:
#                    parameters[ms].append(sp[i])
#                p[ms] = sp[i]
#            parameters['parameter_sets'].append(p)
    for ms in SPECTRAL_MODEL_PARAMETERS_INORDER:
#        if ms=='teff' or ms =='logg' or ms=='z':
#            parameters[ms] = [float(x) for x in parameters[ms]]
        if ms in list(parameters.keys()):
            parameters[ms].sort()
            parameters[ms] = numpy.array(parameters[ms])

#    print(parameters.keys(),parameters)

    if kwargs.get('pandas',False) == True:
        dp = pandas.DataFrame(parameters['parameter_sets'])
 #       for ms in ['teff','logg','z']: dp[ms] = [float(t) for t in dp[ms]]
        return dp
    else:
        return parameters




def loadModel(modelset='btsettl08',instrument='APOGEE-RAW',raw=False,sed=False,*args,**kwargs):
    '''
    Purpose: 
        Loads up a model spectrum based on a set of input parameters. The models may be any one of the following listed below. For parameters between the model grid points, loadModel calls the function `_loadInterpolatedModel()`_.

    .. _`_loadInterpolatedModel()` : api.html#splat_model._loadInterpolatedModel

    Required Inputs:
        :param: **model**: The model set to use; may be one of the following:

            - *nextgen99*: model set from `Allard et al. (1999) <http://adsabs.harvard.edu/abs/2012RSPTA.370.2765A>`_  with effective temperatures of 900 to 1600 K (steps of 100 K); surface gravities of 5.0 and 5.5 in units of cm/s^2; and metallicity fixed to solar (alternate designations: `nextgen`)
            - *cond01*: model set from `Allard et al. (2001) <http://adsabs.harvard.edu/abs/2001ApJ...556..357A>`_  with effective temperatures of 100 to 4000 K (steps of 100 K); surface gravities of 4.0 to 6.0 in units of cm/s^2 (steps of 0.5 dex); and metallicity fixed to solar; with condensate species removed from the photosphere (alternate designation: `cond`)
            - *dusty01*: model set from `Allard et al. (2001) <http://adsabs.harvard.edu/abs/2001ApJ...556..357A>`_  with effective temperatures of 500 to 3000 K (steps of 100 K); surface gravities of 3.5 to 6.0 in units of cm/s^2 (steps of 0.5 dex); and metallicity fixed to solar; with condensate species left in situ (alternate designation: `dusty`)
            - *burrows06*: model set from `Burrows et al. (2006) <http://adsabs.harvard.edu/abs/2006ApJ...640.1063B>`_ with effective temperatures of 700 to 2000 K (steps of 50 K); surface gravities of 4.5 to 5.5 in units of cm/s^2 (steps of 0.1 dex); metallicity of -0.5, 0.0 and 0.5; and either no clouds or grain size 100 microns (fsed = 'nc' or 'f100'). equilibrium chemistry is assumed. Note that this grid is not completely filled and some gaps have been interpolated (alternate designations: `burrows`, `burrows2006`)
            - *btsettl08*: (default) model set from `Allard et al. (2012) <http://adsabs.harvard.edu/abs/2012RSPTA.370.2765A>`_  with effective temperatures of 400 to 2900 K (steps of 100 K); surface gravities of 3.5 to 5.5 in units of cm/s^2 (steps of 0.5 dex); and metallicity of -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.3, and 0.5 for temperatures greater than 2000 K only; cloud opacity is fixed in this model, and equilibrium chemistry is assumed. Note that this grid is not completely filled and some gaps have been interpolated (alternate designations: `btsettled`, `btsettl`, `allard`, `allard12`)
            - *btsettl15*: model set from `Allard et al. (2015) <http://adsabs.harvard.edu/abs/2015A&A...577A..42B>`_  with effective temperatures of 1200 to 6300 K (steps of 100 K); surface gravities of 2.5 to 5.5 in units of cm/s^2 (steps of 0.5 dex); and metallicity fixed to solar (alternate designations: 'allard15','allard2015','btsettl015','btsettl2015','BTSettl2015')
            - *morley12*: model set from `Morley et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...756..172M>`_ with effective temperatures of 400 to 1300 K (steps of 50 K); surface gravities of 4.0 to 5.5 in units of cm/s^2 (steps of 0.5 dex); and sedimentation efficiency (fsed) of 2, 3, 4 or 5; metallicity is fixed to solar, equilibrium chemistry is assumed, and there are no clouds associated with this model (alternate designations: `morley2012`)
            - *morley14*: model set from `Morley et al. (2014) <http://adsabs.harvard.edu/abs/2014ApJ...787...78M>`_ with effective temperatures of 200 to 450 K (steps of 25 K) and surface gravities of 3.0 to 5.0 in units of cm/s^2 (steps of 0.5 dex); metallicity is fixed to solar, equilibrium chemistry is assumed, sedimentation efficiency is fixed at fsed = 5, and cloud coverage fixed at 50% (alternate designations: `morley2014`)
            - *saumon12*: model set from `Saumon et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...750...74S>`_ with effective temperatures of 400 to 1500 K (steps of 50 K); and surface gravities of 3.0 to 5.5 in units of cm/s^2 (steps of 0.5 dex); metallicity is fixed to solar, equilibrium chemistry is assumed, and no clouds are associated with these models (alternate designations: `saumon`, `saumon2012`)
            - *drift*: model set from `Witte et al. (2011) <http://adsabs.harvard.edu/abs/2011A%26A...529A..44W>`_ with effective temperatures of 1700 to 3000 K (steps of 50 K); surface gravities of 5.0 and 5.5 in units of cm/s^2; and metallicities of -3.0 to 0.0 (in steps of 0.5 dex); cloud opacity is fixed in this model, equilibrium chemistry is assumed (alternate designations: `witte`, `witte2011`, `helling`)
            - *madhusudhan11*: model set from `Madhusudhan et al. (2011) <http://adsabs.harvard.edu/abs/2011ApJ...737...34M>`_ with effective temperatures of 600 K to 1700 K (steps of 50-100 K); surface gravities of 3.5 and 5.0 in units of cm/s^2; and metallicities of 0.0 to 1.0 (in steps of 0.5 dex); there are multiple cloud prescriptions for this model, equilibrium chemistry is assumed (alternate designations: `madhusudhan`)
        
    Optional Inputs:
        :param: **teff**: effective temperature of the model in K (e.g. `teff` = 1000)
        :param: **logg**: log10 of the surface gravity of the model in cm/s^2 units (e.g. `logg` = 5.0)
        :param: **z**: log10 of metallicity of the model relative to solar metallicity (e.g. `z` = -0.5)
        :param: **fsed**: sedimentation efficiency of the model (e.g. `fsed` = 'f2')
        :param: **cld**: cloud shape function of the model (e.g. `cld` = 'f50')
        :param: **kzz**: vertical eddy diffusion coefficient of the model (e.g. `kzz` = 2)
        :param: **slit**: slit weight of the model in arcseconds (e.g. `slit` = 0.3)
        :param: **sed**: if set to True, returns a broad-band spectrum spanning 0.3-30 micron (applies only for BTSettl2008 models with Teff < 2000 K)

        :param: **folder**: string of the folder name containing the model set (default = '')
        :param: **filename**: string of the filename of the desired model; should be a space-delimited file containing columns for wavelength (units of microns) and surface flux (F_lambda units of erg/cm^2/s/micron) (default = '')
        :param: **force**: force the filename to be exactly as specified
        :param: **fast**: set to True to do a fast interpolation if needed, only for Teff and logg (default = False)
        :param: **url**: string of the url to the SPLAT website (default = 'http://www.browndwarfs.org/splat/')

    Output:
        A SPLAT Spectrum object of the interpolated model with wavelength in microns and surface fluxes in F_lambda units of erg/cm^2/s/micron.

    Example:

    >>> import splat
    >>> mdl = splat.loadModel(teff=1000,logg=5.0)
    >>> mdl.info()
        BTSettl2008 model with the following parmeters:
        Teff = 1000 K
        logg = 5.0 cm/s2
        z = 0.0
        fsed = nc
        cld = nc
        kzz = eq
        Smoothed to slit width 0.5 arcseconds
    >>> mdl = splat.loadModel(teff=2500,logg=5.0,model='burrows')
        Input value for teff = 2500 out of range for model set burrows06
        Warning: Creating an empty Spectrum object
    '''


# path to model and set local/online
# by default assume models come from local SPLAT directory

#   REMOVED 10/19/2017
#    local = kwargs.get('local',True)
#    online = kwargs.get('online',not local and not checkOnline())
#    local = not online
#    kwargs['local'] = local
#    kwargs['online'] = online
#    kwargs['url']  = kwargs.get('url',SPLAT_URL)
    kwargs['ismodel'] = True
    kwargs['force']   = kwargs.get('force',False)
    kwargs['folder']  = kwargs.get('folder','./')
    runfast           = kwargs.get('runfast',False)
    verbose           = kwargs.get('verbose',False)

    FULL_PATH  = os.path.realpath(__file__)
    BASE, NAME = os.path.split(FULL_PATH)
    ModelPath = ap.LIBRARIES + '/%s/%s/'%(ap.model["grid_name"], ap.data["instrument"])


# has a filename been passed? check first if it is a model set name
# otherwise assume this file is a local file
# and check that the path is correct if its fully provided
# otherwise assume path is inside provided folder keyword
    if len(args) > 0:
        #mset = checkSpectralModelName(args[0])
        mset = 'btsettl08'
        if mset != False: modelset=mset
        else:
            kwargs['filename'] = os.path.normpath(args[0])
            if not os.path.exists(kwargs['filename']):
                kwargs['filename'] = os.path.normpath(kwargs['folder']+os.path.basename(kwargs['filename']))
                if not os.path.exists(kwargs['filename']):
                    raise NameError('\nCould not find model file {} or {}'.format(kwargs['filename'],kwargs['folder']+os.path.basename(kwargs['filename'])))

# check if already read in
            if kwargs['filename'] in list(MODELS_READIN.keys()) and runfast == True:
#            if verbose: print('RUNFAST 1: {}'.format(kwargs['filename']))
                return MODELS_READIN[kwargs['filename']]
            else:
                MODELS_READIN[kwargs['filename']] = Spectrum(**kwargs)
                return MODELS_READIN[kwargs['filename']]


# set up the model set
#    modelset = kwargs.get('model','BTSettl2008')
    modelset = kwargs.get('model',modelset)
    modelset = kwargs.get('set',modelset)
    #mset = checkSpectralModelName(modelset)
    mset = 'btsettl08'
    if mset == False: raise ValueError('Could not find model set {}; possible options are {}'.format(modelset,list(SPECTRAL_MODELS.keys())))
    kwargs['modelset'] = mset

#    kwargs['instrument'] = kwargs.get('instrument','SPEX-PRISM')
    instrument = kwargs.get('instr',instrument)
    if raw == True: instrument = 'RAW'
    if sed == True: instrument = 'SED'
    #inst = checkInstrument(instrument)
    inst = 'APOGEE-RAW'
    if inst != False: instrument = inst
    if instrument not in list(SPECTRAL_MODELS[mset]['instruments']):
        raise ValueError('Models for set {} and instrument {} have not yet been computed; run processModelsToInstrument()'.format(kwargs['modelset'],instrument))
    kwargs['instrument'] = instrument
    kwargs['name'] = kwargs['modelset']+' ('+kwargs['instrument']+')'


# check that model data is available
#    kwargs['folder'] = kwargs.get('folder',os.path.normpath(SPECTRAL_MODELS[kwargs['model']]['folder']+'/'+kwargs['instrument']+'/'))

    kwargs['folder'] = ModelPath
    #kwargs['folder'] = os.path.normpath(SPECTRAL_MODELS[kwargs['modelset']]['instruments'][kwargs['instrument']])
    if not os.path.exists(kwargs['folder']):
        finit = kwargs['folder']
        kwargs['folder'] = os.path.normpath(SPLAT_PATH+SPECTRAL_MODEL_FOLDER+kwargs['modelset']+'/'+kwargs['instrument']+'/')
        if not os.path.exists(kwargs['folder']):
            raise ValueError('\nCould not locate folder {} or {} for model {} and instrument {}; make sure models are properly located'.format(finit,kwargs['folder'],kwargs['modelset'],kwargs['instrument']))

# preset defaults
    mparam = {}
    for ms in SPECTRAL_MODEL_PARAMETERS_INORDER:
        if ms in list(SPECTRAL_MODELS[kwargs['modelset']]['default'].keys()):
            mparam[ms] = kwargs.get(ms,SPECTRAL_MODELS[kwargs['modelset']]['default'][ms])
            if isUnit(mparam[ms]):
                mparam[ms] = (mparam[ms].to(SPECTRAL_MODEL_PARAMETERS[ms]['unit'])).value
    if len(mparam.keys()) == 0:
        raise ValueError('\nDid not have any parameters to set; this is a error in the program')
    for ms in mparam.keys(): kwargs[ms] = mparam[ms]

# generate model filename
    
    filename = os.path.join(kwargs['folder'],kwargs['modelset'])
    #filename = os.path.join(SPECTRAL_MODELS[kwargs['modelset']]['instruments'][kwargs['instrument']],kwargs['modelset'])

    for k in SPECTRAL_MODEL_PARAMETERS_INORDER:
        if k in list(SPECTRAL_MODELS[kwargs['modelset']]['default'].keys()):
            if k in list(mparam.keys()): val = mparam[k] 
            else: val = SPECTRAL_MODELS[mset]['default'][k]
            if SPECTRAL_MODEL_PARAMETERS[k]['type'] == 'continuous':
                kstr = '_{}{:.2f}'.format(SPECTRAL_MODEL_PARAMETERS[k]['prefix'],float(val))
            else:
                kstr = '_{}{}'.format(SPECTRAL_MODEL_PARAMETERS[k]['prefix'],val)
            if k == 'teff': kstr = '_{}{}'.format(SPECTRAL_MODEL_PARAMETERS[k]['prefix'],int(val))
            elif k == 'z': kstr = '_{}{:.2f}'.format(SPECTRAL_MODEL_PARAMETERS[k]['prefix'],float(val)-0.0001)
            filename=filename+kstr
    kwargs['filename'] = filename+'_{}.txt'.format(kwargs['instrument'])

#    kwargs['filename'] = os.path.normpath(kwargs['folder'])+'{}_{:.0f}_{:.1f}_{:.1f}_{}_{}_{}_{}.txt'.\
#        format(kwargs['model'],float(kwargs['teff']),float(kwargs['logg']),float(kwargs['z'])-0.001,kwargs['fsed'],kwargs['cld'],kwargs['kzz'],kwargs['instrument']))

#    if kwargs.get('sed',False):
#        kwargs['filename'] = kwargs['folder']+kwargs['model']+'_{:.0f}_{:.1f}_{:.1f}_nc_nc_eq_sed.txt'.\
#            format(float(kwargs['teff']),float(kwargs['logg']),float(kwargs['z'])-0.001)

# get model parameters
#        parameters = loadModelParameters(**kwargs)
#        kwargs['path'] = kwargs.get('path',parameters['path'])
# check that given parameters are in range
#        for ms in MODEL_PARAMETER_NAMES[0:3]:
#            if (float(kwargs[ms]) < parameters[ms][0] or float(kwargs[ms]) > parameters[ms][1]):
#                raise NameError('\n\nInput value for {} = {} out of range for model set {}\n'.format(ms,kwargs[ms],kwargs['set']))
#        for ms in MODEL_PARAMETER_NAMES[3:6]:
#            if (kwargs[ms] not in parameters[ms]):
#                raise NameError('\n\nInput value for {} = {} not one of the options for model set {}\n'.format(ms,kwargs[ms],kwargs['set']))


# have we already read in? if so just return saved spectrum object
    if kwargs['filename'] in list(MODELS_READIN.keys()):
#        if verbose: print('RUNFAST 2: {}'.format(kwargs['filename']))
        return MODELS_READIN[kwargs['filename']]


# check that folder/set is present either locally or online
# if not present locally but present online, switch to this mode
# if not present at either raise error

# REMOVED THIS 8/30/2017

#    folder = checkLocal(kwargs['folder'])
#    if folder=='':
#        folder = checkOnline(kwargs['folder'])
#        if folder=='':
#            print('\nCould not find '+kwargs['folder']+' locally or on SPLAT website')
#            print('\nAvailable model set options are:')
#            for s in DEFINED_MODEL_SET:
#                print('\t{}'.format(s))
#            raise NameError()
#        else:
#            kwargs['folder'] = folder
#            kwargs['local'] = False
#            kwargs['online'] = True
#    else:
#        kwargs['folder'] = folder

# check if file is present; if so, read it in, otherwise go to interpolated
# locally:
#    if kwargs.get('local',True) == True:
    file = checkLocal(kwargs['filename'])
    if file=='':
        file = checkLocal(kwargs['filename']+'.gz')
        if file=='':
            if kwargs['force']: raise NameError('\nCould not find '+kwargs['filename']+' locally\n\n')
            else: sp = _loadInterpolatedModel(**kwargs)
        else: kwargs['filename'] = kwargs['filename']+'.gz'
#                kwargs['local']=False
#                kwargs['online']=True
#        else:
#    else:
    if file != '':
        #import splat
        sp = ap.Spectrum(**kwargs)
        #print(kwargs)
        #print(dir(sp))
        #print(MODELS_READIN[kwargs['filename']])
        #MODELS_READIN[kwargs['filename']] = sp

# populate model parameters
    setattr(sp,'modelset',kwargs['modelset'])
    #print(sp.modelset)
    #sys.exit()
    setattr(sp,'instrument',kwargs['instrument'])
    for k in list(SPECTRAL_MODELS[kwargs['modelset']]['default'].keys()):
        if k in list(mparam.keys()): setattr(sp,k,mparam[k])
        else: setattr(sp,k,SPECTRAL_MODELS[mset]['default'][k])

# add to read in files
    
    return sp



################################################################################



def _loadModelParameters(*args,**kwargs):
    '''
    Purpose: 
        Assistant routine for `loadModel()`_ that loads in the spectral model grid points.

    .. _`loadModel()` : api.html#splat_model.loadModel

    Required Inputs:
        :param: model: set of models to use; see options in `loadModel()`_ (default = 'BTSettl2008')
    
    Optional Inputs:
        **instrument = 'RAW'**: search specifically for an instrument-designed model 
        **pandas = False**: return a pandas Dataframe


    The parameters for `loadModel()`_ can also be used here.

    Output:
        A dictionary containing the individual parameter values for the grid points in the given model set (not all of these grid points will be filled); this dictionary includs a list of dictionaries containing the individual parameter sets.
    '''

# model set
    modelset = False
    if len(args) > 0: modelset = args[0]
    modelset = kwargs.get('modelset',modelset)
    modelset = kwargs.get('model',modelset)
    modelset = kwargs.get('set',modelset)
    mset = 'btsettl08'#checkSpectralModelName(modelset)
    if mset == False:
        raise NameError('\nInput model set {} not in defined set of models:\n{}\n'.format(modelset,list(SPECTRAL_MODELS.keys())))

# instrument
    instrument = ''
    if len(args) > 1: instrument = args[1]
    instrument = kwargs.get('instrument',instrument)
    instrument = kwargs.get('instr',instrument)
    instr = 'APOGEE-RAW'#checkInstrument(instrument)
    if instr != False: instrument = instr
    #if instrument not in list(SPECTRAL_MODELS[mset]['instruments'].keys()):
    #    ins = list(SPECTRAL_MODELS[mset]['instruments'].keys())
    #    if 'ORIGINAL' in ins: ins.remove('ORIGINAL')
    #    if len(ins) == 0: raise ValueError('\nNo SPLAT-processed models for {}'.format(mset))
    #    instrument = ins[0]

# folder for models
    FULL_PATH  = os.path.realpath(__file__)
    BASE, NAME = os.path.split(FULL_PATH)
    #ModelPath = BASE + '/apogee_tools/libraries/btsettl08_highres/'
    ModelPath = ap.LIBRARIES + '/%s/%s/'%(ap.model["grid_name"], ap.data["instrument"])

    #mfolder = os.path.normpath(SPECTRAL_MODELS[mset]['instruments'][instrument])
    mfolder = ModelPath
    if not os.access(mfolder, os.R_OK):
#        raise NameError('\nInstrument setting {} is not defined for model set {}\n'.format(instrument,mset))
#        mfolder = os.path.normpath(SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/'+mset)
#        if not os.access(mfolder, os.R_OK):
        raise OSError('\nCould not find model folder {}\n'.format(mfolder))

    parameters = {'model': mset, 'instrument': instrument, 'parameter_sets': []}
    for ms in list(SPECTRAL_MODELS[mset]['default'].keys()):
        parameters[ms] = []
#    print(parameters.keys())

# establish parameters from list of filenames
#    if kwargs.get('old',False) == False:
    mfiles = glob.glob(mfolder+'/*.txt')
    if instr == 'RAW' or len(mfiles) == 0: mfiles = glob.glob(mfolder+'/*.gz')
    if len(mfiles) == 0:
        raise ValueError('\nCould not find any model files in {}'.format(mfolder))
    for mf in mfiles:
        p = {'model': mset, 'instrument': instrument}
        sp = numpy.array(os.path.basename(mf).replace('.txt','').replace('.gz','').replace(mset+'_','').replace('_'+instrument,'').split('_'))
        if '' in sp: 
            sp = list(sp)
            sp.remove('')
            sp = numpy.array(sp)
        for ms in SPECTRAL_MODEL_PARAMETERS_INORDER:
            if ms in list(parameters.keys()):
#                print(mf,sp,ms)
                val = sp[[SPECTRAL_MODEL_PARAMETERS[ms]['prefix'] in l for l in sp]][0][len(SPECTRAL_MODEL_PARAMETERS[ms]['prefix']):]
                if SPECTRAL_MODEL_PARAMETERS[ms]['type'] == 'continuous': val = float(val)
                parameters[ms].append(val)
                p[ms] = val
        parameters['parameter_sets'].append(p)

#        if len(sp) >= len(SPECTRAL_MODEL_PARAMETERS_INORDER):
#            p = {'model': mset}
#            for i,ms in enumerate(SPECTRAL_MODEL_PARAMETERS_INORDER):
#                if sp[i] not in parameters[ms]:
#                    parameters[ms].append(sp[i])
#                p[ms] = sp[i]
#            parameters['parameter_sets'].append(p)
    for ms in SPECTRAL_MODEL_PARAMETERS_INORDER:
#        if ms=='teff' or ms =='logg' or ms=='z':
#            parameters[ms] = [float(x) for x in parameters[ms]]
        if ms in list(parameters.keys()):
            parameters[ms].sort()
            parameters[ms] = numpy.array(parameters[ms])

#    print(parameters.keys(),parameters)

    if kwargs.get('pandas',False) == True:
        dp = pandas.DataFrame(parameters['parameter_sets'])
 #       for ms in ['teff','logg','z']: dp[ms] = [float(t) for t in dp[ms]]
        return dp
    else:
        return parameters



################################################################################
################################################################################
################################################################################



def _checkModelParametersInRange(mparam):
# list of model parameters provided
    mp = list(mparam.keys())
    if 'model' not in mp:
        mparam['model'] = 'BTSettl2008'
    if 'instrument' not in mp:
        mparam['instrument'] = 'SPEX-PRISM'
    parameters = _loadModelParameters(**mparam)
    flag = True

# check that given parameters are in model ranges
    for ms in SPECTRAL_MODEL_PARAMETERS_INORDER:
        if SPECTRAL_MODEL_PARAMETERS[ms]['type'] == 'continuous':
#            ms=='teff' or ms=='logg' or ms=='z':
            if ms in mp:
                if (float(mparam[ms]) < numpy.min(parameters[ms]) or float(mparam[ms]) > numpy.max(parameters[ms])):
                    print('\n\nInput value for {} = {} out of range for model set {}\n'.format(ms,mparam[ms],mparam['model']))
                    flag = False
        else:
            if ms in mp:
                if (mparam[ms] not in parameters[ms]):
                    print('\n\nInput value for {} = {} not one of the options for model set {}\n'.format(ms,mparam[ms],mparam['model']))
                    flag = False
    return flag



################################################################################



def _loadInterpolatedModel(*args,**kwargs):
    '''
    Purpose: 
        Generates as spectral model with is interpolated between model parameter grid points. This routine is called by `loadModel()`_, or it can be called on its own.

    .. _`loadModel()` : api.html#splat_model.loadModel

    Required Inputs:
        :param model: set of models to use; see options in `loadModel()`_

    Optional Inputs:
        :param: The parameters for `loadModel()`_ can also be used here.

    Output:
        A SPLAT Spectrum object of the interpolated model with wavelength in microns and surfae fluxes in F_lambda units of erg/cm^2/s/micron.

    Example:

    >>> import splat.model as spmdl
    >>> mdl = spmdl.loadModel(teff=1000,logg=5.0)
    >>> mdl.info()
        BT-Settl (2008) Teff=1000 logg=5.0 [M/H]=0.0 atmosphere model with the following parmeters:
        Teff = 1000 K
        logg = 5.0 dex
        z = 0.0 dex
        fsed = nc
        cld = nc
        kzz = eq

        If you use this model, please cite Allard, F. et al. (2012, Philosophical Transactions of the Royal Society A, 370, 2765-2777)
        bibcode = 2012RSPTA.370.2765A
    '''
# attempt to generalize models to extra dimensions
    mkwargs = kwargs.copy()
    mkwargs['force'] = True
#    mkwargs['local'] = kwargs.get('local',True)

# set up the model set
#    mkwargs['model'] = mkwargs.get('model','BTSettl2008')
#    mkwargs['model'] = mkwargs.get('modelset',mkwargs['model'])
#    mkwargs['model'] = mkwargs.get('set',mkwargs['model'])
#    mkwargs['model'] = checkSpectralModelName(mkwargs['model'])
#    mkwargs['instrument'] = mkwargs.get('instrument','SPEX_PRISM')
#    mkwargs['instrument'] = checkInstrument(mkwargs['instrument'])
#    mkwargs['name'] = mkwargs['model']
    
#    mkwargs['folder'] = SPLAT_PATH+SPECTRAL_MODEL_FOLDER+mkwargs['model']+'/'


#    for ms in SPECTRAL_MODEL_PARAMETERS_INORDER:
#        mkwargs[ms] = kwargs.get(ms,SPECTRAL_MODEL_PARAMETERS[ms]['default'])

# some special defaults
#    if mkwargs['model'] == 'morley12':
#        if mkwargs['fsed'] == 'nc':
#            mkwargs['fsed'] = 'f2'
#    if mkwargs['model'] == 'morley14':
#        if mkwargs['fsed'] == 'nc':
#            mkwargs['fsed'] = 'f5'
#        if mkwargs['cld'] == 'nc':
#            mkwargs['cld'] = 'f50'

# first get model parameters
    if _checkModelParametersInRange(mkwargs) == False:
        raise ValueError('\n\nModel parameter values out of range for model set {}\n'.format(mkwargs['model']))
    
# check that given parameters are in range - RETHINK THIS
#    for ms in SPECTRAL_MODEL_PARAMETERS_INORDER:
#        if ms=='teff' or ms =='logg' or ms=='z':
#            if (float(mkwargs[ms]) < numpy.min(parameters[ms]) or float(mkwargs[ms]) > numpy.max(parameters[ms])):
#                raise ValueError('\n\nInput value for {} = {} out of range for model set {}\n'.format(ms,mkwargs[ms],mkwargs['model']))
#        else:
#            if (mkwargs[ms] not in parameters[ms]):
#                raise ValueError('\n\nInput value for {} = {} not one of the options for model set {}\n'.format(ms,mkwargs[ms],mkwargs['model']))



# FAST METHOD - just calculate a simple weight factor that linearly interpolates between grid points (all logarithmic)

    if kwargs.get('fast',True) == True:
        parameters = _loadModelParameters(mkwargs['model'],mkwargs['instrument'],pandas=True)
        mparams = {}
        mweights = {}
        mgrid = []
        pgrid = []
        plin = []
        for ms in list(SPECTRAL_MODEL_PARAMETERS.keys()):
            if ms in list(parameters.keys()):
                if SPECTRAL_MODEL_PARAMETERS[ms]['type'] == 'discrete': 
                    mparams[ms] = mkwargs[ms]
                    mweights[ms] = 1.
                    plin.append(ms)
                else:
                    l = parameters[parameters[ms] <= float(mkwargs[ms])].sort_values(ms)[ms].iloc[-1]
                    h = parameters[parameters[ms] >= float(mkwargs[ms])].sort_values(ms)[ms].iloc[0]
                    if ms == 'teff':
                        d = numpy.log10(h)-numpy.log10(l)
                        w = (numpy.log10(h)-numpy.log10(float(mkwargs[ms])))/d
                    else:
                        d = h-l
                        w = (h-float(mkwargs[ms]))/d
                    if d == 0.: w = 0.5
                    mparams[ms] = [l,h]
                    mweights[ms] = w
                    mgrid.append([l,h])
                    pgrid.append(ms)

# generate all possible combinations - doing this tediously due to concerns over ordering
        x = numpy.meshgrid(*mgrid)
        a = {}
        weights = numpy.ones(len(x[0].flatten()))
        for i,ms in enumerate(pgrid): 
            a[ms] = x[i].flatten()
            for j,v in enumerate(a[ms]):
                if v == mparams[ms][0]:
                    weights[j] *= mweights[ms]
                else:
                    weights[j] *= (1.-mweights[ms])

# read in models
        models = []
        for i in numpy.arange(len(weights)):
            mparam = copy.deepcopy(mkwargs)
            for ms in list(SPECTRAL_MODEL_PARAMETERS.keys()):
                if ms in list(parameters.keys()):
                    if SPECTRAL_MODEL_PARAMETERS[ms]['type'] == 'discrete': 
                        mparam[ms] = mkwargs[ms]
                    else:
                        mparam[ms] = a[ms][i]
            del mparam['filename']
            models.append(loadModel(**mparam))

# create interpolation
        #print('TEST2:', dir(models[0]))
        #print('TEST3:', models[0].flux)
        #print('TEST3:', models[0].flux.value)
        mflx = []
        for i,w in enumerate(models[0].wave):
            #val = numpy.array([numpy.log10(m.flux.value[i]) for m in models])
            val = numpy.array([numpy.log10(m.flux[i]) for m in models])
            mflx.append(10.**(numpy.sum(val*weights)/numpy.sum(weights)))


# REGULAR METHOD - uses meshgrid & griddata - about 4x slower
    else:

# identify grid points around input parameters
# 3x3 grid for teff, logg, z
        parameters = _loadModelParameters(mkwargs['model'],mkwargs['instrument'])

        tvals = numpy.array([float(p['teff']) for p in parameters['parameter_sets']])
        gvals = numpy.array([float(p['logg']) for p in parameters['parameter_sets']])
        zvals = numpy.array([float(p['z']) for p in parameters['parameter_sets']])
        tdiff = numpy.array([numpy.log10(float(mkwargs['teff']))-numpy.log10(v) for v in tvals])
        gdiff = numpy.array([float(mkwargs['logg'])-v for v in gvals])
        zdiff = numpy.array([float(mkwargs['z'])-v for v in zvals])
        dist = tdiff**2+gdiff**2+zdiff**2

# get closest models in 8 quadrant points
        mparams = []
#    mparam_names = []
        psets = numpy.array(parameters['parameter_sets'])
        for i in numpy.arange(0,2):
            dt = dist[numpy.where(tdiff*((-1)**i)>=0)]
            pt = psets[numpy.where(tdiff*((-1)**i)>=0)]
            gt = gdiff[numpy.where(tdiff*((-1)**i)>=0)]
            zt = numpy.round(zdiff[numpy.where(tdiff*((-1)**i)>=0)]*50.)/50.
            for j in numpy.arange(0,2):
                dg = dt[numpy.where(gt*((-1)**j)>=0)]
                pg = pt[numpy.where(gt*((-1)**j)>=0)]
                zg = zt[numpy.where(gt*((-1)**j)>=0)]
                for k in numpy.arange(0,2):
                    dz = dg[numpy.where(zg*((-1)**k)>=0)]
                    pz = pg[numpy.where(zg*((-1)**k)>=0)]

# if we can't get a quadrant point, quit out
                    if len(pz) == 0: 
#                    print(i,j,k)
#                    print(pg)
#                    print(zg)
#                    print(zg)
                        raise ValueError('\n\nModel parameter values out of range for model set {}\n'.format(mkwargs['model']))

                    pcorner = pz[numpy.argmin(dz)]
                    mparams.append(pz[numpy.argmin(dz)])

# generate meshgrid with slight offset and temperature on log scale
        rng = []
        for ms in SPECTRAL_MODEL_PARAMETERS_INORDER:
            if ms=='teff' or ms =='logg' or ms=='z':
                vals = [float(m[ms]) for m in mparams]
                r = [numpy.min(vals),numpy.max(vals)]
                if numpy.absolute(r[0]-r[1]) < 1.e-3*numpy.absolute(parameters[ms][1]-parameters[ms][0]):
                    r[1] = r[0]+1.e-3*numpy.absolute(parameters[ms][1]-parameters[ms][0])
                if ms == 'teff':
                    r = numpy.log10(r)
                rng.append(r)
        mx,my,mz = numpy.meshgrid(rng[0],rng[1],rng[2])

# read in unique models
        bmodels = {}
        models = []
        mp = copy.deepcopy(mparams[0])
        for i in numpy.arange(len(mx.flatten())):
            mp['teff'] = int(numpy.round(10.**(mx.flatten()[i])))
            mp['logg'] = my.flatten()[i]
            mp['z'] = mz.flatten()[i]
            mstr = '{:d}{:.1f}{:.1f}'.format(mp['teff'],mp['logg'],mp['z'])
            if mstr not in list(bmodels.keys()):
                bmodels[mstr] = loadModel(instrument=mkwargs['instrument'],force=True,**mp)
            models.append(bmodels[mstr])


#    mpsmall = [dict(y) for y in set(tuple(x.items()) for x in mparams)]
#    mpsmall = numpy.unique(numpy.array(mparams))
#    bmodels = []
#    bmodel_names = []
#    for m in mpsmall:
#        bmodels.append(loadModel(instrument=mkwargs['instrument'],force=True,**m))
#        mstr = ''
#        for ms in SPECTRAL_MODEL_PARAMETERS_INORDER: mstr+=str(m[ms])
#        bmodel_names.append(mstr)
#        if kwargs.get('verbose',False): print(m)
#    bmodels = numpy.array(bmodels)
#    bmodel_names = numpy.array(bmodel_names)

# now set up model array in mx,my,mz order
#    mparam_names = []
#    models = []
#    for i,m in enumerate(mparam_names):
#    for i in numpy.arange(len(mx.flatten())):
#        mstr = '{:d}{:.1f}{:.1f}'.format(mx.flatten()[i],my.flatten()[i],mz.flatten()[i])
#        for ms in SPECTRAL_MODEL_PARAMETERS_INORDER: mstr+=str(mparams[-1][ms])
#        models.append(bmodels[numpy.where(bmodel_names==m)][0])
#        mparam_names.append(mstr)

#    models = []
#    for i,m in enumerate(mparam_names):
#        models.append(bmodels[numpy.where(bmodel_names==m)][0])

        if kwargs.get('debug',False):
            print(mx.flatten())
            print([m.teff for m in models])
            print(my.flatten())
            print([m.logg for m in models])
            print(mz.flatten())
            print([m.z for m in models])
#        print(mparams)
#        print(mparam_names)


# final interpolation
        mflx = []
        for i,w in enumerate(models[0].wave):
            #val = numpy.array([numpy.log10(m.flux.value[i]) for m in models])
            val = numpy.array([numpy.log10(m.flux[i]) for m in models])
            mflx.append(10.**(griddata((mx.flatten(),my.flatten(),mz.flatten()),\
                val,(numpy.log10(float(mkwargs['teff'])),float(mkwargs['logg']),float(mkwargs['z'])),'linear')))

    return ap.Spectrum(wave=models[0].wave,flux=mflx*models[0].funit,**mkwargs)



#parameters = {'model': mset, 'instrument': instrument, 'parameter_sets': []}
##for ms in list(SPECTRAL_MODELS[mset]['default'].keys()):
##    parameters[ms] = []
#parameters['teff'] = [3000.]
#parameters['logg'] = [5.0]
#parameters['z']    = [0.]}


gridparam = loadModelParameters('btsettl08','APOGEE-RAW') 
#print(gridparam)
#mdl = loadModel()
mdl = loadModel(teff = 3050, logg = 4.1, z = 0)

plt.plot(mdl.wave, mdl.flux)
plt.show()

