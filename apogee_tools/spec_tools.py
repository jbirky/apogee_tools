import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import h5py

# from apogee_tools.core import *
import apogee_tools as ap
#from .core import *
# import apogee_tools
# from apogee_tools import Spectrum

#Get the path of apogee_tools file
FULL_PATH  = os.path.realpath(__file__)
TOOLS_PATH = os.path.split(FULL_PATH)[0]
BASE = os.path.split(TOOLS_PATH)[0]

AP_PATH = os.environ['APOGEE_DATA']


def HDF5Convert(filename, **kwargs):

    """
    Create an hdf5 file for a spectrum.
    @Jessica Birky

    Input:  filename : string: name of file to be written
            'wls'.   : wavelength array
            'fls'    : flux array
            'sigmas' : uncertainty array

    Output: filename.hdf5
    """
    
    with h5py.File(filename+'.hdf5', 'w') as hf:
        hf.create_dataset('wls', data=kwargs.get('wls'))
        hf.create_dataset('fls', data=kwargs.get('fls'))
        hf.create_dataset('sigmas', data=kwargs.get('sigmas'))
        if 'masks' in kwargs:
            hf.create_dataset('masks', data=kwargs.get('masks'))


class HDF5Interface(object):

    """
    class HDF5Interface(object) has been taken from the Starfish code by Ian Czekala
    see: https://github.com/iancze/Starfish
    Function: Connect to an HDF5 file that stores spectra.
    """

    def __init__(self, **kwargs):
        '''
        :param filename: the name of the HDF5 file
        :type param: string
        :param ranges: optionally select a smaller part of the grid to use.
        :type ranges: dict
        '''
        filename = kwargs.get('hdf5_path', BASE + '/libraries/BTSETTLb_APOGEE.hdf5')
        key_name = kwargs.get('key_name', "t{0:.0f}g{1:.1f}z{2:.1f}")

        self.filename = os.path.expandvars(filename)
        self.key_name = key_name

        # In order to properly interface with the HDF5 file, we need to learn
        # a few things about it

        # 1.) Which parameter combinations exist in the file (self.grid_points)
        # 2.) What are the minimum and maximum values for each parameter (self.bounds)
        # 3.) Which values exist for each parameter (self.points)

        with h5py.File(self.filename, "r") as hdf5:
            self.wl = hdf5["wl"][:]
            self.wl_header = dict(hdf5["wl"].attrs.items())
            self.dv = self.wl_header["dv"]
            self.grid_points = hdf5["pars"][:]

        #determine the bounding regions of the grid by sorting the grid_points
        low = np.min(self.grid_points, axis=0)
        high = np.max(self.grid_points, axis=0)
        self.bounds = np.vstack((low, high)).T
        self.points = [np.unique(self.grid_points[:, i]) for i in range(self.grid_points.shape[1])]

        self.ind = None #Overwritten by other methods using this as part of a ModelInterpolator

    def load_flux(self, parameters):
        '''
        Load just the flux from the grid, with possibly an index truncation.

        :param parameters: the stellar parameters
        :type parameters: np.array

        :raises KeyError: if spectrum is not found in the HDF5 file.

        :returns: flux array
        '''
        print(parameters)
        key = self.key_name.format(*parameters)
        with h5py.File(self.filename, "r") as hdf5:
            try:
                if self.ind is not None:
                    fl = hdf5['flux'][key][self.ind[0]:self.ind[1]]
                else:
                    fl = hdf5['flux'][key][:]
            except KeyError as e:
                raise C.GridError(e)

        #Note: will raise a KeyError if the file is not found.

        return fl

    @property
    def fluxes(self):
        '''
        Iterator to loop over all of the spectra stored in the grid, for PCA.

        Loops over parameters in the order specified by grid_points.
        '''
        for grid_point in self.grid_points:
            yield self.load_flux(grid_point)

    def load_flux_hdr(self, parameters):
        '''
        Just like load_flux, but also returns the header
        '''
        key = self.key_name.format(*parameters)
        with h5py.File(self.filename, "r") as hdf5:
            try:
                hdr = dict(hdf5['flux'][key].attrs)
                if self.ind is not None:
                    fl = hdf5['flux'][key][self.ind[0]:self.ind[1]]
                else:
                    fl = hdf5['flux'][key][:]
            except KeyError as e:
                raise C.GridError(e)

        #Note: will raise a KeyError if the file is not found.

        return (fl, hdr)


def loadGrid(**kwargs):

    """
    Load model grids from .hdf5 files (which were created by Starfish) stored on github. 
    @Jessica Birky

    Input:  'gridPath' : directory to .hdf5 models file

    Output: 'params' : array: [Teff, logg, [Fe/H]]
            'wave'   : wavelength 
            'flux'   : flux 
    """

    # import Starfish
    # from Starfish.grid_tools import HDF5Interface

    # optional key word arguments:
    grid_lib = BASE + '/libraries/BTSETTLb_APOGEE.hdf5'

    gridPath = kwargs.get('gridPath', grid_lib)
    myHDF5 = HDF5Interface(hdf5_path=gridPath)

    if 'params' in kwargs:
        params = kwargs.get('params')

    elif ('teff' in kwargs) and ('logg' in kwargs) and ('fe_h' in kwargs):
        teff = kwargs.get('teff')
        logg = kwargs.get('logg')
        fe_h = kwargs.get('fe_h')

        params = [teff, logg, fe_h]

    wave = myHDF5.wl
    flux = myHDF5.load_flux(np.array(params))

    return params, wave, flux


def _rvShift(wave, **kwargs): 

    """
    Shift wavelenth of spectrum object by specific rv.
    @Jessica Birky

    Input: 'rv' : radial velocity (km/s)

    """
    
    rv = kwargs.get('rv', -80) 

    shift   = 1. + rv/299792.
    rv_wave = wave*shift

    return rv_wave


def rvShiftSpec(sp, **kwargs): 

    """
    Rv shift, given spectrum object, return spectrum object
    @Jessica Birky

    Input: 'rv' : radial velocity (km/s)

    """
    
    rv = kwargs.get('rv', -80) 

    shift   = 1. + rv/299792.
    rv_wave = sp.wave*shift
    rv_spec = ap.Spectrum(wave=rv_wave, flux=sp.flux)

    return rv_spec


def calcScale(sp1, sp2):

    """
    Find the scale factor for the second spectrum to minimize chi-squared value
    @Jessica Birky

    Input:  two spectrum objects

    Output: compute scale factor
    """

    wave1, flux1 = sp1.wave, sp1.flux
    wave2, flux2 = sp2.wave, sp2.flux
    unc = sp1.sigmas

    a = sum(flux2*flux1/unc**2)
    b = sum(flux2**2/unc**2)

    scale = a/b

    return scale


def compareSpectra(sp1, sp2, **kwargs):

    """
    Compute the chi-square fit between two spectra (or model and spectrum)
    @Jessica Birky

    Input:  two spectrum objects

    Output: chi-squared value and scaled spectrum objects
    """

    fit_range = kwargs.get('fit_range', [sp1.wave[0], sp1.wave[-1]])
    fit_scale = kwargs.get('fit_scale', True)
    normalize = kwargs.get('norm', True)

    from scipy.interpolate import interp1d

    #sp2 will be interpolated to the same wavelenth values as sp1
    wave = sp1.wave
    unc = sp1.sigmas
    flux1 = sp1.flux

    #Interpolation function 
    f = interp1d(sp2.wave, sp2.flux, bounds_error=False, fill_value=0.)
    flux2 = f(sp1.wave)

    if normalize == True:
        #Make sure fluxes are normalized
        flux1 = flux1/max(flux1)
        flux2 = flux2/max(flux2)

    #Create a new spectrum object for sp2
    sp2 = ap.Spectrum(wave=sp2.wave, flux=np.array(flux2), params=sp2.params)

    if fit_scale == True:
        #Compute scale factor for 2nd spectrum
        scale = calcScale(sp1, sp2)
    else:
        scale = 1

    #Compute chi-squared value
    flux2 = scale*flux2
    diff = flux1 - flux2
    stat = diff**2/np.array(unc)**2
    chi  = sum(stat)

    sp2 = ap.Spectrum(wave=wave, flux=flux2, params=sp2.params)

    return chi, sp1, sp2


def subtractContinuum(spec, **kwargs):

    """
    Fit polynomial to spectrum and subtract from spectrum
    @Jessica Birky

    Input:  'deg'    : degree of polynomial
            'plot'   : True or False
            'xrange' : wavelegth range to fit continuum to

    Output: 'sub_spec'  : continuum normalized spectrum object
            'continuum' : continuum polynomial 
    """

    deg  = kwargs.get('deg', 20)
    plot = kwargs.get('plot', False)

    wave, flux = spec.wave, spec.flux

    #Cut spectrum to APOGEE wavelength range
    xrange = kwargs.get('xrange', [15200,16940])
    mask = np.where((wave > xrange[0]) & (wave < xrange[-1]))
    wave, flux = wave[mask], flux[mask]
    
    #Normalize flux
    flux = flux/max(flux)

    polynomial = np.poly1d(np.polyfit(wave, flux, deg))
    continuum  = polynomial(np.linspace(wave[0], wave[-1], len(wave)))

    sub_flux = flux - continuum + np.mean(flux)

    #Plot polynomial against spectrum
    if plot == True:
        plt.plot(wave, flux, color='k', alpha=.8)
        plt.plot(wave, continuum, color='c', alpha=.8)
        plt.plot(wave, sub_flux-.6, color='r')
        plt.show()
        plt.close()

    sub_spec = ap.Spectrum(wave=wave, flux=sub_flux, params=spec.params)

    return sub_spec, continuum


def readModels(**kwargs):

    """
    Input desired parameter ranges, return Spectrum objects
    Configured to read grid specified in config.yaml file inside apogee_tools
    """

    parrange = kwargs.get('parrange')
    gridPath = kwargs.get('gridPath')
    subCont  = kwargs.get('subCont', True)
    xrange   = kwargs.get('xrange', [15200,16940])

    teff_range = parrange[0]
    logg_range = parrange[1]
    fe_h_range = parrange[2]

    teffs = np.arange(teff_range[0], teff_range[1]+1, 100)
    loggs = np.arange(logg_range[0], logg_range[1]+.1, .5)
    fe_hs = np.arange(fe_h_range[0], fe_h_range[1]+.1, .5)

    specs = []

    grid_combos = [[i,j,k] for i in teffs for j in loggs for k in fe_hs]

    if subCont == True:
        for grid in grid_combos:
            params, wave, flux = loadGrid(params=grid, gridPath=gridPath)

            #Cut out parts of grid not in specified xrange
            mask = np.where((wave>xrange[0]) & (wave<xrange[1]))
            new_wave = wave[mask]
            new_flux = flux[mask]

            #Create spectrum object for the model spectrum
            sp = ap.Spectrum(wave=new_wave, flux=new_flux, params=grid)
            
            #Subtract continuum and store in new spectrum object
            norm_sp, cont = subtractContinuum(sp, plot=False)

            specs.append(norm_sp)

    else:
        for grid in grid_combos:
            params, wave, flux = loadGrid(params=grid, gridPath=gridPath)

            #Cut out parts of grid not in specified xrange
            mask = np.where((wave>xrange[0]) & (wave<xrange[1]))
            new_wave = wave[mask]
            new_flux = flux[mask]

            sp = ap.Spectrum(wave=new_wave, flux=new_flux, params=grid)
            specs.append(sp)

    return specs


def getModel(**kwargs):

    """
    Obtain just one model, given parameters. Make sure they are in order [Teff, logg, [Fe/H]]
    """

    grid     = kwargs.get('grid', 'BTSETTLb') 
    grid_lib = BASE + '/libraries/' + grid + '_APOGEE.hdf5'

    params   = kwargs.get('params', [3000, 5.0, 0.0])
    gridPath = kwargs.get('gridPath', grid_lib)
    xrange   = kwargs.get('xrange', [15200,16940])
    subCont  = kwargs.get('subCont', False)

    params, m_wave, m_flux = loadGrid(params=params, gridPath=gridPath)

    mask  = np.where((m_wave>xrange[0]) & (m_wave<xrange[1]))
    m_wave = m_wave[mask]
    m_flux = m_flux[mask]

    print(grid+': '+str(params))
    model_spec = ap.Spectrum(wave=m_wave, flux=m_flux, params=params, name=grid+': '+str(params))

    #Subtract continuum from the model
    if subCont == True:
        model_spec, cont = subtractContinuum(model_spec)

    return model_spec


def plotModel(spec, **kwargs):

    #Plot specific grid against spectrum

    #Retrieve grid library from the 'libraries' folder inside apogee_tools
    grid     = kwargs.get('grid', 'BTSETTLb')
    grid_lib = BASE + '/libraries/' + grid + '_APOGEE.hdf5'

    params   = kwargs.get('params')
    gridPath = kwargs.get('gridPath', grid_lib)
    save     = kwargs.get('save', False)
    output   = kwargs.get('output', spec.name+'.pdf')
    xrange   = kwargs.get('xrange', [spec.wave[0],spec.wave[-1]])

    # params, m_wave, m_flux = loadGrid(params=params, gridPath=gridPath)
    m_spec = getModel(params=params, xrange=xrange)
    m_wave, m_flux = m_spec.wave, m_spec.flux

    mask = np.where((spec.wave>xrange[0]) & (spec.wave<xrange[1]))
    wave = spec.wave[mask]
    flux = spec.flux[mask]
    mask2  = np.where((m_wave>xrange[0]) & (m_wave<xrange[1]))
    m_wave = m_wave[mask2]
    m_flux = m_flux[mask2]/max(m_flux[mask2])

    fig = plt.figure(figsize=(24,8))                                                               
    ax  = fig.add_subplot(1,1,1) 

    plt.plot(wave, flux, label=spec.name, color='k', alpha=.9)
    plt.plot(m_wave, m_flux, label=grid+': '+str(params), color='r', alpha=.8)
    plt.xlim(xrange)
    plt.tight_layout()
    plt.legend(loc='lower left')
    plt.xlabel(r"$\lambda [AA]$")
    plt.ylabel(r"$f_\lambda$")
    if save == True:
            plt.savefig(output)
    plt.show()
    plt.close()


def optimizeChi(spec, **kwargs):

    """
    Find the best chi-square fitting model, for a particular grid type and between particular parameter range.
    @Jessica Birky 

    Input:  'grid'     : type of model; right now BTSETTLb or PHOENIX
            'parrange' : parameter range to search between. Possible ranges specified in documentation.
            'gridPath' : default: apogee_tools libraries folder
            ...more: see kwargs list

    Output: Spectrum objects with best chi-squareds
            'best_fits', 'best_chis', 'best_pars'
    """

    #Retrieve grid library from the 'libraries' folder inside apogee_tools
    grid     = kwargs.get('grid', 'BTSETTLb')
    grid_lib = BASE + '/libraries/' + grid + '_APOGEE.hdf5'

    #default parameter ranges depending on grid type
    if grid == 'PHOENIX':
        prange = [[2300, 3500], [2.5, 5.5], [-0.5, 0.5]]
    elif grid == 'BTSETTLb':
        prange = [[2200, 3200], [2.5, 5.5], [-0.5, 0.0]]
    else: 
        prange = [[2200, 3200], [2.5, 5.5], [-0.5, 0.0]]

    parrange = kwargs.get('parrange', prange)
    gridPath = kwargs.get('gridPath', grid_lib)
    xrange   = kwargs.get('xrange', [spec.wave[0],spec.wave[-1]])
    numBest  = kwargs.get('numBest', 5)
    plot     = kwargs.get('plot', True)
    save     = kwargs.get('save', False)
    output   = kwargs.get('output', 'best_fits.pdf')
    subCont  = kwargs.get('subCont', True) #Continuum will be subtracted from the model grids

    grid_specs = readModels(parrange=parrange, gridPath=gridPath, subCont=subCont, xrange=xrange)

    #Cut input spectrum to the specified xrange
    mask = np.where((spec.wave > xrange[0]) & (spec.wave < xrange[1]))
    new_wave = spec.wave[mask]
    new_flux = spec.flux[mask]
    new_sigm = spec.sigmas[mask]
    spec = ap.Spectrum(wave=new_wave, flux=new_flux, sigmas=new_sigm, name=spec.name)

    chi_vals = {}

    #Compute chi values for each grid spectrum; append together to dictionary
    for sp in grid_specs:
        mspec = ap.Spectrum(wave=sp.wave, flux=sp.flux, params=sp.params)

        chi, sp1, sp2 = compareSpectra(spec, mspec)
        chi_vals[sp2] = chi

    #Sort grids in dictionary by lowest chi-squared value
    ranked = []
    values = []
    for w in sorted(chi_vals, key=chi_vals.get):
        ranked.append(w)
        values.append(chi_vals[w])

    best_fits = ranked[0:numBest]
    best_chis = values[0:numBest]
    best_pars = [fit.params for fit in best_fits]

    #Plot and save best fits
    if plot == True:
        nplots = numBest
        fig_dim = [12, 3*nplots]

        fig, axs = plt.subplots(nplots, 1, figsize=fig_dim)
        fig.subplots_adjust(hspace=.2)
        axs = axs.ravel()

        for n in range(nplots):
            axs[n].plot(spec.wave, spec.flux, color='k', alpha=.9, label='chi-squared = '+str(best_chis[n]))
            axs[n].plot(best_fits[n].wave, best_fits[n].flux, color='r', alpha=.7, label=best_fits[n].params)
            
            axs[n].set_xlim(xrange)
            axs[n].set_ylim([.5, 1.4])

            axs[n].legend()
            axs[n].set_ylabel(r"$f_\lambda$")

        # plt.title(spec.name, fontsize=12)
        plt.xlabel(r"$\lambda [AA]$")

        if save == True:
            plt.savefig(output)
        plt.show()
        plt.close()

    return best_fits, best_chis, best_pars


def optimizeRV(sp1, sp2, **kwargs):

    """
    Radial velocity computed by cross correlation with a template spectrum
    Radial velocity of sp1 will be shifted to align with sp2
    @Jessica Birky

    Input:  'sp1': data spectrum object
            'sp2': model/template spectrum object
    
    Output: 'rv'     : cross-correlated radial velocity  
            'cc_sp1' : rv shifted spectrum 1
            'sp2'    : same as input
    """

    from PyAstronomy import pyasl

    rv_lim = kwargs.get('rv_lim', [-150, 150])
    xrange = kwargs.get('xrange', [15200,15700])

    #Scale template sp2 to the data sp1
    chi, sp1, sp2 = compareSpectra(sp1, sp2)

    #Cut spectrum to calculate rv from specific region
    mask = np.where((sp1.wave > xrange[0]) & (sp1.wave < xrange[-1]))

    #Read in wavelength and flux from data and model spectrum objects
    wave1, flux1 = sp1.wave[mask], sp1.flux[mask]
    wave2, flux2 = sp2.wave[mask], sp2.flux[mask]

    #Cross correlate sp1 (data) to sp2 (model/template); 
    #Return radial velocity and cross-correlated spectrum
    rvs, ccs = pyasl.crosscorrRV(wave1, flux1, wave2, flux2, rv_lim[0], rv_lim[1], 30./50., skipedge=300)

    maxind = np.argmax(ccs)
    rv = rvs[maxind]

    #Shift wavelength of first spectrum by calculated radial velocity
    shift_wave = _rvShift(wave1, rv=-rv)

    # fig = plt.figure(figsize=(24,8))
    plt.subplot(2,1,1)
    plt.plot(wave2, flux2, color='k', label='template')
    if rv >= 0:
        plt.plot(wave1, flux1-.3, color='r', label='data rv = '+str(rv))
    else: 
        plt.plot(wave1, flux1-.3, color='b', label='data rv = '+str(rv))
    plt.plot(shift_wave, flux1, color='m', label='data rest frame')

    plt.xlim(xrange)
    plt.ylim([.2, 1.8])
    plt.legend(loc='upper left')
    plt.tight_layout()

    plt.subplot(2,1,2)
    plt.plot(rvs, ccs, label='Cross-correlation function')

    plt.legend(loc='upper left')
    plt.show()
    plt.close()

    #Turn cross-correlated spectrum into spectrum object
    cc_sp1 = Spectrum(wave=shift_wave, flux=flux1, sigmas=sp1.sigmas, name=sp1.name)

    return rv, cc_sp1, sp2


def smoothVSINI(mspec, **kwargs):

    """
    Add vsini to spectrum using PyAstronomy.
    @Jessica Birky

    Input:  'limb'  : limb darkening coefficient
            'vsini' : rotational velocity (km/s)
            'xlim'  : specify a wavelength region to perform over
            'plot'  : plot broadened spectrum

    Output: 'rot_spec' : broadened spectrum object
    """

    from PyAstronomy import pyasl
    from scipy.interpolate import interp1d

    #Input limb-darkening coefficient, vsini, plotting x range
    limb  = kwargs.get('limb', 0.6)
    vsini = kwargs.get('vsini', 1)
    xlim  = kwargs.get('xlim')

    #Read in wavelength and flux from model spectrum object
    m_wave, m_flux, params = mspec.wave, mspec.flux, mspec.params

    #Create evenly spaced wavelength array
    bounds = [m_wave[0], m_wave[-1]]
    npts   = len(m_wave)
    n_wave = np.linspace(bounds[0], bounds[1], npts)
    f = interp1d(m_wave, m_flux, bounds_error=False, fill_value=0.)
    n_flux = f(n_wave)

    #Perform rotational broadening using PyAstronomy; return flux
    rflux = pyasl.rotBroad(n_wave, n_flux, limb, vsini)

    #Save broadened spectrum to new spectrum object
    rot_spec = Spectrum(wave=n_wave, flux=rflux, params=params)

    if kwargs.get('plot', False) == True:
        plt.plot(m_wave, m_flux, color='k', label=params, alpha=.9)
        plt.plot(m_wave, rflux-.4, color='g', label='vsini = '+str(vsini), alpha=.8)
        plt.legend(loc='lower left')
        plt.xlabel(r"$\lambda [AA]$")
        plt.ylabel(r"$f_\lambda$")
        if 'xlim' in kwargs:
            plt.xlim(xlim)
        plt.show()
        plt.close()

    return rot_spec


def optimizeVSINI(spec, mspec, **kwargs):

    #Find best visini fit
    #I don't remember if this actually works

    step   = kwargs.get('step', 1)
    limits = kwargs.get('limits', [1, 50])
    iter_vals = np.arange(limits[0], limits[1], step)
    xrange = kwargs.get('xrange', [15200,16940])

    #Read in wavelength and flux from model spectrum object
    m_wave, m_flux, params = mspec.wave, mspec.flux, mspec.params

    #Cut model spectrum to range of data spectrum
    mask = np.where((m_wave > xrange[0]) & (m_wave < xrange[1]))
    m_wave, m_flux = m_wave[mask], m_flux[mask]
    mspec = Spectrum(wave=m_wave, flux=m_flux, params=params)

    chi_vals = {}

    #Test each vsini value, evaluate by chi-square
    for vsini in iter_vals:
        rot_spec = smoothVSINI(mspec, vsini=vsini)
        rot_spec = Spectrum(wave=rot_spec.wave, flux=rot_spec.flux, params=rot_spec.params, vsini=vsini)
        chi, spec, mdl = compareSpectra(spec, rot_spec)
        chi_vals[mdl] = chi

    #Sort grids in dictionary by lowest chi-squared value
    ranked, values = [], []
    for w in sorted(chi_vals, key=chi_vals.get):
        ranked.append(w)
        values.append(chi_vals[w])

    vsinis = [r.vsini for r in ranked]
    for i in range(len(ranked)):
        print(vsinis[i], values[i])

    #Find best fitting vsini spectrum, lowest chi-squared
    lowest = np.where(values == min(values))[0][0]
    best_chisq = values[lowest]
    best_vsini = vsinis[lowest]
    best_mspec = ranked[lowest]

    plt.subplot(2,1,1)
    plt.plot(spec.wave, spec.flux, color='k', alpha=.8, label=str(spec.name))
    plt.plot(best_mspec.wave, best_mspec.flux, color='g', label='vsini = '+str(best_mspec.vsini)+' km/s')
   
    plt.xlim(xrange)
    plt.ylim([.6, 1.4])
    plt.legend(loc='upper left')
    plt.tight_layout()

    plt.subplot(2,1,2)
    # l1, l2, [list(x) for x in zip(*sorted(zip(values, vsinis), key=lambda pair: pair[0]))]
    plt.plot(vsinis, values, label='vsini vs. chi')

    plt.legend(loc='upper left')
    plt.show()
    plt.close()

    return best_chisq, best_vsini, best_mspec

