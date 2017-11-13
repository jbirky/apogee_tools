import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import h5py
from scipy.interpolate import interp1d
from scipy.integrate import trapz
import apogee_tools as ap

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
        # print(parameters)
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


def integralResample(xh, yh, xl, nsamp=100):

    '''
    Function from SPLAT. See: https://github.com/aburgasser/splat

    :Purpose: A 1D integral smoothing and resampling function that attempts to preserve total flux. Usese
    scipy.interpolate.interp1d and scipy.integrate.trapz to perform piece-wise integration

    Required Inputs:

    :param xh: x-axis values for "high resolution" data
    :param yh: y-axis values for "high resolution" data
    :param xl: x-axis values for resulting "low resolution" data, must be contained within high resolution and have fewer values

    Optional Inputs:

    :param nsamp: Number of samples for stepwise integration

    Output:

    y-axis values for resulting "low resolution" data

    :Example:
    >>> # a coarse way of downsampling spectrum
    >>> import splat, numpy
    >>> sp = splat.Spectrum(file='high_resolution_spectrum.fits')
    >>> w_low = numpy.linspace(numpy.min(sp.wave.value),numpy.max(sp.wave.value),len(sp.wave.value)/10.)
    >>> f_low = splat.integralResample(sp.wave.value,sp.flux.value,w_low)
    >>> n_low = splat.integralResample(sp.wave.value,sp.noise.value,w_low)
    >>> sp.wave = w_low*sp.wave.unit
    >>> sp.flux = f_low*sp.flux.unit
    >>> sp.noise = n_low*sp.noise.unit
    '''

    # check inputs
    if xl[0] < xh[0] or xl[-1] > xh[-1]: raise ValueError('\nLow resolution x range {} to {} must be within high resolution x range {} to {}'.format(xl[0],xl[-1],xh[0],xh[-1]))
    if len(xl) > len(xh): raise ValueError('\nTarget x-axis must be lower resolution than original x-axis')

    # set up samples
    xs = [np.max([xh[0],xl[0]-0.5*(xl[1]-xl[0])])]
    for i in range(len(xl)-1): xs.append(xl[i]+0.5*(xl[i+1]-xl[i]))
    xs.append(np.min([xl[-1]+0.5*(xl[-1]-xl[-2]),xh[-1]]))

    f = interp1d(xh,yh)

    # integral loop
    ys = []
    for i in range(len(xl)):
        dx = np.linspace(xs[i],xs[i+1],nsamp)
        ys.append(trapz(f(dx),x=dx)/trapz(np.ones(nsamp),x=dx))

    return ys

