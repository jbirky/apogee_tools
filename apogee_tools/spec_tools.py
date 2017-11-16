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


def calcScale(sp1, sp2):

    """
    Find the scale factor for the second spectrum to minimize chi-squared value
    @Jessica Birky

    Input:  two spectrum objects

    Output: compute scale factor
    """

    wave1, flux1 = sp1.wave, sp1.flux
    # wave2, flux2 = sp2.wave, sp2.flux
    unc = sp1.sigmas

    # interpolate spectrum 2 to the array size of spectrum 1
    wave2 = wave1
    flux2 = np.interp(wave1, sp2.wave, sp2.flux)

    a = sum(flux2*flux1/unc**2)
    b = sum(flux2**2/unc**2)

    scale = a/b

    return scale


def compareSpectra(sp1, sp2, **kwargs):

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
        
    flux1_vals = [x for x in flux1 if str(x) != 'nan']
    flux1 = flux1/max(flux1_vals)
    flux1 = np.array([0 if str(x) == 'nan' else x for x in flux1])

    flux2_vals = [x for x in flux2 if str(x) != 'nan']
    flux2 = flux2/max(flux2_vals)
    flux2 = np.array([0 if str(x) == 'nan' else x for x in flux2])

    #Create a new spectrum object for sp2
    sp2 = ap.Spectrum(wave=sp2.wave, flux=flux2, params=sp2.params)

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
