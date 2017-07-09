from __future__ import print_function, division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from operator import itemgetter
import itertools
from itertools import groupby, chain
import os
import h5py
import json

from astropy.table import Table
from astropy.io import fits, ascii
from astropy import units as u

from apogee_tools.spec_tools import _rvShift

#Get the path of apogee_tools file
FULL_PATH  = os.path.realpath(__file__)
BASE, NAME = os.path.split(FULL_PATH)

AP_PATH = os.environ['APOGEE_DATA']

class Spectrum():

    """ 
    Spectrum class for reading in apogee fits files; includes features for aspcap, apStar, and apVisit.
    Read in a file, or write in your own parameters.
    @Jessica Birky, Christian Aganze

    Input: 'id'   : 2MASS name. If downloaded the spectrum will be loaded from 'APOGEE_DATA' environmental variable
           'type' : aspcap, apStar, or apVisit
    
    Spectrum object contains parameters: 
        wave, flux, sky, noise (or sigmas), name, continuum, params, apModel (aspcap model)
    """

    def __init__(self, **kwargs):  

        self.d_type = kwargs.get('type')

        self.name   = kwargs.get('id')
        spec_id = self.name

        if self.d_type == 'aspcap':

            file = '%s/aspcap_data/aspcapStar-r6-l30e.2-%s.fits' %(AP_PATH, spec_id)

            """ ASPCAP file info:
            HDU0: The Primary Header
            HDU1: Spectrum array
            HDU2: Error array
            HDU3: Best fit spectrum
            HDU4: ASPCAP data table
            Information found at https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/LOCATION_ID/aspcapStar.html 
            """

            openFile = fits.open(file)

            self.HDU0 = openFile[0]
            self.HDU1 = openFile[1]
            self.HDU2 = openFile[2]
            self.HDU3 = openFile[3]
            self.HDU4 = openFile[4]

            master    = self.HDU0
            spectra   = self.HDU1
            error     = self.HDU2
            best_fit  = self.HDU3
            aspcap_dt = self.HDU4
        
            #conversion from pixel to wavelength, info available in the hdu header
            crval = spectra.header['CRVAL1']
            cdelt = spectra.header['CDELT1']
            wave  = np.array(pow(10, crval+cdelt*np.arange(spectra.header['NAXIS1']))/10000)*u.micron #microns
            
            # convert fluxes from  (10^-17 erg/s/cm^2/Ang) to  ( erg/s/cm^2/Mircon)
            spectras = [1e-13*np.array(f)*u.erg/u.s/u.centimeter**2/u.micron for f in spectra.data]
                    
            self.wave    = 10000*np.array(wave.value)
            self.flux    = np.array(spectra.data)
            self.sigmas  = np.array(error.data)
            self.apModel = np.array(best_fit.data)

            self.avgFlux = np.mean(self.flux)
            self.stdFlux = np.std(self.flux)
            
            #Find outlier flux and bad pixels
            self.smooth_flux = self.flux       
            self.smooth_flux[self.smooth_flux <= 1-.6*self.stdFlux] = 0

            #Mask outliers
            mask = np.where(self.smooth_flux == 0)
            
            self.wave    = np.delete(self.wave, list(mask))
            self.flux    = np.delete(self.flux,list(mask))
            self.sigmas  = np.delete(self.sigmas, list(mask))
            self.noise   = self.sigmas
            self.apModel = np.delete(self.apModel, list(mask))

            #Obtain aspcap parameters
            self.params = aspcap_dt.data['PARAM'] 

        elif self.d_type == 'apstar':

            file = '%s/apstar_data/apStar-r6-%s.fits' %(AP_PATH, spec_id)

            """ APSTAR file info:
            HDU0: master header with target information
            HDU1: spectra: combined and individual
            HDU2: error spectra
            HDU3: mask spectra
            HDU4: sky spectra
            HDU5: sky error spectra
            HDU6: telluric spectra
            HDU7: telluric error spectra
            HDU8: table with LSF coefficients
            HDU9: table with RV/binary information
            Information found at https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/TELESCOPE/LOCATION_ID/apStar.html 
            """

            openFile = fits.open(file)

            self.HDU0 = openFile[0]
            self.HDU1 = openFile[1]
            self.HDU2 = openFile[2]
            self.HDU3 = openFile[3]
            self.HDU4 = openFile[4]
            self.HDU5 = openFile[5]
            self.HDU6 = openFile[6]
            self.HDU7 = openFile[7]
            self.HDU8 = openFile[8]
            self.HDU9 = openFile[9]

            master   = self.HDU0
            spectra  = self.HDU1
            error    = self.HDU2
            mask     = self.HDU3
            sky      = self.HDU4
            sky_err  = self.HDU5
            telluric = self.HDU6
            tell_err = self.HDU7
            lsf_coef = self.HDU8
            rv_info  = self.HDU9
        
            #conversion from pixel to wavelength, info available in the hdu header
            crval = spectra.header['CRVAL1']
            cdelt = spectra.header['CDELT1']
            wave  = np.array(pow(10, crval+cdelt*np.arange(spectra.header['NAXIS1']))/10000)*u.micron #microns
            
            # convert fluxes from  (10^-17 erg/s/cm^2/Ang) to  ( erg/s/cm^2/Mircon)
            spectras = [1e-13*np.array(f)*u.erg/u.s/u.centimeter**2/u.micron for f in spectra.data]
                    
            self.wave    = 10000*np.array(wave.value)
            self.flux    = np.array(spectra.data)
            self.sigmas  = np.array(error.data)
            self.sky     = np.array(sky.data)
 
            self.avgFlux = np.mean(self.flux)
            self.stdFlux = np.std(self.flux)

            mask = mask.data
            
            self.wave    = np.delete(self.wave, list(mask))
            self.flux    = np.delete(self.flux,list(mask))
            self.sigmas  = np.delete(self.sigmas, list(mask))
            self.noise   = self.sigmas

        elif self.d_type == 'apvisit':

            self.file  = '%s/apstar_data/apVisit-r6-%s.fits' %(AP_PATH, spec_id)

            """ APVISIT file info:
            HDU0: master header with target information
            HDU1: spectra: combined and individual
            HDU2: error
            HDU3: mask
            HDU4: wavelength array
            HDU5: sky
            HDU6: sky error
            HDU7: telluric
            HDU8: telluric error
            HDU9: wavelength coeficients
            HDU10: LSF coeficients
            Information found at: https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/TELESCOPE/PLATE_ID/MJD5/apVisit.html
            """

            openFile = fits.open(self.file)
            
            self.HDU0  = openFile[0]
            self.HDU1  = openFile[1]
            self.HDU2  = openFile[2]
            self.HDU3  = openFile[3]
            self.HDU4  = openFile[4]
            self.HDU5  = openFile[5]
            self.HDU6  = openFile[6]
            self.HDU7  = openFile[7]
            self.HDU8  = openFile[8]
            self.HDU9  = openFile[9]
            self.HDU10 = openFile[10]

            hdu     = self.HDU0
            flux    = self.HDU1
            error   = self.HDU2
            wave    = self.HDU3
            sky     = self.HDU5
            sky_err = self.HDU6

            #Combine the data from the three chips into one list
            self.wave = np.array(list(wave[0]) + list(wave[1]) + list(wave[2]))
            self.flux = np.array(list(flux[0]) + list(flux[1]) + list(flux[2]))

            self.sigmas = np.array(error.data)
            self.noise  = self.sigmas
            self.sky    = np.array(sky.data)

        
        else:
            #Save spectrum values into the spectrum object class
            self.wave   = kwargs.get('wave', [])
            self.flux   = kwargs.get('flux', [])
            self.sky    = kwargs.get('sky', [])
            self.sigmas = kwargs.get('sigmas', [0 for i in range(len(self.wave))])
            self.noise  = self.sigmas
            self.model  = kwargs.get('model', [])
            self.name   = kwargs.get('name', Spectrum)  
            self.params = kwargs.get('params', [])   
            self.vsini  = kwargs.get('vsini', []) 
            self.type   = 'input'   
        
        
    def plot(self, **kwargs):

        xrange = kwargs.get('xrange', [self.wave[0],self.wave[-1]])
        rv     = kwargs.get('rv', 0)
        items  = kwargs.get('items', ['spec'])
        save   = kwargs.get('save', False)
        output = kwargs.get('output', self.name + '.pdf')
        
        rv_wave = _rvShift(self.wave, rv=rv)
        
        fig = plt.figure(figsize=(12,5))                                                               
        ax  = fig.add_subplot(1,1,1) 

        #Plot masked spectrum
        if ('spectrum' in items) or ('spec' in items):
            plt.plot(rv_wave, self.flux, color='k', alpha=.9, label=self.name)

        #Plot spectrum noise
        if 'noise' in items:
            plt.plot(wave, self.sigmas, color='c', alpha=.6)

        #Plot continuum
        if ('cont' in items) or ('continuum' in items):
            plt.plot(rv_wave, self.cont, color='m', alpha=.8)

        #Plot aspcap model 
        if 'apModel' in items:
            plt.plot(wave, self.apModel, color='r', alpha=.8, label='ASPCAP Model')
        
        #Plot read in model
        if 'model' in items:
            plt.plot(wave, self.model, color='r', alpha=.8, label='Model')

        #Plot piece-wise model segments in different colors
        if 'model' in kwargs:
            colors = ['m', 'b', 'g', 'y', '#ffa500', 'r']
            model  = kwargs.get('model')
            labels = kwargs.get('labels')

            for i in range(len(model)):
                plt.plot(model[i]['wl'], model[i]['model'], label=labels[i],color=colors[i], alpha=.8)
        
        plt.legend(loc='upper right')
        
        major_ticks = np.arange(15100, 17000, 200)
        minor_ticks = np.arange(15100, 17000, 50)
        ax.set_xticks(major_ticks)                                                       
        ax.set_xticks(minor_ticks, minor=True) 
        
        plt.xlim(xrange)     
    
        plt.xlabel(r'$\lambda$ [$\mathring{A}$]')
        plt.ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]')
        plt.tight_layout()

        if save == True:
            plt.savefig(output)

        plt.show()
        plt.close()

    # Add more functions here to easily manipulate spectrum object (like splat):

    def shift_rv(self):

        """
        Shift wavelenth of spectrum object by specific rv.
        @Jessica Birky

        Input: 'rv' : radial velocity (km/s)
        """

        rv = kwargs.get('rv', -80) 
        self.wave = _rvShift(wave, rv=rv)


    def broaden(self):

        """
        Add vsini to spectrum using PyAstronomy.
        @Jessica Birky

        Input:  'limb'  : limb darkening coefficient
                'vsini' : rotational velocity (km/s)
                'xlim'  : specify a wavelength region to perform over
                'plot'  : plot broadened spectrum
        """

        #Input limb-darkening coefficient, vsini, plotting x range
        limb  = kwargs.get('limb', 0)
        vsini = kwargs.get('vsini', 0)
        xlim  = kwargs.get('xlim')
        plot  = kwargs.get('plot', False)

        broad_spec = smoothVSINI(self, limb=limb, vsini=vsini, xlim=xlim, plot=plot)

        self.wave = broad_spec.wave
        self.flux = broad_spec.flux

