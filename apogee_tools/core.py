from __future__ import print_function, division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
rc('font', family='serif')
from operator import itemgetter
import itertools
from itertools import groupby, chain
import os
import h5py
import json

from astropy.table import Table
from astropy.io import fits, ascii
from astropy import units as u

import apogee_tools as ap
from apogee_tools.spec_tools import _rvShift
#from libraries import features

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

            file_dr13 = '%s/aspcap_data/aspcapStar-r6-l30e.2-%s.fits' %(AP_PATH, spec_id)
            file_dr14 = '%s/aspcap_data/aspcapStar-r8-l31c.2-%s.fits' %(AP_PATH, spec_id)

            """ ASPCAP file info:
            HDU0: The Primary Header
            HDU1: Spectrum array
            HDU2: Error array
            HDU3: Best fit spectrum
            HDU4: ASPCAP data table
            Information found at https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/LOCATION_ID/aspcapStar.html 
            """

            if (os.path.exists(file_dr13) == False) & (os.path.exists(file_dr14) == False):
                ap.download(self.name, type='aspcap')

            try:
                openFile = fits.open(file_dr14)
            except:
                openFile = fits.open(file_dr13)

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

            file_dr13 = '%s/apstar_data/apStar-r6-%s.fits' %(AP_PATH, spec_id)
            file_dr14 = '%s/apstar_data/apStar-r8-%s.fits' %(AP_PATH, spec_id)

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

            if (os.path.exists(file_dr13) == False) & (os.path.exists(file_dr14) == False):
                ap.download(self.name, type='apstar')

            try:
                openFile = fits.open(file_dr14)
            except:
                openFile = fits.open(file_dr13)

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

            visit = kwargs.get('visit', 1)
            self.file  = '%s/apvisit_data/apVisit-%s-%s.fits' %(AP_PATH, spec_id, visit)

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

            if os.path.exists(self.file) == False:
                ap.download(self.name, type='apvisit')

            hdu   = fits.open(self.file)
            flux  = hdu[1].data
            error = hdu[2].data
            wave  = hdu[4].data

            #Combine the data from the three chips into one list
            self.wave = np.array(list(wave[0]) + list(wave[1]) + list(wave[2]))
            self.flux = np.array(list(flux[0]) + list(flux[1]) + list(flux[2]))

            self.sigmas = np.array(error.data)
            self.noise  = self.sigmas
            # self.sky    = np.array(sky.data)
        
        elif self.d_type == 'ap1d':

            """ AP1D file info:
            HDU0: master header
            HDU1: image (ADU) [FLOAT]
            HDU2: error (ADU) [FLOAT]
            HDU3: flag mask [INT*2]
            HDU4: wavelength array
            HDU5: wavelength coefficients array
            Information found at: https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/red/MJD5/ap1D.html
            """

            visit = kwargs.get('visit', 1)
    
            ap1d_dir = AP_PATH + '/{}_data/' .format(self.d_type)
            
            #look up fiber number for given 2MASS name and its visit number
            ap_id, plates, mjds, fibers = ap.searchVisits(id_name=self.name)
            fiber = fibers[visit-1] #-1 for python zero indexing
            
            #get names of files for each chip
            ChipA = ap1d_dir + 'ap1d-{}-{}-a.fits'.format(self.name, visit)
            ChipB = ap1d_dir + 'ap1d-{}-{}-b.fits'.format(self.name, visit)
            ChipC = ap1d_dir + 'ap1d-{}-{}-c.fits'.format(self.name, visit)
            
            if (os.path.exists(ChipA) == False) or (os.path.exists(ChipB) == False) or (os.path.exists(ChipC) == False) :
                ap.download(self.name, type='ap1d', visit=visit, frame=1)

            #get wavelength, flux, and error arrays of spectrum
            hdu_a = fits.open(ChipA)
            t1    = hdu_a[1].data[fiber] # assuming the first extension is a table
            wave1 = hdu_a[4].data[fiber]
            flux1 = t1*1e-17
            err1  = hdu_a[2].data[fiber]
            
            hdu_b = fits.open(ChipB)
            t2    = hdu_b[1].data[fiber] # assuming the first extension is a table
            wave2 = hdu_b[4].data[fiber]
            flux2 = t2*1e-17
            err2  = hdu_b[2].data[fiber]
            
            hdu_c = fits.open(ChipC)
            t3    = hdu_c[1].data[fiber] # assuming the first extension is a table
            wave3 = hdu_c[4].data[fiber]
            flux3 = t3*1e-17
            err3  = hdu_c[2].data[fiber]
            
            #concatenate arrays from 3 chips and reverse lists so that lowest wavelength is first
            self.wave = np.concatenate([wave1, wave2, wave3])[::-1]
            self.flux = np.concatenate([flux1, flux2, flux3])[::-1]
            self.sigmas = np.concatenate([err1, err2, err3])[::-1]
            self.noise  = self.sigmas

        else:
            #Save spectrum values into the spectrum object class
            self.wave   = kwargs.get('wave', [])
            self.flux   = kwargs.get('flux', [])
            self.sky    = kwargs.get('sky', [])
            self.sigmas = kwargs.get('sigmas', [0 for i in range(len(self.wave))])
            self.noise  = self.sigmas
            self.model  = kwargs.get('model', [])
            self.name   = kwargs.get('name', 'spectrum')  
            self.params = kwargs.get('params', [])   
            self.vsini  = kwargs.get('vsini', []) 
            self.type   = 'input'   
        
        
    def plot(self, **kwargs):

        xrange = kwargs.get('xrange', [self.wave[0], self.wave[-1]])
        yrange = kwargs.get('yrange', [min(self.flux)-.2, max(self.flux)+.2])
        rv     = kwargs.get('rv', 0)
        items  = kwargs.get('items', ['spec'])
        title  = kwargs.get('title')
        save   = kwargs.get('save', False)
        output = kwargs.get('output', str(self.name) + '.pdf')
        
        rv_wave = _rvShift(self.wave, rv=rv)
        
        fig = plt.figure(figsize=(16,4))                                                               
        ax  = fig.add_subplot(1,1,1) 

        #Plot masked spectrum
        if ('spectrum' in items) or ('spec' in items):
            plt.plot(rv_wave, self.flux, color='k', alpha=.8, linewidth=1, label=self.name)

        #Plot spectrum noise
        if 'noise' in items:
            plt.plot(self.wave, self.sigmas, color='c', linewidth=1, alpha=.6)

        #Plot continuum
        if ('cont' in items) or ('continuum' in items):
            plt.plot(rv_wave, self.cont, color='m', linewidth=1, alpha=.8)

        #Plot aspcap model 
        if 'apModel' in items:
            plt.plot(self.wave, self.apModel, color='r', alpha=.8, linewidth=1, label='ASPCAP Model')
        
        #Plot read in model
        if 'model' in items:
            plt.plot(self.wave, self.model, color='r', alpha=.8, linewidth=1, label='Model')

        #Plot and label atomic lines
        if 'lines' in items:
            line_list = ap.lines
            line_names = line_list.keys()

            for lines in line_names:
                for feature in line_list[lines]:

                    if (feature <= xrange[1]) & (feature >= xrange[0]):
                        # determine position of the line and label based on pixel of the spectrum
                        xpos = min(self.wave, key=lambda x:abs(x - feature))
                        index = list(self.wave).index(xpos)
                        ypos = self.flux[index]
                        plot_ypos_min = (ypos - yrange[0] -.15)/(yrange[1] - yrange[0])
                        plot_ypos_max = (ypos - yrange[0] -.1)/(yrange[1] - yrange[0])

                        plt.axvline(x=feature, ymin=plot_ypos_min, ymax=plot_ypos_max, linewidth=1, color = 'g')
                        plt.text(feature, ypos-.2, lines, rotation=90, ha='center', color='b', fontsize=8)

        #Plot and highlight molecular bands
        # if 'bands' in items:

        #Plot piece-wise model segments in different colors
        if 'model' in kwargs:
            colors = ['m', 'b', 'g', 'y', '#ffa500', 'r']
            model  = kwargs.get('model')
            labels = kwargs.get('labels')

            for i in range(len(model)):
                plt.plot(model[i]['wl'], model[i]['model'], label=labels[i],color=colors[i], alpha=.8)
        
        plt.legend(loc='upper right', fontsize=12)
        
        
        plt.xlim(xrange)
        plt.ylim(yrange)    
    
        minor_locator = AutoMinorLocator(5)
        ax.xaxis.set_minor_locator(minor_locator)
        # plt.grid(which='minor') 
    
        plt.xlabel(r'$\lambda$ [$\mathring{A}$]', fontsize=18)
        plt.ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]', fontsize=18)
        if title != None:
            plt.title(title, fontsize=20)
        plt.tight_layout()

        if save == True:
            plt.savefig(output)

        plt.show()
        plt.close()

    # Add more functions here to easily manipulate spectrum object (like splat):

    def mask(self, **kwargs):

        """
        Mask all pixels that are out of the specified sigma cutoff range specified.

        Input: 'sigma'        : [lower cuttoff, upper cutoff]
               'pixel_buffer' : [lower mask pixel buffer, upper mask pixel buffer]
        """

        sigma = kwargs.get('sigma', [2,1])
        pixel_buffer = kwargs.get('pixel_buffer', [0,2])

        fmean = np.mean(self.flux)
        fstd  = np.std(self.flux)

        #Find outlier flux and bad pixels 
        cut_low  = np.where(self.flux <= fmean - sigma[0]*fstd)[0]
        cut_high = np.where(self.flux >= fmean + sigma[1]*fstd)[0]

        group_low_cut = []
        group_high_cut = []

        for k, g in itertools.groupby(enumerate(cut_low), lambda x:x[0]-x[1]):
            group_low_cut.append(list(map(itemgetter(1), g)))

        for k, g in itertools.groupby(enumerate(cut_high), lambda x:x[0]-x[1]):
            group_high_cut.append(list(map(itemgetter(1), g)))

        cuts = []
        for pixels in group_low_cut:
            pix1, pixn = pixels[0], pixels[-1]
            for b in range(pixel_buffer[0]):
                pixels.append(pix1 - (b+1))
                pixels.append(pixn + (b+1))
            cuts.append(pixels)

        for pixels in group_high_cut:
            pix1, pixn = pixels[0], pixels[-1]
            for b in range(pixel_buffer[1]):
                pixels.append(pix1 - (b+1))
                pixels.append(pixn + (b+1))
            cuts.append(pixels)

        cuts = list(itertools.chain(*cuts))
        self.flux[cuts] = np.nan


    def shift_rv(self, **kwargs):

        """
        Shift wavelenth of spectrum object by specific rv.
        @Jessica Birky

        Input: 'rv' : radial velocity (km/s)
        """

        rv = kwargs.get('rv', -80) 
        self.wave = _rvShift(self.wave, rv=rv)


    def broaden(self, **kwargs):

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
