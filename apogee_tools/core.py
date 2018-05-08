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

#Get the path of apogee_tools file
FULL_PATH  = os.path.realpath(__file__)
BASE, NAME = os.path.split(FULL_PATH)

AP_PATH = os.environ['APOGEE_DATA']


class Spectrum:

    """
    Master Spectrum class.
    """

    def __init__(self, **kwargs):

            self.wave  = kwargs.get('wave', [])
            self.flux  = kwargs.get('flux', [])
            self.error = kwargs.get('error', [0 for i in range(len(self.wave))])
            self.mask  = kwargs.get('mask', [])
            self.model = kwargs.get('model', [])
            self.sky   = kwargs.get('sky', [])
            self.name  = kwargs.get('name', 'Spectrum')  
            self.param = kwargs.get('param', [])   

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


    def plot(self, **kwargs):

        xrange = kwargs.get('xrange', [self.wave[0], self.wave[-1]])
        yrange = kwargs.get('yrange', [min(self.flux)-.2, max(self.flux)+.2])
        rv     = kwargs.get('rv', 0)
        items  = kwargs.get('items', ['spec'])
        title  = kwargs.get('title')
        save   = kwargs.get('save', False)
        output = kwargs.get('output', str(self.name) + '.pdf')

        highlight = kwargs.get('highlight')
        hcolor    = kwargs.get('hcolor', 'b')
        
        rv_wave = ap.rvShift(self.wave, rv=rv)
        
        fig = plt.figure(figsize=(16,4))                                                               
        ax  = fig.add_subplot(1,1,1) 

        # Plot masked spectrum
        if ('spectrum' in items) or ('spec' in items):
            plt.plot(rv_wave, self.flux, color='k', alpha=.8, linewidth=1, label=self.name)

        # Plot spectrum error
        if 'error' in items:
            plt.plot(self.wave, self.error, color='c', linewidth=1, alpha=.6)
        
        # Plot read in model
        if 'model' in items:
            plt.plot(self.wave, self.model, color='r', alpha=.8, linewidth=1, label='Model')

        # Highlight specified bands
        if 'highlight' in kwargs:
            for h in highlight:
               plt.axvspan(h[0], h[1], color=hcolor, alpha=0.1)

        # Plot and label atomic lines
        if ('lines' in items) or ('line_list' in kwargs):
            line_list = kwargs.get('line_list', ap.lines)
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

        # Plot multiple dictionaries of lines in different colors
        if 'line_lists' in kwargs:
            line_lists = kwargs.get('line_lists', [ap.lines])
            list_labels = kwargs.get('list_labels', ['list'+str(i) for i in range(len(line_lists))])

            colors = ['r', 'b', 'g', 'p', 'c', 'y']
            cindex = 0

            for line_list in line_lists:
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

                            plt.axvline(x=feature, ymin=plot_ypos_min, ymax=plot_ypos_max, linewidth=1, color=colors[cindex])
                            plt.text(feature, ypos-.2, lines, rotation=90, ha='center', color='k', fontsize=8)

                cindex += 1

        
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


class Apogee(Spectrum):

    """ 
    APOGEE Spectrum class for reading in apogee fits files; includes features for aspcap, apStar, and apVisit.
    Read in a file, or write in your own parameters.
    @Jessica Birky, Christian Aganze

    Input: 'id'   : 2MASS name. If downloaded the spectrum will be loaded from 'APOGEE_DATA' environmental variable
           'type' : aspcap, apStar, or apVisit
    """

    def __init__(self, **kwargs):  

        self.d_type = kwargs.get('type')

        # Get designation of 2MASS object
        input_id  = kwargs.get('id')            
        name = ap.formatDesignation(input_id)

        if self.d_type == 'aspcap':

            file_dr14 = '%s/aspcap_data/aspcapStar-r8-l31c.2-%s.fits' %(AP_PATH, name)

            if (os.path.exists(file_dr14) == False):
                ap.download(name, type='aspcap')

            openFile = fits.open(file_dr14)

            crval = openFile[1].header['CRVAL1']
            cdelt = openFile[1].header['CDELT1']

            wave  = np.array(pow(10, crval+cdelt*np.arange(openFile[1].header['NAXIS1'])))    
            flux  = np.array(openFile[1].data)
            error = np.array(openFile[2].data)
            model = np.array(openFile[3].data)
            param = openFile[4].data['PARAM'] 

            super().__init__(wave=wave, flux=flux, error=error, model=model, param=param, name=name)

        elif self.d_type == 'apstar':

            file_dr14 = '%s/apstar_data/apStar-r8-%s.fits' %(AP_PATH, name)

            if (os.path.exists(file_dr14) == False):
                ap.download(name, type='apstar')

            openFile = fits.open(file_dr14)
        
            #conversion from pixel to wavelength, info available in the hdu header
            crval = openFile[1].header['CRVAL1']
            cdelt = openFile[1].header['CDELT1']

            wave  = np.array(pow(10, crval+cdelt*np.arange(openFile[1].header['NAXIS1'])))           
            flux  = np.array(openFile[1].data)[0]
            error = np.array(openFile[2].data)

            super().__init__(wave=wave, flux=flux, error=error, name=name)

        elif self.d_type == 'apvisit':

            visit = kwargs.get('visit', 1)
            file  = '%s/apvisit_data/apVisit-%s-%s.fits' %(AP_PATH, name, visit)

            if os.path.exists(file) == False:
                ap.download(name, type='apvisit')

            openFile = fits.open(file)

            flux  = openFile[1].data
            error = openFile[2].data
            wave  = openFile[4].data

            # Combine the data from the three chips into one list
            wave = np.array(list(wave[0]) + list(wave[1]) + list(wave[2]))
            flux = np.array(list(flux[0]) + list(flux[1]) + list(flux[2]))
            error = np.array(error.data)

            super().__init__(wave=wave, flux=flux, error=error, name=name)
        
        elif self.d_type == 'ap1d':

            visit = kwargs.get('visit', 1)
    
            ap1d_dir = AP_PATH + '/{}_data/' .format(self.d_type)
            
            # Look up fiber number for given 2MASS name and its visit number
            ap_id, plates, mjds, fibers = ap.searchVisits(id_name=name)
            fiber = fibers[visit-1] #-1 for python zero indexing
            
            #get names of files for each chip
            ChipA = ap1d_dir + 'ap1d-{}-{}-a.fits'.format(name, visit)
            ChipB = ap1d_dir + 'ap1d-{}-{}-b.fits'.format(name, visit)
            ChipC = ap1d_dir + 'ap1d-{}-{}-c.fits'.format(name, visit)
            
            if (os.path.exists(ChipA) == False) or (os.path.exists(ChipB) == False) or (os.path.exists(ChipC) == False):
                ap.download(name, type='ap1d', visit=visit, frame=1)

            # Get wavelength, flux, and error arrays of spectrum
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
            
            # Concatenate arrays from 3 chips and reverse lists so that lowest wavelength is first
            wave = np.concatenate([wave1, wave2, wave3])[::-1]
            flux = np.concatenate([flux1, flux2, flux3])[::-1]
            error = np.concatenate([err1, err2, err3])[::-1]

            super().__init__(wave=wave, flux=flux, error=error, name=name)


class Nirspec(Spectrum):

    def __init__(self, **kwargs): 

        name = kwargs.get('id')
        path = kwargs.get('path', AP_PATH)
        order = kwargs.get('order')

        file = path + '/' + name + '_' + str(order) + '_all.fits'

        try:
            openFile = fits.open(file, ignore_missing_end=True)
        except:
            print(file, 'not found.')

        wave  = openFile[0].data
        flux  = openFile[1].data
        error = openFile[2].data

        super().__init__(wave=wave, flux=flux, error=error, name=name)

        
class ModelGrid(Spectrum):

    def __init__(self, **kwargs): 

        self.type = kwargs.get('type', 'BTSETTL')

        # read in model
        # super().__init__(wave=wave, flux=flux, error=error, name=name)