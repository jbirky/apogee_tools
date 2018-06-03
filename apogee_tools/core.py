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

# Unit standards
DEFAULT_WAVE_UNIT = u.micron
DEFAULT_FLUX_UNIT = u.erg/u.s/u.cm/u.cm/u.micron
DEFAULT_SED_UNIT  = u.erg/u.s/u.cm/u.cm

SPECTRAL_MODELS = {\
#    'gaia': {'folder': SPLAT_PATH+SPECTRAL_MODEL_FOLDER+'/gaia/', 'name': 'AMES GAIA', 'citation': 'Hauschildt et al. (1999)', 'bibcode': '1999ApJ...525..871H', 'altname': ['nextgen,hauschildt,hauschildt99,hauschildt1999'], 'rawfolder': HOME_FOLDER+'/models/phoenix/nextgen/fullres/', 'default': {'teff': 2000., 'logg': 5.0, 'z': 0.0}}, \
#    'btnextgen': {'instruments': {}, 'name': 'BT NextGen', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['nextgen-bt','btnextgen'], 'default': {'teff': 3000., 'logg': 5.0, 'z': 0.0, 'enrich': 0.}}, \
    'btsettl08': {'instruments': {'APOGEE-RAW'}, 'name': 'BTSettl08', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['allard','allard12','allard2012','btsettl','btsettled','btsettl08','btsettl2008','BTSettl2008'], 'default': {'teff': 2000., 'logg': 5.0, 'z': 0., 'enrich': 0.}}, \
#    'btsettl15': {'instruments': {}, 'name': 'BTSettl15', 'citation': 'Allard et al. (2015)', 'bibcode': '2015A&A...577A..42B', 'altname': ['allard15','allard2015','btsettl015','btsettl2015','BTSettl2015'],  'default': {'teff': 1500., 'logg': 5.0, 'z': 0.}}, \
#    'burrows06': {'instruments': {}, 'name': 'Burrows 2006', 'citation': 'Burrows et al. (2006)', 'bibcode': '2006ApJ...640.1063B', 'altname': ['burrows','burrows2006'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'cld': 'nc'}}, \
#    'cond01': {'instruments': {}, 'name': 'AMES Cond', 'citation': 'Allard et al. (2001)', 'bibcode': '2001ApJ...556..357A', 'altname': ['cond','cond-ames','amescond'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
#    'drift': {'instruments': {}, 'name': 'Drift', 'citation': 'Witte et al. (2011)', 'bibcode': '2011A&A...529A..44W', 'altname': ['witte','witte11','witte2011','helling'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.}}, \
#    'dusty01': {'instruments': {}, 'name': 'AMES Dusty', 'citation': 'Allard et al. (2001)', 'bibcode': '2001ApJ...556..357A', 'altname': ['dusty','dusty-ames','amesdusty'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.0}}, \
#    'madhusudhan11': {'instruments': {}, 'name': 'Madhusudhan 2011', 'citation': 'Madhusudhan et al. (2011)', 'bibcode': '2011ApJ...737...34M', 'altname': ['madhu','madhusudhan','madhu11','madhu2011','madhusudhan2011'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.,'cld': 'ae60', 'kzz': 'eq','fsed': 'eq'}}, \
#    'morley12': {'instruments': {}, 'name': 'Morley 2012', 'citation': 'Morley et al. (2012)', 'bibcode': '2012ApJ...756..172M', 'altname': ['morley','morley2012'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0., 'fsed': 'f5'}}, \
#    'morley14': {'instruments': {}, 'name': 'Morley 2014', 'citation': 'Morley et al. (2014)', 'bibcode': '2014ApJ...787...78M', 'altname': ['morley2014'], 'default': {'teff': 300., 'logg': 5.0, 'z': 0., 'fsed': 'f5', 'cld': 'h50'}}, \
#    'nextgen99': {'instruments': {}, 'name': 'Phoenix NextGen', 'citation': 'Hauschildt et al. (1999)', 'bibcode': '1999ApJ...525..871H', 'altname': ['nextgen,hauschildt,hauschildt99,hauschildt1999'], 'default': {'teff': 2000., 'logg': 5.0, 'z': 0.0}}, \
#    'saumon12': {'instruments': {}, 'name': 'Saumon 2012', 'citation': 'Saumon et al. (2012)', 'bibcode': '2012ApJ...750...74S', 'altname': ['saumon','saumon2012'], 'default': {'teff': 1000., 'logg': 5.0, 'z': 0.}}, \
#    'btcond': {'instruments': {}, 'name': 'BT Cond', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['dusty-cond','bt-cond'], 'rawfolder': '/Volumes/splat/models/btcond/ORIGINAL/', 'default': {'teff': 1500., 'logg': 5.0, 'z': 0.0, 'enrich': 0.0}}, \
#    'btdusty': {'instruments': {}, 'name': 'BT Dusty', 'citation': 'Allard et al. (2012)', 'bibcode': '2012RSPTA.370.2765A', 'altname': ['dusty-bt','bt-dusty'], 'rawfolder': '/Volumes/splat/models/btdusty/ORIGINAL/', 'default': {'teff': 1500., 'logg': 5.0, 'z': 0.0}}, \
}

class Spectrum():

    """
    Master Spectrum class.
    """

    def __init__(self, *args, **kwargs):

        self.ismodel  = kwargs.get('ismodel',False)
        self.wunit    = kwargs.get('wunit',DEFAULT_WAVE_UNIT)
        self.funit    = kwargs.get('funit',DEFAULT_FLUX_UNIT)
        self.filename = kwargs.get('file','')
        self.filename = kwargs.get('filename',self.filename)
        self.d_type   = kwargs.get('type')

        self.wave  = np.array(kwargs.get('wave', []))
        self.flux  = np.array(kwargs.get('flux', []))
        self.error = np.array(kwargs.get('error', [0 for i in range(len(self.wave))]))
        self.model = np.array(kwargs.get('model', []))
        self.sky   = np.array(kwargs.get('sky', []))
        self.name  = kwargs.get('name', 'Spectrum')  
        self.param = np.array(kwargs.get('param', [])) 

        self.mean_flux = np.nanmean(self.flux)
        self.std_flux  = np.nanstd(self.flux)  

        # process arguments
        # option 1: a filename is given
        if len(args) > 0:
            if isinstance(args[0],str):
                self.filename = args[0]

        # option 2: a spectrum ID is given
            elif isinstance(args[0],int) or isinstance(args[0],float) or isinstance(args[0],np.int64) or isinstance(args[0],np.float64):
                self.idkey = int(args[0])
                try:
                    sdb = keySpectrum(self.idkey)
                    if isinstance(sdb,bool) == False:
                        self.filename = sdb['DATA_FILE'].iloc[0]
                except:
                    print('Warning: problem reading in spectral database')

        # option 3: arrays are given - interpret as wave, flux, and optionally noise
        if len(args) > 1:
            if (isinstance(args[0],list) or isinstance(args[0],np.ndarray)) and (isinstance(args[1],list) or isinstance(args[1],np.ndarray)):
                kwargs['wave'] = kwargs.get('wave',args[0])
                kwargs['flux'] = kwargs.get('flux',args[1])
        if len(args) > 2:
            if isinstance(args[2],list) or isinstance(args[2],np.ndarray):
                kwargs['noise'] = kwargs.get('noise',args[2])

        if len(kwargs.get('wave','')) > 0 and len(kwargs.get('flux','')) > 0:
            self.wave = kwargs['wave']
            self.flux = kwargs['flux']
            if len(kwargs.get('noise','')) > 0:
                self.noise = kwargs['noise']
            else:
                self.noise = np.zeros(len(self.wave))
        # some extras
            others = ['pixel','mask','flag','flags','model']
            for o in others:
                if len(kwargs.get(o,'')) > 0:
                    setattr(self,o,kwargs[o])

        # read in file (NOTE: this overrules passing wave, flux, noise arrays)
        elif self.filename != '':

        # set up parameters
            mkwargs = {}
            #mkwargs['instrument']=self.instrument
            #mkwargs['instrument_mode']=self.instrument_mode
            mkwargs['filename']=self.filename
            self.simplefilename = os.path.basename(self.filename)
#            self.file = self.filename
            self.name = kwargs.get('name',self.simplefilename)

            # folder is by default the current directory
            mkwargs['folder'] = kwargs.get('folder','./')
            mkwargs['wunit'] = self.wunit
            mkwargs['funit'] = self.funit


        # breakouts for specific instruments
            if (kwargs.get('APOGEE') == True or kwargs.get('apogee') == True or kwargs.get('instrument','SPEX-PRISM').upper() == 'APOGEE') and self.filename != '':
                rs = _readAPOGEE(self.filename,**kwargs)
                self.instrument = 'APOGEE'
                # for k in list(rs.keys()): setattr(self,k.lower(),rs[k])
                self.history.append('Spectrum successfully loaded')
        # create a copy to store as the original
                self.original = copy.deepcopy(self)

            elif (kwargs.get('BOSS',False) == True or kwargs.get('boss',False) == True or kwargs.get('eboss',False) == True or kwargs.get('EBOSS',False) == True or kwargs.get('instrument','SPEX-PRISM').upper() == 'BOSS' or kwargs.get('instrument','SPEX-PRISM').upper() == 'EBOSS') and self.filename != '':
                rs = _readBOSS(self.filename)
                # for k in list(rs.keys()): setattr(self,k.lower(),rs[k])
                self.wunit = kwargs.get('wunit',u.Angstrom)
                self.funit = kwargs.get('funit',u.erg/(u.cm**2 * u.s * u.Angstrom))
                self.history.append('Spectrum successfully loaded')
        # create a copy to store as the original
                self.original = copy.deepcopy(self)

            else:
                rs = readSpectrum(self.filename,**mkwargs)
            for k in list(rs.keys()): setattr(self,k.lower(),rs[k])

        # information on model
        if self.ismodel == True:
            self.teff = kwargs.get('teff',np.nan)
            self.logg = kwargs.get('logg',np.nan)
            self.z = kwargs.get('z',np.nan)
            self.fsed = kwargs.get('fsed',np.nan)
            self.cld = kwargs.get('cld',np.nan)
            self.kzz = kwargs.get('kzz',np.nan)
            self.slit = kwargs.get('slit',np.nan)
            self.model = kwargs.get('model','')
            #mset = checkSpectralModelName(self.model) # dirty fix
            mset = 'btsettl08'
            if mset != False: 
                self.model = mset
                for k in list(SPECTRAL_MODELS[mset].keys()):
                    setattr(self,k.lower(),SPECTRAL_MODELS[mset][k])
            self.name = self.model+' model'
            self.shortname = self.name
            self.fscale = 'Surface'
            self.published = 'Y'

        else:
            #Save spectrum values into the spectrum object class
            self.wave   = kwargs.get('wave', [])
            self.flux   = kwargs.get('flux', [])
            self.sky    = kwargs.get('sky', [])
            self.sigmas = kwargs.get('sigmas', [0 for i in range(len(self.wave))])
            self.noise  = self.sigmas
            self.ivar   = kwargs.get('ivar', [])
            self.model  = kwargs.get('model', [])
            self.name   = kwargs.get('name', 'spectrum')  
            self.params = kwargs.get('params', [])   
            self.vsini  = kwargs.get('vsini', []) 
            self.type   = 'input'   
      

    def mask(self, **kwargs):

        """
        Mask all pixels that are out of the specified sigma cutoff range specified.

        Input: 'sigma'        : [lower cuttoff, upper cutoff]
               'pixel_buffer' : [lower mask pixel buffer, upper mask pixel buffer]
        """

        sigma = kwargs.get('sigma', [-2,1])
        pixel_buffer = kwargs.get('pixel_buffer', [0,0])

        #Find outlier flux and bad pixels 
        cut_low  = np.where(self.flux < (self.mean_flux + sigma[0]*self.std_flux))[0]
        cut_high = np.where(self.flux > (self.mean_flux + sigma[1]*self.std_flux))[0]

        group_low_cut = []
        group_high_cut = []

        for k, g in itertools.groupby(enumerate(cut_low), lambda x:x[0]-x[1]):
            group_low_cut.append(list(map(itemgetter(1), g)))

        for k, g in itertools.groupby(enumerate(cut_high), lambda x:x[0]-x[1]):
            group_high_cut.append(list(map(itemgetter(1), g)))

        cuts = []
        for pixels in group_low_cut:
            try: # prevent issues at the edges
                pix1, pixn = pixels[0], pixels[-1]
                for b in range(pixel_buffer[0]):
                    pixels.append(pix1 - (b+1))
                    pixels.append(pixn + (b+1))
                cuts.append(pixels)
            except:
                pass

        for pixels in group_high_cut:
            try: # prevent issues at the edges
                pix1, pixn = pixels[0], pixels[-1]
                for b in range(pixel_buffer[1]):
                    pixels.append(pix1 - (b+1))
                    pixels.append(pixn + (b+1))
                cuts.append(pixels)
            except:
                pass

        cuts = list(itertools.chain(*cuts))
        self.flux[cuts] = np.nan


    def plot(self, **kwargs):

        # Plot specifications
        xrange = kwargs.get('xrange', [self.wave[0], self.wave[-1]])
        yrange = kwargs.get('yrange', [self.mean_flux - 3*self.std_flux, self.mean_flux + 3*self.std_flux])
        items  = kwargs.get('items', ['spec'])
        style  = kwargs.get('style')
        title  = kwargs.get('title')
        save   = kwargs.get('save', False)
        output = kwargs.get('output', str(self.name) + '.pdf')
        
        # Shift by radial velocity
        rv      = kwargs.get('rv', 0)
        rv_wave = ap.rvShift(self.wave, rv=rv)
        
        # Set up plot
        fig = plt.figure(figsize=kwargs.get('figsize', [16,4]))                                                               
        ax  = fig.add_subplot(1,1,1) 

        # Plot masked spectrum
        if ('spectrum' in items) or ('spec' in items):
            if style == 'step':
                plt.step(rv_wave, self.flux, color='k', alpha=.8, linewidth=1, label=self.name, where='mid')
            else:
                plt.plot(rv_wave, self.flux, color='k', alpha=.8, linewidth=1, label=self.name)

        # Plot spectrum error
        if 'error' in items:
            if style == 'step':
                plt.step(self.wave, self.error, color='c', linewidth=1, alpha=.6, label='Error', where='mid')
            else:
                plt.plot(self.wave, self.error, color='c', linewidth=1, alpha=.6, label='Error')
        
        # Plot read in model
        if 'model' in items:
            if style == 'step':
                plt.step(self.wave, self.model, color='r', alpha=.8, linewidth=1, label='Model', where='mid')
            else:
                plt.plot(self.wave, self.model, color='r', alpha=.8, linewidth=1, label='Model')

        # Input additional spectrum objects to plot
        if 'objects' in kwargs:
            objects = kwargs.get('objects')
            obj_style = kwargs.get('obj_style')

            for obj in objects:
                if obj_style == 'step':
                    plt.step(obj.wave, obj.flux, label=obj.name)
                else:
                    plt.plot(obj.wave, obj.flux, label=obj.name)


        # Highlight specified bands
        highlight = kwargs.get('highlight')
        hcolor    = kwargs.get('hcolor', 'b')

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

                        if kwargs.get('line_style', 'short') == 'short':
                            plot_ypos_min = (ypos - yrange[0] -.15)/(yrange[1] - yrange[0])
                            plot_ypos_max = (ypos - yrange[0] -.1)/(yrange[1] - yrange[0])

                            plt.axvline(x=feature, ymin=plot_ypos_min, ymax=plot_ypos_max, linewidth=1, color = 'g')
                        
                        else:
                            plt.axvline(x=feature, linewidth=.3, color = 'g')

                        plt.text(feature, ypos-.2, lines, rotation=90, ha='center', color='b', fontsize=8)

        # Plot multiple dictionaries of lines in different colors
        if 'line_lists' in kwargs:
            line_lists = kwargs.get('line_lists')
            list_labels = kwargs.get('list_labels', ['list'+str(i) for i in range(len(line_lists))])

            line_colors = kwargs.get('line_colors', ['r', 'g', 'b', 'v', 'c', 'y'])
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

                            if kwargs.get('line_style', 'short') == 'short':
                                plot_ypos_min = (ypos - yrange[0] -.15)/(yrange[1] - yrange[0])
                                plot_ypos_max = (ypos - yrange[0] -.1)/(yrange[1] - yrange[0])

                                plt.axvline(x=feature, ymin=plot_ypos_min, ymax=plot_ypos_max, linewidth=1, color=line_colors[cindex])
                            
                            else:
                                plt.axvline(x=feature, linewidth=.3, color=line_colors[cindex])

                            plt.text(feature, ypos-.2, lines, rotation=90, ha='center', color='k', fontsize=8)
                cindex += 1

        # Plot plain vertical lines
        if 'vert_lines' in kwargs:
            vert_lines = kwargs.get('vert_lines')
            vcolor = kwargs.get('vcolor', 'g')

            for lines in vert_lines:
                lines = np.array(lines)
                range_indices = np.where((lines > xrange[0]) & (lines < xrange[1]))[0]
                lines = np.array(lines)[range_indices]

                ysize = yrange[1] - yrange[0]

                for line in lines:
                    plt.axvline(x=line, linewidth=.3, color=vcolor)
                    plt.text(line, yrange[0] + .03*ysize, round(line,2), rotation=90, ha='right', va='bottom', color='k', fontsize=10)

        if 'sigma_levels' in kwargs:
            sigma_levels = kwargs.get('sigma_levels')

            for level in sigma_levels:
                plt.axhline(y=self.mean_flux + level*self.std_flux, linestyle='--', c=np.random.rand(3,), label=r'$%s \sigma$'%(str(level)), linewidth=1)

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


def readSpectrum(verbose=False,*args,**kwargs):
    '''
    .. DOCS: will come back to this one
    '''
    # keyword parameters
    folder = kwargs.get('folder','')
    catchSN = kwargs.get('catchSN',True)
    kwargs['model'] = False
    instrument = kwargs.get('instrument','')
    #inst = checkInstrument(instrument)
    inst = 'APOGEE-RAW'
    if inst != False: instrument = inst
    # local = kwargs.get('local',True)

    # filename
    file = kwargs.get('file','')
    file = kwargs.get('filename',file)
    if (len(args) > 0):
        file = args[0]
    kwargs['filename'] = file
    kwargs['model'] = False

    # a filename must be passed
    if file == '':
        raise NameError('\nNo filename passed to readSpectrum')

    # first pass: check if file is local
    # if online == False:
    if not os.path.exists(os.path.normpath(file)):
        nfile = folder+os.path.basename(kwargs['filename'])
        if not os.path.exists(os.path.normpath(nfile)):
            if verbose==True: print('Cannot find '+kwargs['filename']+' locally, trying online\n\n')
        else:
            file=nfile
            kwargs['filename'] = file

    # instrument specific reads
    if instrument.upper()=='APOGEE': output = _readAPOGEE(file,**kwargs)
    # elif instrument.upper()=='BOSS': output = _readBOSS(file,**kwargs)
    # elif instrument.upper()=='LDSS3': output = _readIRAF(file,**kwargs)
    # elif instrument.upper()=='FIRE': output = _readFIRE(file,**kwargs)
    # elif instrument.upper()=='MAGE': output = _readMAGE(file,**kwargs)
    # elif instrument.upper()=='WFC3': output = _readWFC3(file,**kwargs)
    # elif instrument.upper()=='KAST-RED' or instrument.upper()=='KAST-BLUE': output = _readKAST(file,**kwargs)

    #   determine which type of file
    ftype = file.split('.')[-1]
    zipflag = ''

    # gzip compressed file - unzip and then rezip later
    if ftype == 'gz':
        zipflag = ftype
        file = file.replace('.'+ftype,'')
        with open(os.path.normpath(file), 'wb') as f_out, gzip.open(os.path.normpath(file+'.gz'), 'rb') as f_in:
            shutil.copyfileobj(f_in, f_out)

    # bz2 compressed file - unzip and then rezip later
    if ftype == 'bz2':
        zipflag = ftype
        file = file.replace('.'+ftype,'')
        with open(os.path.normpath(file), 'wb') as f_out, bz2.open(os.path.normpath(file+'.bz2'), 'rb') as f_in:
            shutil.copyfileobj(f_in, f_out)

    # fits file
    if (ftype == 'fit' or ftype == 'fits'):
   # df = fits.open(file)
   # with fits.open(file, ignore_missing_end=True) as data:
        with fits.open(os.path.normpath(file),ignore_missing_end=True) as data:
            data.verify('silentfix+warn')
            if 'NAXIS3' in list(data[0].header.keys()):
                d = np.copy(data[0].data[0,:,:])
            else:
                d =  np.copy(data[0].data)
            header = data[0].header

    # ascii file
    else:
        try:
            d = np.genfromtxt(os.path.normpath(file), comments='#', unpack=False, \
                missing_values = ('NaN','nan'), filling_values = (np.nan)).transpose()
        except ValueError:
            d = np.genfromtxt(os.path.normpath(file), comments=';', unpack=False, \
                 missing_values = ('NaN','nan'), filling_values = (np.nan)).transpose()
        header = fits.Header()      # blank header

    # delete file if this was an online read
   # if online and not local and os.path.exists(os.path.basename(file)):
   #     os.remove(os.path.normpath(os.path.basename(file)))

    # remove file if this was a zipped file
    if zipflag != '':
        os.remove(os.path.normpath(file))

    # assign arrays to wave, flux, noise
    if len(d[:,0]) > len(d[0,:]): d = d.transpose()  # array is oriented wrong

    # SDSS format for wavelength scale - in header and log format
    if kwargs.get('sdss',False) == True or (kwargs.get('waveheader',False) == True and kwargs.get('wavelog',False) == True):
        flux = d[0,:]
        if 'CRVAL1' in list(data[0].header.keys()) and 'CDELT1' in list(data[0].header.keys()):
            wave = 10.**(np.linspace(float(data[0].header['CRVAL1']),float(data[0].header['CRVAL1'])+len(flux)*float(data[0].header['CDELT1']),num=len(flux)))
        else: raise ValueError('\nCannot find CRVAL1 and CDELT1 keywords in header of fits file {}'.format(file))
        if len(d[:,0]) > 1:
            noise = d[1,:]
        else:
            noise = np.zeros(len(flux))
            noise[:] = np.nan

    #  wavelength scale in header and linear format
    elif (kwargs.get('waveheader',False) == True and kwargs.get('wavelinear',False) == True):
        flux = d[0,:]
        if 'CRVAL1' in list(data[0].header.keys()) and 'CDELT1' in list(data[0].header.keys()):
            wave = np.linspace(float(data[0].header['CRVAL1']),float(data[0].header['CRVAL1'])+len(flux)*float(data[0].header['CDELT1']),num=len(flux))
        else: raise ValueError('\nCannot find CRVAL1 and CDELT1 keywords in header of fits file {}'.format(file))
        if len(d[:,0]) > 1:
            noise = d[1,:]
        else:
            noise = np.zeros(len(flux))
            noise[:] = np.nan

    # wavelength is explicitly in data array 
    else:
        wave = d[0,:]
        flux = d[1,:]
        if len(d[:,0]) > 2:
            noise = d[2,:]
        else:
            noise = np.zeros(len(flux))
            noise[:] = np.nan

    output = {'wave': wave, 'flux': flux, 'noise': noise, 'header': header}


    # make sure arrays are numpy arrays
    output['wave']  = np.array(output['wave'])
    output['flux']  = np.array(output['flux'])
    output['noise'] = np.array(output['noise'])

    # fix places where noise is claimed to be 0
    w = np.where(output['noise'] == 0.)
    output['noise'][w] = np.nan

# fix nans in flux
#    w = np.where(np.isnan(flux) == True)
#    flux[w] = 0.

    # remove all parts of spectrum that are nans
    w = np.where(np.logical_and(np.isnan(output['wave']) == False,np.isnan(output['flux']) == False))
    output['wave']  = output['wave'][w]
    output['flux']  = output['flux'][w]
    output['noise'] = output['noise'][w]


    # fix to catch badly formatted files where noise column is S/N
    if catchSN == True:
          w = np.where(output['flux'] > np.nanmedian(output['flux']))
          if (np.nanmedian(output['flux'][w]/output['noise'][w]) < 1.):
              output['noise'] = output['flux']/output['noise']
              w = np.where(np.isnan(output['noise']))
              output['noise'][w] = np.nanmedian(output['noise'])

# add in instrument specific information
#    if inst != False:
#        for k in list(INSTRUMENTS[inst].keys()): output[k] = INSTRUMENTS[inst][k]  

# clean up
#    if url != '' and not local:
#        os.remove(os.path.basename(TMPFILENAME))
    if 'wunit' not in list(output.keys()): output['wunit'] = kwargs.get('wunit',DEFAULT_WAVE_UNIT)
    if 'funit' not in list(output.keys()): output['funit'] = kwargs.get('funit',DEFAULT_FLUX_UNIT)
    return output


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

            # Correct wavelength by -80 km/s 
            # wave = ap.rvShift(wave, rv=80)

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