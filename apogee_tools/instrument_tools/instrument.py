import numpy as np
from astropy.io import fits
from astropy import units as u
import os


def formatDesignation(input_id):

	# Make sure id is in the proper 2MASS format
	if '+' in input_id:
		spec_id = '2M' + input_id.split('+')[0][-8:] + '+' + input_id.split('+')[1]
	elif '-' in input_id:
		spec_id = '2M' + input_id.split('-')[0][-8:] + '-' + input_id.split('-')[1]
	else:
		print('Designation improperly formated. Must be form "__00034394+8606422".')

	return spec_id


def getShortname(input_id):

    """
    Return shortname of a spectrum designation. 
    Ex: >> getShortname('2M03425325+2326495')
        '2M0342+2326'

    """
    
    name = ap.formatDesignation(input_id)
    
    return name[0:6] + name[10:15]


class apogeeFile():

    """ 
    Stores entire APOGEE file info.
    """

    def __init__(self, **kwargs):  

        self.d_type = kwargs.get('type')

        # Get designation of 2MASS object
        input_id  = kwargs.get('id')            
        self.name = ap.formatDesignation(input_id)

        if self.d_type == 'aspcap':

            file_dr14 = '%s/aspcap_data/aspcapStar-r8-l31c.2-%s.fits' %(AP_PATH, spec_id)

            """ ASPCAP file info:
            HDU0: The Primary Header
            HDU1: Spectrum array
            HDU2: Error array
            HDU3: Best fit spectrum
            HDU4: ASPCAP data table
            Information found at https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/LOCATION_ID/aspcapStar.html 
            """

            if (os.path.exists(file_dr14) == False):
                ap.download(self.name, type='aspcap')

            openFile = fits.open(file_dr14)

            self.HDU0 = openFile[0]
            self.HDU1 = openFile[1]
            self.HDU2 = openFile[2]
            self.HDU3 = openFile[3]
            self.HDU4 = openFile[4]

            #conversion from pixel to wavelength, info available in the hdu header
            crval = self.HDU1.header['CRVAL1']
            cdelt = self.HDU1.header['CDELT1']
            wave  = np.array(pow(10, crval+cdelt*np.arange(self.HDU1.header['NAXIS1']))/10000)*u.micron #microns

            # convert fluxes from  (10^-17 erg/s/cm^2/Ang) to  ( erg/s/cm^2/Mircon)
            spectras = [1e-13*np.array(f)*u.erg/u.s/u.centimeter**2/u.micron for f in self.HDU1.data]
                    
            self.wave    = 10000*np.array(wave.value)
            self.flux    = np.array(self.HDU1.data)
            self.error   = np.array(self.HDU2.data)
            self.model	 = np.array(self.HDU3.data)

            #Obtain aspcap parameters
            self.params = self.HDU4.data['PARAM'] 

        elif self.d_type == 'apstar':

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

            if (os.path.exists(file_dr14) == False):
                ap.download(self.name, type='apstar')

            openFile = fits.open(file_dr14)

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
        
            #conversion from pixel to wavelength, info available in the hdu header
            crval = self.HDU1.header['CRVAL1']
            cdelt = self.HDU1.header['CDELT1']
            wave  = np.array(pow(10, crval+cdelt*np.arange(self.HDU1.header['NAXIS1']))/10000)*u.micron #microns
            
            # convert fluxes from  (10^-17 erg/s/cm^2/Ang) to  ( erg/s/cm^2/Mircon)
            spectras = [1e-13*np.array(f)*u.erg/u.s/u.centimeter**2/u.micron for f in self.HDU1.data]
                    
            self.wave    = 10000*np.array(wave.value)
            self.flux    = np.array(self.HDU1.data)
            self.error   = np.array(self.HDU2.data)
            self.sky     = np.array(self.HDU4.data)
            self.mask 	 = np.array(self.HDU3.data)

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

            openFile = fits.open(self.file)

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
            self.HDU10 = openFile[10]

            flux  = self.HDU1
            error = self.HDU2
            wave  = self.HDU4

            #Combine the data from the three chips into one list
            self.wave = np.array(list(wave[0]) + list(wave[1]) + list(wave[2]))
            self.flux = np.array(list(flux[0]) + list(flux[1]) + list(flux[2]))

            self.error = np.array(error.data)
        
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
            self.error = np.concatenate([err1, err2, err3])[::-1]

        