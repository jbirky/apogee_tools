import numpy as np
import matplotlib.pyplot as plt
import os, sys
import glob
from astropy.io import fits
from astropy.utils.data import download_file
import apogee_tools as ap


def get_urls_from_header(path, frame_id):

    """
    This functions reads the header of a visit file
    and returns 1D-visit spectrum urls for different frames
    @Christian Aganze
    """

    #read visit spectrum header if local, obtain mjd and extract the "FRAME[N]" header keyword
    header = fits.open(path)[0].header
    
    frame_ids = None #frame ids are 1, 2, 3, ...
    frames = [] #strings : FRAM1, FRAME2, ...
    
    #option to obtain all the frames
    if frame_id =='all':
        #get all 12 frames
        frame_ids = ['FRAME'+str(i) for i in np.arange(12)]

    #option to obtain one frame
    if isinstance(frame_id, int) and (frame_id !='all'):
        frame_ids = ['FRAME'+str(frame_id)]
        
    #option to obtain specific frames
    if isinstance(frame_id, list):
        fram_ids = ['FRAME'+str(i) for i in frame_id]
    
    #get frame ids, ignore key errors 
    for frame_id in frame_ids:
        try:
            frame = str(header[frame_id])
            #- add a "0" to the front if it is 7 digits = NUM
            if len(frame) == 7: frame='0'+frame
            frames.append(frame)
        except KeyError:
            continue 
    
    #get mjd
    mjd=header['MJD5']
    
    #format urls
    urls=[['https://data.sdss.org/sas/dr14/apogee/spectro/redux/r8/red/'+ str(mjd)+\
    '/ap1D-'+str(band)+'-'+str(frame)+'.fits' for band in ['a', 'b', 'c']] for frame in frames]
    
    return dict(zip(frame_ids, urls))


def get_1dspec_urls(apogee_id, **kwargs):

    """
    This function retrieves 1D visit spectra urls given an apogee id,
    a visit number (i.e, first visit, second visit etc.) and a frame id
    @Christian Aganze
    
    nvisit: should be an integer, default is 'all' or a list 
    frame_id: should be an integer, default is  'all'
    """

    #download visit spectrum
    print('Downloading apVisit files for {}'.format(apogee_id))
    ap_path = kwargs.get('dir', os.environ['APOGEE_DATA'])

    ap.download(apogee_id, type='apvisit',  dir=ap_path+'/apvisit_data/') #apogee_tools should not download if the file exists
    
    nvisit=kwargs.get('nvisit','all')
    frame_id=kwargs.get('frame_id', 'all')
    
    #initialize paths
    vis_spec_paths=None
    
    #if the user wants to obtain all the visit spectra
    if nvisit =='all':
        nvisit =''
        
    #if the user wants to obtain one visit spectrum
    if isinstance(nvisit, str) or isinstance(nvisit, int) and (nvisit !='all'):
        p=os.environ['APOGEE_DATA']+'/apvisit_data/apVisit-'+apogee_id+'-'+str(nvisit)+'*'+'.fits'
        vis_spec_paths=glob.glob(p)
    
    #if the user wants to obtain a specific list of visit spectra
    if isinstance(nvisit, list):
        ps=[os.environ['APOGEE_DATA']+'/apvisit_data/apVisit-'+apogee_id+'-'+str(vis)+'*'+'.fits'
            for vis in nvisit]
        vis_spec_paths=np.concatenate([glob.glob(p) for p in ps])
    
    #get formatted urls by reading each specific visit header file
    url_list=[get_urls_from_header(path, frame_id) for path in vis_spec_paths]
    
    #some formatting 
    visit_numbers=[x.split('.fits')[0].split('-')[-1] for x in vis_spec_paths]
    urls= {'visit'+visit: url for (visit, url) in zip(visit_numbers , url_list)}
    
    return urls


def coadd_spectra(spectra, errorArray):

	"""
	@Chris Theissen
	"""

	# Do a masked, uncertainty weighted average over all the spectra
	AveragedFluxes = np.ma.average(np.array(spectra), axis=0, weights = 1./np.array(errorArray)**2.)
	AveragedErrors = np.ma.sqrt(1. / np.sum(1./np.array(errorArray)**2., axis=0))

	return AveragedFluxes, AveragedErrors


def coadd_epoch(star_id, epoch=1, dr='dr14', show_progess=False):

	"""
	@Chris Theissen
	"""

	# Get all APOGEE MJDs for each visit for the object
	ap_id, plates, mjds, fibers = ap.searchVisits(id_name=star_id)

	# Get the list of urls for each spectrum from each visit
	spectralist = ap.utils.ap1d.get_1dspec_urls(star_id)

	# Check if visit (epoch) is valid
	if visit > len(mjds): 
		raise Exception('Epoch %s is invalid. Only %s epochs available for %s.'%(visit, len(mjds), star_id))

	# Which fiber is it?
	fiber = fibers[visit-1]

	# Initialize empty arrays for all the spectra and errors
	spectraA = []
	errorsA  = []
	spectraB = []
	errorsB  = []
	spectraC = []
	errorsC  = []

	# Get the spectra from each from of the epoch
	for j in spectralist['visit%s'%visit]:

		# Download the 1-D spectra
		fitsA = fits.open(download_file(spectralist['visit%s'%visit]['%s'%j][0], cache=False, show_progress=show_progress))
		fitsB = fits.open(download_file(spectralist['visit%s'%visit]['%s'%j][1], cache=False, show_progress=show_progress))
		fitsC = fits.open(download_file(spectralist['visit%s'%visit]['%s'%j][2], cache=False, show_progress=show_progress))
		
		# Chip A
		t1     = fitsA[1].data[fiber] 
		wave1  = fitsA[4].data[fiber]
		order1 = np.unravel_index(np.argsort(wave1, axis=None), wave1.shape)
		flux1  = t1#*1e-17
		err1   = fitsA[2].data[fiber][order1]
		mask1  = fitsA[3].data[fiber][order1]
		wave1  = wave1[order1]

		fluxA = np.ma.array(flux1, mask=mask1)
		errA  = np.ma.array(err1,  mask=mask1)

		spectraA.append(fluxA)
		errorsA.append(errA)

		# Chip B
		t2     = fitsB[1].data[fiber] 
		wave2  = fitsB[4].data[fiber]
		order2 = np.unravel_index(np.argsort(wave2, axis=None), wave2.shape)
		flux2  = t2#*1e-17
		err2   = fitsB[2].data[fiber][order2]
		mask2  = fitsB[3].data[fiber][order2]
		wave2  = wave2[order2]

		fluxB = np.ma.array(flux2, mask=mask2)
		errB  = np.ma.array(err2,  mask=mask2)

		spectraB.append(fluxB)
		errorsB.append(errB)

		# Chip C
		t3     = fitsC[1].data[fiber] 
		wave3  = fitsC[4].data[fiber]
		order3 = np.unravel_index(np.argsort(wave3, axis=None), wave3.shape)
		flux3  = t3#*1e-17
		err3   = fitsC[2].data[fiber][order3]
		mask3  = fitsC[3].data[fiber][order3]
		wave3  = wave3[order3]

		fluxC = np.ma.array(flux3, mask=mask3)
		errC  = np.ma.array(err3,  mask=mask3)

		spectraC.append(fluxC)
		errorsC.append(errC)

	# Average the spectra for the epoch
	FinalFluxA, FinalErrorA = coadd_spectra(spectraA, errorsA)
	FinalFluxB, FinalErrorB = coadd_spectra(spectraB, errorsB)
	FinalFluxC, FinalErrorC = coadd_spectra(spectraC, errorsC)

	# return (for each band) flux, error, and wavelengths
	return FinalFluxA, FinalErrorA, wave1, FinalFluxB, FinalErrorB, wave2, FinalFluxC, FinalErrorC, wave3

