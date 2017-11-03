import numpy as np
import os
import glob
import apogee_tools as ap
from astropy.io import fits
import matplotlib.pyplot as plt


def get_urls_from_header(path, frame_id):

    """
    This functions reads the header of a visit file
    and returns 1D-visit spectrum urls for different frames
    """
    
    #read visit spectrum header if local, obtain mjd and extract the "FRAME[N]" header keyword
    header=fits.open(path)[0].header
    
    frame_ids=None #frame ids are 1, 2, 3, ...
    frames=[] #strings : FRAM1, FRAME2, ...
    
    #option to obtain all the frames
    if frame_id =='all':
        #get all 12 frames
        frame_ids=['FRAME'+str(i) for i in np.arange(12)]

    #option to obtain one frame
    if isinstance(frame_id, int) and (frame_id !='all'):
        frame_ids=['FRAME'+str(frame_id)]
        
    #option to obtain specific frames
    if isinstance(frame_id, list):
        fram_ids=['FRAME'+str(i) for i in frame_id]
    
    #get frame ids, ignore key errors 
    for frame_id in frame_ids:
        try:
            frame=str(header[frame_id])
            #- add a "0" to the front if it is 7 digits = NUM
            if len(frame) ==7: frame='0'+frame
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
    
    nvisit: should be an integer, default is 'all' or a list 
    frame_id: should be an integer, default is  'all'
    """

    #download visit spectrum
    ap.download(apogee_id, type='apvisit',  dir=os.environ['APOGEE_DATA']+'/apVisit_data/') #apogee_tools should not download if the file exists
    
    nvisit=kwargs.get('nvisit','all')
    frame_id=kwargs.get('frame_id', 'all')
    
    #initialize paths
    vis_spec_paths=None
    
    #if the user wants to obtain all the visit spectra
    if nvisit =='all':
        nvisit =''
        
    #if the user wants to obtain one visit spectrum
    if isinstance(nvisit, str) or isinstance(nvisit, int) and (nvisit !='all'):
        p=os.environ['APOGEE_DATA']+'/ap1d_data/ap1d-'+apogee_id+'-'+str(nvisit)+'*'+'.fits'
        vis_spec_paths=glob.glob(p)
    
    #if the user wants to obtain a specific list of visit spectra
    if isinstance(nvisit, list):
        ps=[os.environ['APOGEE_DATA']+'/ap1d_data/ap1d-'+apogee_id+'-'+str(vis)+'*'+'.fits'
            for vis in nvisit]
        vis_spec_paths=np.concatenate([glob.glob(p) for p in ps])
    
    #get formatted urls by reading each specific visit header file
    url_list=[get_urls_from_header(path, frame_id) for path in vis_spec_paths]
    
    #some formatting 
    visit_numbers=[x.split('.fits')[0].split('-')[-1] for x in vis_spec_paths]
    urls= {'visit'+visit: url for (visit, url) in zip(visit_numbers , url_list)}
    
    return urls


def downloadAp1d(url, destination):
    os.system("wget {} -P {}".format(str(url), str(destination)))
    return 

