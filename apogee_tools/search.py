import os
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.io import fits, ascii
from astroquery.simbad import Simbad

db_path = os.environ['APOGEE_DATA'] 
AP_PATH = db_path

data_releases = {'dr10':'allStar-v304.fits', 'dr11':'allStar-v402.fits', \
    'dr12':'allStar-v603.fits', 'dr13':'allStar-l30e.2.fits', 'dr14':'allStar-l31c.2.fits'} 

def searchStars(**kwargs):

    """
    Searching DR12-DR13 Databases
    https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/allStar.html
    @Christian Aganze

    Input:  'id_name' : 2MASS name
            'rel'     : data release to search 
    (note: I realized the latest data release is a compilation of all previous data releases, so DR10-DR12 are unnecessary to search)

    Output: 'apogee_id', 'location_id' : needed to download spectra
    """

    searchval = kwargs.get('id_name')
    release   = kwargs.get('rel', 'dr14')

    #Read the database file
    hdu        = fits.open(AP_PATH + '/' + data_releases[release])
    keys       = hdu[1].header
    search_par = kwargs.get('SearchPar', 'APOGEE_ID')
    t          = hdu[1].data[search_par]
    condition  = np.where(t==searchval)[0]
    ap_id      = hdu[1].data['APOGEE_ID'][condition]
    asp_id     = hdu[1].data['ASPCAP_ID'][condition]
    loc_id     = hdu[1].data['LOCATION_ID'][condition]

    print(ap_id, loc_id)

    return ap_id, loc_id[0]


def searchVisits(**kwargs):

    """
    Search visits data file allVisit-l30e.2.fits, which should be stored in the folder set as environmental variable 'APOGEE_DATA'
    @Jessica Birky

    Input:  'id_name' : 2MASS name

    Output: 'ap_id', 'plates', 'mjds', 'fibers' : info needed to download visit spectra
    """

    ap_id = kwargs.get('id_name')

    #Read the database file
    hdu = fits.open(AP_PATH + '/allVisit-l30e.2.fits')

    keys = hdu[1].header
    data = hdu[1].data

    search = data['APOGEE_ID']
    pos = np.where(search == ap_id)

    plates = data['PLATE'][pos]
    mjds = data['MJD'][pos]
    fibers = data['FIBERID'][pos]

    return ap_id, plates, mjds, fibers

    
def download(star_id, **kwargs):

    """
    Download star by 2MASS name
    @Jessica Birky

    Input:  'star_id' : 2MASS name
            'type'    : aspcap, apstar, or apvisit
    """

    dr = kwargs.get('rel', 'dr14')

    # Datatypes to download
    d_type = kwargs.get('type','aspcap').lower()

    # Get download directory; create if it doesn't already exist
    default_dir = AP_PATH + '/%s_data/' %(d_type)
    dl_dir = kwargs.get('dir', default_dir)

    if not os.path.exists(dl_dir):
        os.makedirs(dl_dir)

    key = {'dr13': ['r6','l30e','l30e.2'], 'dr14': ['r8','l31c','l31c.2']}

    if d_type == 'aspcap':

        fname = 'aspcapStar-%s-%s-%s.fits'%(key[dr][0], key[dr][2], star_id)

        if fname not in os.listdir(dl_dir):

            #Look up location id from allStar-l30e.2.fits file
            ap_id, loc_id = searchStars(id_name=star_id, rel=dr)

            if len(ap_id) != 0:  

                print('Downloading %s from %s'%(ap_id, dr))
                os.system("wget https://data.sdss.org/sas/%s/apogee/spectro/redux/%s/stars/%s/%s/%s/%s -P %s" %(dr, key[dr][0], key[dr][1], key[dr][2], loc_id, fname, dl_dir)) 

            else:
                print(star_id + ' does not exist in ' + dr)
        else:
            print('Already have file for ' + star_id)

    if d_type == 'apstar': 

        #Look up location id from allStar-l30e.2.fits file
        ap_id, loc_id = searchStars(id_name=star_id)

        fname = 'apStar-r6-'+star_id+'.fits'
        
        if len(ap_id) != 0:
            if fname not in os.listdir(dl_dir):
                os.system("wget https://data.sdss.org/sas/%s/apogee/spectro/redux/%s/stars/apo25m/%s/%s -P %s" %(dr, key[dr][0], loc_id, fname, dl_dir))
        else:
            print(star_id + ' does not exist in ' + dr)

    if d_type == 'apvisit': #Can only download DR13 now

        #This will download the data in two formats: 
            #(1) An hdf5 file which stores all visits
            #(2) The individual fits files for every visit

        ap_id, plates, mjds, fibers = searchVisits(id_name=star_id)
        nVisits = len(plates)

        if len(ap_id) != 0:
            #Download individual visit .fits files
            for v in range(nVisits):
                #Look up plate, mjd and fiber numbers from allVisit-l30e.2.fits file
                plate, mjd, fiber = plates[v], mjds[v], fibers[v]

                fname = 'apVisit-r6-%s-%s-%s.fits' %(plate, mjd, fiber)

                if fname not in os.listdir(dl_dir):
                    os.system("wget -O %sapVisit-%s-%s.fits https://data.sdss.org/sas/%s/apogee/spectro/redux/%s/apo25m/%s/%s/%s" %(dr, key[dr][0], dl_dir, ap_id, v, plate, mjd, fname))
        else:
            print(star_id + ' does not exist in ' + dr)


def multiParamSearch(**kwargs):

    """
    Search allStar file for stars in specific ranges for a parameter (or multiple parameters)
    @Jessica Birky

    Input:  'par'    : names of parameters to filter by (ex: 'TEFF', 'LOGG', 'M_H'...)
            'select' : range of values for each parameter (must be in same order as 'par')
            'dir'    : directory where you want to save search table

    Output: data table of stars matching search criteria; printed and saved as csv
    """

    search_par = kwargs.get('par', ['TEFF'])
    select     = kwargs.get('select', [[3000, 4500]])
    save       = kwargs.get('save', True)
    release    = kwargs.get('rel', ['dr14'])

    # default save directory is in your APOGEE_DATA environmental variable folder
    if not os.path.exists(AP_PATH + '/tables'):
        os.makedirs(AP_PATH + '/tables')
    output_dir = kwargs.get('dir', AP_PATH + '/tables')

    data_releases = {'dr10':'allStar-v304.fits', 'dr11':'allStar-v402.fits', \
    'dr12':'allStar-v603.fits', 'dr13':'allStar-l30e.2.fits'}

    database = db_path + '/allStar-l30e.2.fits'
    data = Table(fits.open(database)[1].data)

    for i in range(len(search_par)):
        p = data[search_par[i]]
        if isinstance(select[i][0], type('k')):
            # matching by strings
            condition = [j for j in range(0, len(data[search_par[i]]))  if data[search_par[i]][j] in select[i]]

        else:
            print(" ALL stars with ", search_par[i], "between", select[i])
            condition = np.where((p>=select[i][0] ) & (p<=select[i][1]))[0]
        data = data[condition]

    # Put only interesting values into the data table: 2MASS name and aspcap parameters
    data_dict  = {'2MASS_ID':data['APOGEE_ID'], 'TEFF':data['TEFF'], 'LOGG':data['LOGG'], 'M_H':data['M_H']}
    # data_dict = {'2MASS_ID':data['APOGEE_ID'], 'RA':data['RA'], 'DEC':data['DEC']}
    data_table = pd.DataFrame(data=data_dict)

    # Save data frame to csv file and save to the 'output' directory specified by keyword argument
    tmin, tmax = select[0][0], select[0][1]
    lmin, lmax = select[1][0], select[1][1]
    mmin, mmax = select[2][0], select[2][1]

    fname = 'teff_%s_%s_logg_%s_%s_feh_%s_%s.csv'%(tmin, tmax, lmin, lmax, mmin, mmax)
    if save == True:
        save_dir = '%s/%s' %(output_dir, fname)
        data_table.to_csv(save_dir)
        print('Table saved to ', save_dir)

    return data


def returnSimbadParams(id_list):

    """
    Enter list of 2MASS IDs, and return stellar parameters from simbad
    @Jessica Birky

    Input:  'id_list' : list of 2MASS IDs

    Output: 'result' : table of stellar parameters obtained from simbad
    """

    customSimbad = Simbad()
    customSimbad.add_votable_fields('fe_h', 'rv_value', 'rvz_error', 'sptype', 'rot')
    customSimbad.get_votable_fields()

    objects = []
    no_entry = 0
    
    for s in id_list:
        ID = str(s)
        obj = customSimbad.query_object(ID)

        try:
            spt = str(obj['SP_TYPE']).split('-\n')[1]
            d = {'ID':str(s), 'SPT':spt, 'TEFF':obj['Fe_H_Teff'], 'LOGG':obj['Fe_H_log_g'], \
                 'FE_H':obj['Fe_H_Fe_H']} 
            df = pd.DataFrame(data=d)
            objects.append(df)
        except:
            no_entry += 1

    print('No simbad values for ' + no_entry + ' objects.')

    result = pd.concat(objects)

    return result




