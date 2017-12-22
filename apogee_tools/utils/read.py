import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import h5py
import apogee_tools as ap

#Get the path of apogee_tools file
FULL_PATH  = os.path.realpath(__file__)
BASE, NAME = os.path.split(FULL_PATH)

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