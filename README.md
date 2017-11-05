# apogee_tools

Forward modeling pipeline for [APOGEE](http://www.sdss.org/dr13/irspec/) VLM star spectra.

Authors:
* Christian Aganze (UCSD)
* Jessica Birky (UCSD)
* Adam Burgasser, PI (UCSD)
* Dino Chih-Chun Hsu (UCSD)
* Elizabeth Moreno (Guanajuato)
* Chris Theissen (BU)

This code also borrows from two other sources, see:
* [Starfish](https://github.com/iancze/Starfish) - Ian Czekala
* [apogee](https://github.com/jobovy/apogee) - Jo Bovy

## Code Setup

Dependencies:
* astropy
* astroquery
* numpy 
* pandas
* scipy
* matplotlib
* PyAstronomy

Create a new directory to store your data files:

	$ mkdir apogee_data/

Then download the APOGEE data info file for DR14:

	$ cd apogee_data/

	$ https://data.sdss.org/sas/dr14/apogee/spectro/redux/r8/allStar-l31c.2.fits --no-check-certificate
	$ https://data.sdss.org/sas/dr14/apogee/spectro/redux/r8/allVisit-l31c.2.fits --no-check-certificate

Download the apogee_tools code and then set up the following environmental variables in your `.bash_profile` or `.bashrc`:

	export PATH=$PATH:'/Users/path_to/apogee_tools'
	export APOGEE_DATA=/Users/path_to/apogee_data

## Downloading and reading APOGEE data files

To download APOGEE spectrum by 2MASS name and data type 'aspcap', or 'apstar':

	import apogee_tools as ap
	ap.download('2M03425325+2326495', type='aspcap')
	ap.download('2M03425325+2326495', type='apstar')

For data type 'apvisit' or 'ap1d': 

	ap.download('2M03425325+2326495', type='apvisit')
	ap.download('2M03425325+2326495', type='ap1d', visit=1, frame=1)

note: `type='apvisit'` will download the spectra for all visits observed, while `type='ap1d'` will download only the visit specified (and if not specified, will default to visit=1, frame=1).

For information on APOGEE data files, see the following:
* [aspcap](https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/LOCATION_ID/aspcapStar.html) - combined, continuum normalized spectra
* [apStar](https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/TELESCOPE/LOCATION_ID/apStar.html) - combined spectra
* [apVisit](https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/TELESCOPE/PLATE_ID/MJD5/apVisit.html) - individual raw visit spectra with telluric correction
* [ap1D](https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/red/MJD5/ap1D.html) - individual raw visit spectra with NO telluric correction

Also for info about the allStar file (such as aspcap pipeline parameters and photometry for all of the sources), see: [allStar](https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/allStar.html)

Once the data for a source has been downloaded, read aspcap or apStar files by specifying the 2MASS name and data type:

	data = ap.Spectrum(id='2M03425325+2326495', type='aspcap')

Or for single visit spectrum, indicate the index of the visit number at the end:

	data = ap.Spectrum(id='2M03425325+2326495', type='apvisit', visit=1)


## Basic Tools

**Search APOGEE catalog**

Example search--will search the `allStar-l30e.2.fits` you downloaded:

	params = ['TEFF', 'LOGG', 'M_H']
	ranges = [[-10000,4000], [0,5], [-2,2]]
	source_table = ap.multiParamSearch(par=params, select=ranges, save_dir='/path_to/')

**Read in a model grid**

Read in a model, specifying the parameters `[Teff, logg, [Fe/H]]`, grid type (listed below), and wavelength range `xrange`. Models sampled to APOGEE resolution are contained in the libraries folder of this package, and span the following parameter ranges: `PHOENIX: [[2500, 5500], [0.0, 5.5], [-1.0, 1.0]]`, `BTSETTL (CIFIST 2011b & 2015): [[2200, 3200], [2.5, 5.5], [-0.5, 0.0]]`. To use grids outside of these ranges, download the libraries from the links below, create an `.hdf5` file using Starfish, and add it to the `libraries` folder.

	mdl = ap.getModel(params=[3200, 5.0, 0.0], grid='BTSETTLb', xrange=[15200,16940])

Grid types:
* [PHOENIX](http://phoenix.astro.physik.uni-goettingen.de/) (Husser et. al. 2013)
<!-- * [BTSETTL](https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/) (Allard et. al. 2010) - CIFIST 2015 -->
* [BTSETTLb](https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011b/)  (Allard et. al. 2010) - CIFIST 2011b

**Plot data**

Some plotting examples:

	data = ap.Spectrum(id='2M03290406+3117075', type='aspcap')

	# plot spectrum
	data.plot()

	# plot aspcap model and noise:
	data.plot(items=['spec', 'apModel', 'noise'], save=True)

	# plot indentified lines (from Souto 2016):
	data.plot(items=['spec', 'lines'], xrange=[15200,15500], yrange=[.6,1.2])

Compare two spectra; return `chi` (chi-squared value between data and mdl), `norm_data` (`data` spectrum normalized), and `scaled_mdl` (`mdl` which has been scaled to `data`):

	chi, norm_data, scaled_mdl = ap.compareSpectra(data, mdl)


## Modeling Tools

**Synthesize a model**

First specify a dictionary of stellar parameters:

	params = {'teff': 3051, 'logg': 5.2, 'z': -0.25, 'vsini': 10., 'rv': -12, 'alpha': 0.2}

Read in some data you are creating a model for:

	ap.download('2M01195227+8409327', type='ap1d', visit=1, frame=1)
	data = ap.Spectrum(id='2M01195227+8409327', type='ap1d', visit=1)

Look up the spectrum's fiber number:

	ap_id, plates, mjds, fibers = ap.searchVisits(id_name='2M01195227+8409327')

Synthesize a model:

	mdl = ap.makeModel(params=params, fiber=fibers[0], plot=True, xrange=[15678,15694])

