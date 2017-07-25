# apogee_tools

Forward modeling [APOGEE](http://www.sdss.org/dr13/irspec/) spectra.

Authors:
* Jessica Birky (UCSD)
* Adam Burgasser (UCSD)
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

Or if you don't have access, download the files for DR13:

	$ wget https://data.sdss.org/sas/dr13/apogee/spectro/redux/r6/allStar-l30e.2.fits --no-check-certificate
	$ wget https://data.sdss.org/sas/dr13/apogee/spectro/redux/r6/allVisit-l30e.2.fits --no-check-certificate

Download the apogee_tools code and then set up the following environmental variables in your `.bash_profile` or `.bashrc`:

	export PATH=$PATH:'/Users/path_to/apogee_tools'
	export APOGEE_DATA=/Users/path_to/apogee_data

## Downloading and reading APOGEE data files

To download APOGEE spectrum by 2MASS name and data file type (aspcap, apStar):

	import apogee_tools as ap
	ap.download('2M03290406+3117075', type='aspcap')

For information on APOGEE data files, see the following:
* [aspcap](https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/LOCATION_ID/aspcapStar.html) - combined, continuum normalized spectra
* [apStar](https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/TELESCOPE/LOCATION_ID/apStar.html) - combined spectra
* [apVisit](https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/TELESCOPE/PLATE_ID/MJD5/apVisit.html) - individual raw spectra

Also infomation about the allStar file (such as parameters and photometry for all of the sources), see: [allStar](https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/allStar.html)

Once the data for a source has been downloaded, read aspcap or apStar files by specifying the 2MASS name and data type:

	data = ap.Spectrum(id='2M03290406+3117075', type='aspcap')

<!-- Or for single visit spectrum, indicate the index of the visit number at the end:

	data = ap.Spectrum(id='2M03290406+3117075-0', type='apvisit') -->


## Tools

**Search APOGEE catalog**

Example search--will search the allStar-l30e.2.fits you downloaded:

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

<!-- Scale one spectrum to another, where `sp1` and `sp2` are spectrum objects:

	scale = ap.calcScale(sp1, sp2)

**Radial Velocity shift**

Perform velocity shift the wavelength of a spectrum:

	data = ap.Spectrum(id='2M03290406+3117075', type='aspcap')
	rv = -80 # Radial velocity in km/s
	rv_wave = ap._rvShift(data.wave, rv=rv)

	# Create spectrum object with rv shift
	spec = ap.Spectrum(wave=rv_shift, flux=data.flux, sigmas=data.sigmas, name=data.name)
	spec.plot()

Calculate radial velocity by cross correlating with a template:

	wave_rng = [15200,15700]

	# Read in data
	data = ap.Spectrum(id='2M03290406+3117075', type='aspcap')

	# Read in a model template to cross-correlate to
	mdl = ap.getModel(params=[3000, 5.0, 0.0], grid='BTSETTLb', xrange=wave_rng, subCont=True)

	# Return radial velocity
	rv, sp1, sp2 = ap.optimizeRV(data, mdl, xrange=wave_rng)

**Rotational Broadening**

Add rotational broadening to a model using the [PyAstronomy](http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/rotBroad.html) routine:

	mdl = ap.getModel(params=[3000, 5.0, 0.0], grid='BTSETTLb', xrange=wave_rng)
	rot_mdl = ap.smoothVSINI(mdl, vsini=15, xlim=[15200,15500], plot=True)
 -->

<!--  Remove large file from git:
 git filter-branch --force --index-filter 'git rm --cached -r --ignore-unmatch oops.iso' --prune-empty --tag-name-filter cat -- --all
rm -rf .git/refs/original/
git reflog expire --expire=now --all
git gc --prune=now
git gc --aggressive --prune=now -->

