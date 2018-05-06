Instrument Tools
================

APOGEE Data
-----------

Downloading and reading APOGEE data files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To download APOGEE spectrum by 2MASS name and data type ``aspcap``, or ``apstar``:

.. code-block:: python

	import apogee_tools as ap
	ap.download('2M03425325+2326495', type='aspcap')
	ap.download('2M03425325+2326495', type='apstar')

For data type 'apvisit' or 'ap1d': 

.. code-block:: python

	ap.download('2M03425325+2326495', type='apvisit')
	ap.download('2M03425325+2326495', type='ap1d', visit=1, frame=1)

Bote: ``type='apvisit'`` will download the spectra for all visits observed, while ``type='ap1d'`` will download only the visit specified (and if not specified, will default to visit=1, frame=1).

For information on APOGEE data files, see the following:

* `aspcap <https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/LOCATION_ID/aspcapStar.html>`_ - combined, continuum normalized spectra
* `apStar <https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/TELESCOPE/LOCATION_ID/apStar.html>`_ - combined spectra
* `apVisit <https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/TELESCOPE/PLATE_ID/MJD5/apVisit.html>`_ - individual raw visit spectra with telluric correction
* `ap1D <https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/red/MJD5/ap1D.html>`_ - individual raw visit spectra with NO telluric correction

Also for info about the allStar file (such as aspcap pipeline parameters and photometry for all of the sources), see: `allStar <https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/ASPCAP_VERS/RESULTS_VERS/allStar.html>`_.

Once the data for a source has been downloaded, read aspcap or apStar files by specifying the 2MASS name and data type:

.. code-block:: python

	data = ap.Spectrum(id='2M03425325+2326495', type='aspcap')

Or for single visit spectrum, indicate the index of the visit number at the end:

.. code-block:: python

	data = ap.Spectrum(id='2M03425325+2326495', type='apvisit', visit=1)


Search the APOGEE catalog
~~~~~~~~~~~~~~~~~~~~~~~~~

Example search--will search the ``allStar-l30e.2.fits`` you downloaded:

.. code-block:: python

	params = ['TEFF', 'LOGG', 'M_H']
	ranges = [[-10000,4000], [0,5], [-2,2]]
	source_table = ap.multiParamSearch(par=params, select=ranges, dir='/path_to/')

Look up aspcap parameters in ``allStar-l30e.2.fits`` for specific list of 2MASS IDs:

.. code-block:: python

	tm_ids = ['2M01195227+8409327']
	ap_dict = ap.returnAspcapTable(tm_ids, params=['TEFF', 'LOGG', 'M_H', 'SNR'], save=False)


Plot data
~~~~~~~~~

Some plotting examples:

.. code-block:: python

	data = ap.Spectrum(id='2M03290406+3117075', type='aspcap')

	# plot spectrum
	data.plot()

	# plot aspcap model and noise:
	data.plot(items=['spec', 'apModel', 'noise'], save=True)

	# plot indentified lines (from Souto 2016):
	data.plot(items=['spec', 'lines'], xrange=[15200,15500], yrange=[.6,1.2])

Mask outlying flux
~~~~~~~~~~~~~~~~~~

Specify number of standard deviations above and below the mean of the flux to cut (``sigma = [lower cuttoff, upper cutoff]``), and the number pixels to buffer each side of the cut (``pixel_buffer = [lower mask pixel buffer, upper mask pixel buffer]``):

.. code-block:: python

	data.mask(sigma=[3,2], pixel_buffer=[0,3])

Chi-squared comparison
~~~~~~~~~~~~~~~~~~~~~~

Compare two spectra; return ``chi`` (chi-squared value between data and mdl), ``norm_data`` (``data`` spectrum normalized), and ``scaled_mdl`` (``mdl`` which has been scaled to ``data``):

.. code-block:: python

	chi, norm_data, scaled_mdl = ap.compareSpectra(data, mdl)


NIRSPEC Data
------------

More info coming soon.


Adding New Instruments
----------------------

More info coming soon.