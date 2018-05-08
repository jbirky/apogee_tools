Modelling Tools
===============

Reading in Model Grids
----------------------

Read in a model, specifying the parameters ``[Teff, logg, [Fe/H]]``, grid type (listed below), and wavelength range ``xrange``. Models sampled to APOGEE resolution are contained in the libraries folder of this package, and span the following parameter ranges: ``PHOENIX: [[2500, 5500], [0.0, 5.5], [-1.0, 1.0]]``, ``BTSETTL (CIFIST 2011b & 2015): [[2200, 3200], [2.5, 5.5], [-0.5, 0.0]]``. To use grids outside of these ranges, download the libraries from the links below, create an ``.hdf5`` file using Starfish, and add it to the ``libraries`` folder.

.. code-block:: python

	mdl = ap.getModel(params=[3200, 5.0, 0.0], grid='BTSETTL', xrange=[15200,16940])

Grid types:

* `PHOENIX <http://phoenix.astro.physik.uni-goettingen.de/>`_ (Husser et. al. 2013)
* `BTSETTL <https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011b/>`_ (Allard et. al. 2010) - CIFIST 2011b


Synthesize a model
------------------

First specify a dictionary of stellar parameters:

.. code-block:: python

	params = {'teff': 3051, 'logg': 5.2, 'z': -0.25, 'vsini': 10., 'rv': -12, 'alpha': 0.2}

Read in some data you are creating a model for:

.. code-block:: python

	ap.download('2M01195227+8409327', type='ap1d', visit=1, frame=1)
	data = ap.Apogee(id='2M01195227+8409327', type='ap1d', visit=1)

Look up the spectrum's fiber number:

.. code-block:: python

	ap_id, plates, mjds, fibers = ap.searchVisits(id_name='2M01195227+8409327')

Synthesize a model: (with resolution options: 23k, 50k, and 300k)

.. code-block:: python

	mdl = ap.makeModel(params=params, fiber=fibers[0], plot=True, xrange=[15678,15694], res='300k')