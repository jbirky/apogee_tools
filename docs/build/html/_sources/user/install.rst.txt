Installation
============

Dependencies:
-------------

* `astropy <http://www.astropy.org/>`_
* `astroquery <https://astroquery.readthedocs.io/en/latest/>`_
* `emcee <http://emcee.readthedocs.io/en/stable/user/install.html>`_
* `numpy <http://www.numpy.org/>`_
* `matplotlib <http://matplotlib.org/>`_
* `pandas <https://pandas.pydata.org/pandas-docs/stable/install.html>`_
* `PyAstronomy <https://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaCDoc/installingPyA.html>`_
* `scipy <https://www.scipy.org/install.html>`_


APOGEE Data Setup
-----------------

Create a new directory to store your data files:

	$ mkdir apogee_data/

Then download the APOGEE data info file for DR14:

	$ cd apogee_data/

	$ wget https://data.sdss.org/sas/dr14/apogee/spectro/redux/r8/allStar-l31c.2.fits --no-check-certificate
	$ wget https://data.sdss.org/sas/dr14/apogee/spectro/redux/r8/allVisit-l31c.2.fits --no-check-certificate

Download the apogee_tools code and then set up the following environmental variables in your ``.bash_profile`` or ``.bashrc``:

	export PATH=$PATH:'/Users/path_to/apogee_tools'
	export APOGEE_DATA=/Users/path_to/apogee_data


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

