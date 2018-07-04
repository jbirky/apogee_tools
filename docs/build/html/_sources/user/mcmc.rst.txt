MCMC Fitting 
============

.. warning:: These functions are under construction. 

Setup
-----

1. Copy the ``config.yaml`` and ``run.py`` from the main directory to an external folder.

2. Edit your configuration script ``config.yaml``, which should look something like below.

3. In your new directory run ``python run.py`` in terminal.


Configuration
-------------

.. literalinclude:: ../../../config.yaml


Pre-MCMC Testing
----------------

To test to make sure all of the modeling modules are working, run the following command in terminal::

	$ python run.py make_model

which should return something like::

	[25.732014894485474s] MCMC initialization step complete.

	##################################################
	Making model: teff=3500 logg=4.5 fe_h=0.0 rv=-4.77 vsini=5.79 alpha=1.0

	[0.07615256309509277s] Interpolated model
	[0.0025053024291992188s] Shifted radial velocity
	[0.0032796859741210938s] Applied vsini broadening
	[0.05470013618469238s] Convolved telluric model
	[0.08379793167114258s] Applied LSF broadening 


.. image:: images/make_model.png

To test by eye, that your initial MCMC parameters are some close to the data::

	$ python run.py test_fit


Running the MCMC
----------------

Run the MCMC::

	$ python run.py mcmc

Plot the outputs::

	$ python run.py walkers
	$ python run.py corner
