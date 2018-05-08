Analysis Tools
==============

Plotting Features
-----------------

.. todo:: More info coming soon.

Atomic and Molecular Lines
--------------------------

Search atomic and molecular lines from the following databases:

**NIST** (`download1 <https://physics.nist.gov/PhysRefData/ASD/lines_form.html>`_): atomic lines

**HITEMP** (`download2 <ftp://cfa-ftp.harvard.edu/pub/HITEMP-2010/)>`_, `paper2 <https://www.cfa.harvard.edu/atmosphere/publications/2010-HITEMP-JQSRT-111.pdf>`_): molecular lines H2O, CO2, CO, NO, and OH

**APOGEE** (`download3 <https://zenodo.org/record/32629#.Vi0XBBCrSfS>`_), `paper3 <https://arxiv.org/abs/1502.04080>`_): atomic and molecular lines

Usage
~~~~~

.. code-block:: python

	ap.listLibraries()

	ap.listSpecies('APOGEE_ATOMS')

Search a wavelegth region for certain species of atoms/molecules (returns a dictionary):

.. code-block:: python

	lines = apsearchLines(species=['OH', 'Fe I'], libraries=['NIST', 'APOGEE_ATOMS', 'APOGEE_MOLEC'], range=[15200,15300])