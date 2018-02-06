Model Files:
------------

**cannon_phoenix:** 2nd order polynomial coefficients, trained on each pixel for interpolating phoenix model grids. Can be interpolated to produce a spectrum using ap.interpolateGrid(), with resolution options: 23k, 50k, 300k, and 500k. See [TheCannon](https://github.com/annayqho/TheCannon) for details on training.

**pwv_R300k_airmass1.0** and **pwv_R300k_airmass1.5**: High resolution (300k) telluric models. See the [source](ftp://ftp.eso.org/pub/dfs/pipelines/skytools/telluric_libs) and [publication](http://adsabs.harvard.edu/abs/2014A%26A...568A...9M).

**BTSETTLb_APOGEE.hdf5**: 23k resolution [BTSETTLb](https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011b/)  (Allard et. al. 2010) models (CIFIST 2011b). Ranges {Teff, logg, [Fe/H]} = {[2500, 5500], [0.0, 5.5], [-1.0, 1.0]}, and can be read in using ap.getModel().

**PHOENIX_APOGEE.hdf5**: 23k resolution [PHOENIX](http://phoenix.astro.physik.uni-goettingen.de/) (Husser et. al. 2013) models. Ranges {Teff, logg, [Fe/H]} = {[2200, 3200], [2.5, 5.5], [-0.5, 0.0]}, and can be read in using ap.getModel().

**lw_solartrans_apogee.csv**: low resolution telluric model, prepared by Adam.