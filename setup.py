import os
import sys
from setuptools import setup

import apogee_tools

setup(
    name="apogee_tools",
    version="0.1",
    author="Jessica Birky",
    author_email="jbirky@ucsd.edu",
    url="http://apogee-tools.readthedocs.io",
    license="MIT",
    description=("Forward modeling framework for fitting atmospheric models to stellar spectra."),
    long_description=open("README.md").read(),
    install_requires=["numpy", "pandas", "astropy", "astroquery", "emcee", "matplotlib", "pandas", "PyAstronomy", "scipy"],
)