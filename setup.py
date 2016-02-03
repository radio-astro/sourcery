#!/usr/bin/env python

import os
from setuptools import setup, find_packages
import Sourcery

setup(name="sourcery",
    version=Sourcery.__version__,
    description="Tools for creating high fidelity source catalogues from radio interferometric datasets",
    author="Lerato Sebokolodi",
    author_email="mll.sebokolodi@gmail.com",
    url="https://github.com/radio-astro/sourcery",
    packages=find_packages(),
    requires=["setuptools", "numpy", "matplotlib", "scipy", "astlib", "pyfits", "tigger"],
    scripts=["Sourcery/bin/" + i for i in os.listdir("Sourcery/bin")], 
    license="GPL2",
    classifiers=[],
 )
