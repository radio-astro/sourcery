#!/usr/bin/env python

import os
from setuptools import setup
import Sourcery

setup(name="sourcery",
    version=Sourcery.__version__,
    description=".",
    author="Lerato Sebokolodi",
    author_email="Lerato Sebokolodi <mll.sebokolodi@gmail.com>",
    url="https://github.com/sphemakh/sourcery",
    packages=["Sourcery"],
    requires=["numpy", "matplotlib", "scipy", "astlib", "pyfits", "tigger"],
    scripts=["Sourcery/bin/" + i for i in os.listdir("Sourcery/bin")], 
    licence="This program should come with the GNU General Public Licence. "\
            "If not, find it at http://www.gnu.org/licenses/old-licenses/gpl-2.0-standalone.html",
    classifiers=[],
     )
