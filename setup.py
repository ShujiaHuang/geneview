#!/usr/bin/env python
# Copyright (C) 2016 Shujia Huang <huangshujia9@gmail.com>
import os
# temporarily redirect config directory to prevent matplotlib importing
# testing that for writeable directory which results in sandbox error in
# certain easy_install versions
os.environ["MPLCONFIGDIR"] = "."

DESCRIPTION = "Geneview: genomics data visualization"
LONG_DESCRIPTION = """\
Geneview is a library for making attractive and informative genomics graphics in Python. It is built on top of matplotlib and tightly integrated with the PyData stack, including support for numpy data structures.

Some of the features that geneview offers are

- Functions for visualizing univariate and bivariate distributions or for comparing them between subsets of data

- Tools that fit and visualize linear regression models for different kinds of independent and dependent variables

- High-level abstractions for structuring grids of plots that let you easily build complex visualizations
"""

DISTNAME = 'geneview'
MAINTAINER = 'Shujia Huang'
MAINTAINER_EMAIL = 'huangshujia9@gmail.com'
URL = 'https://github.com/ShujiaHuang/geneview'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/ShujiaHuang/geneview'
VERSION = '0.0.1.1'


try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup, find_packages

def check_dependencies():
    install_requires = []

    # Just make sure dependencies exist, I haven't rigorously
    # tested what the minimal versions that will work are
    # (help on that would be awesome)
    try:
        import numpy
    except ImportError:
        install_requires.append('numpy')
    #try:
    #    import scipy
    #except ImportError:
    #    install_requires.append('scipy')
    try:
        import matplotlib
    except ImportError:
        install_requires.append('matplotlib')

    return install_requires

if __name__ == "__main__":

    install_requires = check_dependencies()

    setup(name=DISTNAME,
          author=MAINTAINER,
          author_email=MAINTAINER_EMAIL,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          long_description=LONG_DESCRIPTION,
          license=LICENSE,
          url=URL,
          version=VERSION,
          download_url=DOWNLOAD_URL,
          install_requires=install_requires,
          #packages=['geneview', 'geneview.tests'],
          packages=find_packages(),
          classifiers=[
             'Intended Audience :: Science/Research',
             'Programming Language :: Python :: 2.7',
             'Programming Language :: Python :: 3.3',
             'Programming Language :: Python :: 3.4',
             'License :: OSI Approved :: BSD License',
             'Topic :: Scientific/Engineering :: Visualization',
             'Topic :: Multimedia :: Graphics',
             'Operating System :: POSIX',
             'Operating System :: Unix',
             'Operating System :: MacOS'],
          )
