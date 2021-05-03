# Copyright (C) 2016-2021 Shujia Huang <huangshujia9@gmail.com>
# temporarily redirect config directory to prevent matplotlib importing
# testing that for writeable directory which results in sandbox error in
# certain easy_install versions
import os
from argparse import Namespace

DESCRIPTION = "Geneview: A python package for genomics data visualization."
LONG_DESCRIPTION = """\
Geneview is a library for making attractive and informative genomics graphics in Python. 
It is built on top of matplotlib and tightly integrated with the PyData stack, including 
support for numpy and pandas data structures.

Some of the features that geneview offers are: 

- Functions for visualizing univariate and bivariate distributions or for comparing them between subsets of data.
- Tools that fit and visualize linear regression models for different kinds of independent and dependent variables.
- High-level abstractions for structuring grids of plots that let you easily build complex visualizations.
"""

meta = Namespace(
    __DISTNAME__     = "geneview",
    __AUTHOR__       = "Shujia Huang",
    __AUTHOR_EMAIL__ = "huangshujia9@gmail.com",
    __URL__          = "https://github.com/ShujiaHuang/geneview",
    __LICENSE__      = "BSD (3-clause)",
    __DOWNLOAD_URL__ = "https://github.com/ShujiaHuang/geneview",
    __VERSION__      = "0.0.7",
)

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup, find_packages


if __name__ == "__main__":

    setup(name=meta.__DISTNAME__,
          author=meta.__AUTHOR__,
          author_email=meta.__AUTHOR_EMAIL__,
          maintainer=meta.__AUTHOR__,
          maintainer_email=meta.__AUTHOR_EMAIL__,
          description=DESCRIPTION,
          long_description=LONG_DESCRIPTION,
          license=meta.__LICENSE__,
          url=meta.__URL__,
          version=meta.__VERSION__,
          download_url=meta.__DOWNLOAD_URL__,
          install_requires=[
              "numpy",
              "scipy",
              "pandas",
              "matplotlib",
              "seaborn",
          ],
          packages=find_packages(),
          classifiers=[
             "Intended Audience :: Science/Research",
             "Programming Language :: Python :: 3.7",
             "Programming Language :: Python :: 3.8",
             "License :: OSI Approved :: BSD License",
             "Topic :: Scientific/Engineering :: Bio-Informatics",
             "Topic :: Scientific/Engineering :: Visualization",
             "Topic :: Multimedia :: Graphics",
             "Operating System :: POSIX",
             "Operating System :: Unix",
             "Operating System :: MacOS"],
          )
