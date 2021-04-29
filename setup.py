# Copyright (C) 2016 Shujia Huang <huangshujia9@gmail.com>
import os
# temporarily redirect config directory to prevent matplotlib importing
# testing that for writeable directory which results in sandbox error in
# certain easy_install versions
#os.environ["MPLCONFIGDIR"] = "."

DESCRIPTION = "Geneview: A python package for genomics data visualization"
LONG_DESCRIPTION = """\
Geneview is a library for making attractive and informative genomics graphics in Python. It is built on top of matplotlib and tightly integrated with the PyData stack, including support for numpy data structures.

Some of the features that geneview offers are

- Functions for visualizing univariate and bivariate distributions or for comparing them between subsets of data

- Tools that fit and visualize linear regression models for different kinds of independent and dependent variables

- High-level abstractions for structuring grids of plots that let you easily build complex visualizations
"""

DISTNAME = "geneview"
MAINTAINER = "Shujia Huang"
MAINTAINER_EMAIL = "huangshujia9@gmail.com"
URL = "https://github.com/ShujiaHuang/geneview"
LICENSE = "BSD (3-clause)"
DOWNLOAD_URL = "https://github.com/ShujiaHuang/geneview"
VERSION = "0.0.5"


try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup, find_packages


if __name__ == "__main__":

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
