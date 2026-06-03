# Installation

<cite>
**Referenced Files in This Document**
- [setup.py](file://setup.py)
- [requirements.txt](file://requirements.txt)
- [README.md](file://README.md)
- [geneview/__init__.py](file://geneview/__init__.py)
- [_cluster.py](file://geneview/algorithm/_cluster.py)
- [release.py](file://release.py)
</cite>

## Table of Contents
1. [Introduction](#introduction)
2. [System Requirements](#system-requirements)
3. [Install from PyPI (Recommended)](#install-from-pypi-recommended)
4. [Install from Source](#install-from-source)
5. [Install Development Version](#install-development-version)
6. [Optional Dependencies](#optional-dependencies)
7. [Dependency Management](#dependency-management)
8. [Verification Steps](#verification-steps)
9. [Troubleshooting Common Issues](#troubleshooting-common-issues)

## Introduction
This section provides step-by-step installation instructions and guidance for managing dependencies for GeneView. It covers the recommended installation via PyPI, installing from source, and installing development versions. It also documents system requirements, optional dependencies, and troubleshooting tips to help you resolve common installation issues.

## System Requirements
- Python version: Python 3.7, 3.8, or 3.9
- Operating systems: POSIX-compliant systems including Unix and macOS

These compatibility requirements are declared in the package metadata.

**Section sources**
- [setup.py:52-63](file://setup.py#L52-L63)

## Install from PyPI (Recommended)
Install the latest released version of GeneView using pip. This command installs GeneView along with all required dependencies.

- Command: pip install geneview

What this installs:
- numpy
- scipy
- pandas
- matplotlib
- seaborn

Notes:
- The README explicitly states that installing the released version via pip will install all dependencies.
- The setup script declares the same set of dependencies in install_requires.

**Section sources**
- [README.md:18-27](file://README.md#L18-L27)
- [setup.py:44-50](file://setup.py#L44-L50)
- [requirements.txt:1-6](file://requirements.txt#L1-L6)

## Install from Source
You can install GeneView directly from the source repository. This method is useful if you want to build the package locally or inspect the source before installation.

Typical steps:
- Clone or download the repository
- Navigate to the repository root
- Build the distribution
- Install the built package

Reference:
- The release script demonstrates building a source distribution and uploading to PyPI, indicating the standard Python packaging workflow.

**Section sources**
- [release.py:19-21](file://release.py#L19-L21)

## Install Development Version
To install the development version from the repository, you can use pip to install directly from the repository URL. This allows you to work with the latest changes before they are officially released.

Typical steps:
- Use pip to install from the repository URL
- Optionally specify a branch or commit hash for a specific development state

Reference:
- The release script shows how distributions are prepared and uploaded, which implies a standard development workflow using pip and setup.py.

**Section sources**
- [release.py:13-21](file://release.py#L13-L21)

## Optional Dependencies
While the core installation includes numpy, scipy, pandas, matplotlib, and seaborn, GeneView can optionally use additional libraries to enhance functionality:

- statsmodels: Used by certain plotting functions for advanced kernel density estimation capabilities. If present, GeneView prefers statsmodels for greater flexibility; otherwise, it falls back to scipy.
- fastcluster: Improves performance for hierarchical clustering on large matrices. If available, GeneView will use fastcluster; otherwise, it will use scipy.

Benefits:
- statsmodels enables advanced KDE features and improves accuracy for certain statistical visualizations.
- fastcluster significantly speeds up hierarchical clustering computations for large datasets.

Detection and fallback:
- GeneView checks for optional dependencies at runtime and gracefully falls back to scipy when optional packages are not installed.

**Section sources**
- [README.md:324-334](file://README.md#L324-L334)
- [_cluster.py:66-92](file://geneview/algorithm/_cluster.py#L66-L92)

## Dependency Management
- Required dependencies are declared in setup.py and mirrored in requirements.txt.
- At runtime, GeneView imports core modules such as matplotlib and exposes higher-level plotting functions. This ensures that the underlying visualization stack is initialized properly.

Key points:
- install_requires lists numpy, scipy, pandas, matplotlib, and seaborn.
- requirements.txt contains the same five packages.
- geneview/__init__.py imports matplotlib and sets default font rendering preferences.

**Section sources**
- [setup.py:44-50](file://setup.py#L44-L50)
- [requirements.txt:1-6](file://requirements.txt#L1-L6)
- [geneview/__init__.py:1-15](file://geneview/__init__.py#L1-L15)

## Verification Steps
After installation, verify that GeneView is ready to use:

- Import the package in Python and check that core plotting functions are available.
- Run a simple plotting example to ensure matplotlib and the GeneView plotting functions initialize correctly.
- Confirm that optional dependencies (statsmodels, fastcluster) are recognized if you have installed them.

Example verification commands (conceptual):
- Import geneview and call a basic plotting function
- Verify that no import errors occur for core modules (numpy, pandas, matplotlib, scipy, seaborn)

Note: The repository’s README includes quick-start examples demonstrating imports and basic plotting, which you can adapt for verification.

**Section sources**
- [README.md:28-74](file://README.md#L28-L74)

## Troubleshooting Common Issues
- Missing system dependencies for compiled extensions:
  - Some dependencies (e.g., scipy, matplotlib) may require system-level libraries during compilation. Ensure your system has the necessary build tools and libraries installed before attempting installation.
- Conflicts with existing packages:
  - If you encounter version conflicts, consider using a virtual environment to isolate GeneView and its dependencies.
- Optional dependency warnings:
  - If statsmodels or fastcluster are not installed, GeneView will still function using scipy. A warning may be issued when fastcluster is preferred for performance on large datasets.

Guidance:
- Use a virtual environment to avoid conflicts with system or user-installed packages.
- Reinstall with pip if dependency resolution fails, ensuring you have the latest pip version.

**Section sources**
- [_cluster.py:88-92](file://geneview/algorithm/_cluster.py#L88-L92)