"""Build and validate a geneview release for PyPI.

This script builds both the source distribution (sdist) and the wheel for the
version declared in ``setup.py``, then validates the artifacts with
``twine check``. Uploading is intentionally left as a separate, explicit step
so a release is never published by accident (a PyPI version can never be
re-used once uploaded).

Usage::

    # 1. Build + validate the current version's artifacts
    python release.py

    # 2. Inspect the printed command, then upload manually when ready
    twine upload dist/geneview-<version>.tar.gz dist/geneview-<version>-*.whl

Author: Shujia Huang
Date: 2021-04-30 (modernized 2026)
"""
import importlib.util
import sys
from subprocess import check_call

# Load the ``meta`` namespace from setup.py to learn the distribution
# name and the version currently being released.
spec = importlib.util.spec_from_file_location("_", "./setup.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

distname = module.meta.__DISTNAME__
version = module.meta.__VERSION__

# Build both sdist and wheel with the PEP 517 front-end.
check_call([sys.executable, "-m", "build"])

# Validate only the artifacts for the version being released so stale files
# from earlier builds in dist/ do not cause false failures.
sdist = f"dist/{distname}-{version}.tar.gz"
wheel = f"dist/{distname}-{version}-py3-none-any.whl"
check_call(["twine", "check", sdist, wheel])

print(f"\nBuilt and validated {distname} {version}.")
print("To publish to PyPI, run:\n")
print(f"    twine upload {sdist} {wheel}\n")
