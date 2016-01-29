"""
This module contains miscellaneous functions for ``geneview``.
"""
from __future__ import print_function, division

def is_numeric(s):
    """
    It's a useful function for checking if a string is a numeric.

    This function could identify any kinds of numeric. e.g. '0.1', '10', '-2.',
    2.5, etc

    Parameters
    ----------
    s : int, float or string.
        The input could be any kind of single value except the scalable
        type of python like 'list()', 'dict()', 'set()'

    Notes
    -----
        http://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-float-in-python
    """
    try:
        float(s)
        return True
    except ValueError:
        return False
