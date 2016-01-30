"""
This module contains miscellaneous functions for ``geneview``.
"""
from __future__ import print_function, division

def chr_id_cmp(a, b):
    """
    Sorted the chromosome by the order.
    
    Parameters
    ----------
        a, b : string or int. 
            a and b are the chromosomes' id. They could be 'chr1', 'chr2'
            or '1' and '2'.
    
    Returns
    -------
    Must be one of the number in [-1, 0, 1] just as the same with the 
    python build-in function: ``cmp()``
    """
    a = a.lower().replace('_', '')
    b = b.lower().replace('_', '')
    achr = a[3:] if a.startswith('chr') else a
    bchr = b[3:] if b.startswith('chr') else b

    try:
        return cmp(int(achr), int(bchr))
    except ValueError:
        if achr.isdigit() and not bchr.isdigit(): return -1
        if bchr.isdigit() and not achr.isdigit(): return 1
        # X Y
        return cmp(achr, bchr) 


def is_numeric(s):
    """
    It's a useful function for checking if a data is a numeric.

    This function could identify any kinds of numeric: e.g. '0.1', '10', '-2.',
    2.5, etc

    Parameters
    ----------
    s : int, float or string.
        The input could be any kind of single value except the scalable
        type of python like 'list()', 'dict()', 'set()'

    Returns
    -------
        A boolean

    Notes
    -----
        http://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-float-in-python
    """
    try:
        float(s)
        return True
    except ValueError:
        return False
