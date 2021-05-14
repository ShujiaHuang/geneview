"""
This module contains miscellaneous functions for ``geneview``.
"""


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
        A boolean. Ture if `s` is numeric else False

    Examples
    --------

    >>> from geneview.utils import is_numeric
    >>> is_numeric('a')
    False
    >>> is_numeric(1)
    True
    >>> is_numeric(3.14)
    True

    Notes
    -----
        http://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-float-in-python
    """
    try:
        float(s)
        return True
    except ValueError:
        return False
