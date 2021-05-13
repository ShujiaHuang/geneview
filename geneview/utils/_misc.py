"""
This module contains miscellaneous functions for ``geneview``.
"""
import operator
import pandas as pd
import numpy as np
from scipy import stats


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

    Notes
    -----
        http://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-float-in-python
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def is_integer(s):
    """
    This function could identify any kinds of integer

    Parameters
    ----------
        s : int, float or string.
            The input could be any kind of single value except the scalable
            type of python like 'list()', 'dict()', 'set()'

    Returns
    -------
        A boolean. Ture if `s` is integer value else False 
    """
    if is_numeric(s):
        return True if '.' not in s else False
    else:
        return False


def iqr(a):
    """Calculate the IQR for an array of numbers."""
    a = np.asarray(a)
    q1 = stats.scoreatpercentile(a, 25)
    q3 = stats.scoreatpercentile(a, 75)
    return q3 - q1


def _kde_support(data, bw, gridsize, cut, clip):
    """Establish support for a kernel density estimate."""
    support_min = max(data.min() - bw * cut, clip[0])
    support_max = min(data.max() + bw * cut, clip[1])
    return np.linspace(support_min, support_max, gridsize)


def categorical_order(values, order=None):
    """Return a list of unique data values.

    Determine an ordered list of levels in ``values``.

    Parameters
    ----------
    values : list, array, Categorical, or Series of pandas
        Vector of "categorical" values
    order : list-like, optional
        Desired order of category levels to override the order determined
        from the ``values`` object.

    Returns
    -------
    order : list
        Ordered list of category levels not including null values.

    """
    if order is None:
        if hasattr(values, "categories"):
            order = values.categories
        else:
            try:
                order = values.cat.categories
            except (TypeError, AttributeError):
                try:
                    order = values.unique()
                except AttributeError:
                    order = pd.unique(values)
                try:
                    np.asarray(values).astype(np.float)
                    order = np.sort(order)
                except (ValueError, TypeError):
                    order = order
        order = filter(pd.notnull, order)
    return list(order)


def freedman_diaconis_bins(a):
    """Calculate number of hist bins using Freedman-Diaconis rule."""
    # From http://stats.stackexchange.com/questions/798/
    a = np.asarray(a)
    h = 2 * iqr(a) / (len(a) ** (1 / 3))
    # fall back to sqrt(a) bins if iqr is 0
    if h == 0:
        return int(np.sqrt(a.size))
    else:
        return int(np.ceil((a.max() - a.min()) / h))
