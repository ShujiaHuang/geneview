"""
Utility functions for getting datasets from geneview online dataset repository.
"""
import os

import csv
import pandas as pd

try:
    from six.moves.urllib.request import urlopen, urlretrieve
except ImportError:
    from ..ext.six.moves.urllib.request import urlopen, urlretrieve

from ._misc import is_numeric

def get_dataset_names():
    """Report available example datasets, useful for reporting issues."""
    # delayed import to not demand bs4 unless this function is actually used
    from bs4 import BeautifulSoup
    http = urlopen('https://github.com/ShujiaHuang/geneview-data/')
    gh_list = BeautifulSoup(http)

    return [l.text.replace('.csv', '')
            for l in gh_list.find_all("a", {"class": "js-directory-link"})
            if l.text.endswith('.csv')]


def load_dataset(name, cache=True, data_home=None, **kws):
    """Load a dataset from the online repository (requires internet).

    Parameters
    ----------
    name : str
        Name of the dataset (`name`.csv on
        https://github.com/ShujiaHuang/geneview-data).  You can obtain list of
        available datasets using :func:`get_dataset_names`

    cache : boolean, optional
        If True, then cache data locally and use the cache on subsequent calls

    data_home : string, optional
        The directory in which to cache data. By default, uses ~/geneview_data/

    kws : dict, optional
        Passed to pandas.read_csv

    Examples
    --------
    Load the preview data of GOYA

        >>> import geneview as gv
        >>> goya_preview = gv.util.load_dataset('GOYA_preview')

    """
    path = "https://github.com/ShujiaHuang/geneview-data/raw/master/{0}.csv"
    full_path = path.format(name)

    if cache:
        cache_path = os.path.join(_get_data_home(data_home),
                                  os.path.basename(full_path))
        if not os.path.exists(cache_path):
            urlretrieve(full_path, cache_path)
        full_path = cache_path

    """
    data = []
    with open(full_path) as f:
        f_csv = csv.reader(f)
        headers = next(f_csv)
        # keep the first element as string in dataset(.csv format)
        data = [[row[0]] + map(_tr, row[1:]) for row in f_csv] 
    df = pd.DataFrame(data, columns=headers)
    """
    df = pd.read_csv(full_path, **kws)
    if df.iloc[-1].isnull().all():
        df = df.iloc[:-1]

    return df


def _tr(s):
    """ transform numeric to be float data. """
    return float(s) if is_numeric(s) else s


def _get_data_home(data_home=None):
    """Return the path of the geneview data directory.

    This is used by the ``load_dataset`` function.

    If the ``data_home`` argument is not specified, the default location
    is ``~/geneview-data``.

    Alternatively, a different default location can be specified using the
    environment variable ``GENEVIEW_DATA``.
    """
    if data_home is None:
        data_home = os.environ.get('GENEVIEW_DATA',
                                   os.path.join('~', 'geneview-data'))
    data_home = os.path.expanduser(data_home)
    if not os.path.exists(data_home):
        os.makedirs(data_home)

    return data_home
