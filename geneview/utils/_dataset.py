"""
Utility functions for getting datasets from geneview online dataset repository.
"""
import os
import re
import pandas as pd
from urllib.request import urlopen, urlretrieve


def get_dataset_names():
    """Report available example datasets, useful for reporting issues."""

    url = "https://github.com/ShujiaHuang/geneview-data"
    with urlopen(url) as resp:
        html = resp.read()

    pat = r"/ShujiaHuang/geneview-data/blob/master/(\w*).csv"
    datasets = re.findall(pat, html.decode())
    return datasets


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
    Load the preview data of GWAS

        >>> from geneview.utils import load_dataset, get_dataset_names
        >>> names = get_dataset_names()  # Get the name list of dataset
        >>> gwas_data = load_dataset("gwas")  # load GWAS data
        >>> file_path = load_dataset("bcftools.vcf.stats")  # only return a path name
    """

    path_csv = "https://raw.githubusercontent.com/ShujiaHuang/geneview-data/master/{0}.csv"
    path_full = "https://raw.githubusercontent.com/ShujiaHuang/geneview-data/master/{0}"
    path_name = path_full.format(name) if "." in name else path_csv.format(name)

    if cache:
        cache_path = os.path.join(_get_data_home(data_home), os.path.basename(path_name))
        if not os.path.exists(cache_path):
            urlretrieve(path_name, cache_path)
        path_name = cache_path

    if path_name.endswith(".csv"):
        df = pd.read_csv(path_name, **kws)
        if df.iloc[-1].isnull().all():
            df = df.iloc[:-1]
        return df
    else:
        return path_name


def _get_data_home(data_home=None):
    """Return the path of the geneview data directory.

    This is used by the ``load_dataset`` function.

    If the ``data_home`` argument is not specified, the default location
    is ``~/geneview-data``.

    Alternatively, a different default location can be specified using the
    environment variable ``GENEVIEW_DATA``.
    """
    if data_home is None:
        data_home = os.environ.get('GENEVIEW_DATA', os.path.join('~', 'geneview-data'))
    data_home = os.path.expanduser(data_home)
    if not os.path.exists(data_home):
        os.makedirs(data_home)

    return data_home
