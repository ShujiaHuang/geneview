"""
Utility functions for getting datasets from geneview online dataset repository.
"""
from ..ext.six.moves.urllib.request import urlopen, urlretrieve

def get_dataset_names():
    """Report available example datasets, useful for reporting issues."""
    # delayed import to not demand bs4 unless this function is actually used
    from bs4 import BeautifulSoup
    http = urlopen('https://github.com/ShujiaHuang/geneview-data/')
    gh_list = BeautifulSoup(http)

    return [l.text.replace('.csv', '')
            for l in gh_list.find_all("a", {"class": "js-directory-link"})
            if l.text.endswith('.csv')]
