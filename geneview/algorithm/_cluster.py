"""Functions to cluster data.

Author: Shujia Huang
Date: 2021-04-30 11:50:26
"""
import warnings
import numpy as np
import pandas as pd

try:
    from scipy.cluster import hierarchy
    _no_scipy = False
except ImportError:
    _no_scipy = True

from ..utils import deprecate_positional_args


class _Dendrogram(object):
    """Agglomerative hierarchical clustering by scipy.cluster.hierarchy."""

    def __init__(self, data, linkage, method, metric, axis):
        """Calculating the relationships between the columns of data.

        Parameters
        ----------
        data : pandas.DataFrame
            Rectangular data

        axis : int, optional
            Which axis to use to calculate linkage. 0 is rows, 1 is columns.
        """
        if axis == 1:
            data = data.T

        if isinstance(data, pd.DataFrame):
            array = data.values
        else:
            array = np.asarray(data)
            data = pd.DataFrame(array)

        self.array = array
        self.data = data
        self.axis = axis

        self.shape = self.data.shape
        self.metric = metric
        self.method = method

        if linkage is None:
            # call linkage function.
            self.linkage = self.calculated_linkage
        else:
            self.linkage = linkage

        self.dendrogram = self.calculate_dendrogram()
        # self.dependent_coord = self.dendrogram["dcoord"]
        # self.independent_coord = self.dendrogram["icoord"]

    def _calculate_linkage_scipy(self):
        linkage = hierarchy.linkage(self.array,
                                    method=self.method,
                                    metric=self.metric)
        return linkage

    def _calculate_linkage_fastcluster(self):
        import fastcluster
        # Fastcluster has a memory-saving vectorized version, but only
        # with certain linkage methods, and mostly with euclidean metric
        # vector_methods = ("single", "centroid", "median", "ward")
        euclidean_methods = ("centroid", "median", "ward")
        euclidean = self.metric == "euclidean" and self.method in euclidean_methods
        if euclidean or self.method == "single":
            return fastcluster.linkage_vector(self.array,
                                              method=self.method,
                                              metric=self.metric)
        else:
            linkage = fastcluster.linkage(self.array, method=self.method,
                                          metric=self.metric)
            return linkage

    @property
    def calculated_linkage(self):

        try:
            return self._calculate_linkage_fastcluster()
        except ImportError:
            if np.product(self.shape) >= 10000:
                msg = ("Clustering large matrix with scipy. Installing "
                       "`fastcluster` may give better performance.")
                warnings.warn(msg)

        return self._calculate_linkage_scipy()

    @property
    def reordered_index(self):
        """Indices of the matrix, reordered by the dendrogram"""
        return self.dendrogram["leaves"]

    def calculate_dendrogram(self):
        """Calculates a dendrogram based on the linkage matrix
        Made a separate function, not a property because don't want to
        recalculate the dendrogram every time it is accessed.
        Returns
        -------
        dendrogram : dict
            Dendrogram dictionary as returned by scipy.cluster.hierarchy
            .dendrogram. The important key-value pairing is
            "reordered_ind" which indicates the re-ordering of the matrix
        """
        return hierarchy.dendrogram(self.linkage, no_plot=True, color_threshold=-np.inf)


@deprecate_positional_args
def hierarchical_cluster(
        data=None, linkage=None, method="average", metric="euclidean", axis=1
):
    """A dendrogram cluster (or call hierarchical cluster method) to calculate
    the relationships within a matrix.

    Parameters
    ----------
    data : pandas.DaraFrame.
        Rectangular data
    linkage : numpy.array, optional
        Linkage matrix which has been calculate by other cluster method.
    method : str, optional
        Linkage method to use. Anything valid for
        scipy.cluster.hierarchy.linkage
    metric : str, optional
        Distance metric. Anything valid for scipy.spatial.distance.pdist
    axis : int, optional
        Which axis to use to calculate cluster. 0 is rows, 1 is columns.
        default: 1.

    Returns
    -------
    dendropgram : _Dendrogram
        A dendrogram object.

    """
    if _no_scipy:
        raise RuntimeError("hierarchical cluster requires scipy to be installed.")

    return _Dendrogram(data=data, linkage=linkage, method=method, metric=metric, axis=axis)

