"""Tests for geneview.algorithm module (hierarchical clustering).

Author: Shujia Huang
"""
import pytest
import numpy as np
import pandas as pd

from geneview.algorithm import hierarchical_cluster
from geneview.algorithm._cluster import _Dendrogram


class TestDendrogram:
    """Tests for the _Dendrogram class."""

    def test_basic_clustering_dataframe(self):
        """Test clustering on a simple DataFrame with row axis."""
        np.random.seed(42)
        df = pd.DataFrame(np.random.rand(10, 3))
        dendro = _Dendrogram(data=df, linkage=None, method="average",
                             metric="euclidean", axis=0)
        assert dendro.data is not None
        assert dendro.array.shape == (10, 3)
        assert dendro.shape == (10, 3)
        assert dendro.linkage is not None
        assert dendro.dendrogram is not None
        assert len(dendro.reordered_index) == 10

    def test_clustering_column_axis(self):
        """Test clustering along column axis (axis=1)."""
        np.random.seed(42)
        df = pd.DataFrame(np.random.rand(5, 8))
        dendro = _Dendrogram(data=df, linkage=None, method="average",
                             metric="euclidean", axis=1)
        # When axis=1, data is transposed internally
        assert dendro.reordered_index is not None
        assert len(dendro.reordered_index) == 8

    def test_clustering_numpy_array(self):
        """Test clustering on a raw numpy array."""
        np.random.seed(42)
        arr = np.random.rand(12, 4)
        dendro = _Dendrogram(data=arr, linkage=None, method="average",
                             metric="euclidean", axis=0)
        assert isinstance(dendro.data, pd.DataFrame)
        assert dendro.array.shape == (12, 4)
        assert len(dendro.reordered_index) == 12

    def test_reordered_index_is_permutation(self):
        """Verify reordered_index is a valid permutation of row indices."""
        np.random.seed(42)
        n = 20
        df = pd.DataFrame(np.random.rand(n, 5))
        dendro = _Dendrogram(data=df, linkage=None, method="average",
                             metric="euclidean", axis=0)
        reordered = dendro.reordered_index
        assert sorted(reordered) == list(range(n))

    def test_different_linkage_methods(self):
        """Test various linkage methods supported by scipy."""
        np.random.seed(42)
        df = pd.DataFrame(np.random.rand(15, 3))
        for method in ["single", "complete", "average", "weighted", "ward"]:
            metric = "euclidean" if method == "ward" else "euclidean"
            dendro = _Dendrogram(data=df, linkage=None, method=method,
                                 metric=metric, axis=0)
            assert len(dendro.reordered_index) == 15

    def test_precomputed_linkage(self):
        """Test passing precomputed linkage matrix."""
        from scipy.cluster.hierarchy import linkage
        np.random.seed(42)
        arr = np.random.rand(10, 3)
        pre_linkage = linkage(arr, method="average", metric="euclidean")
        dendro = _Dendrogram(data=arr, linkage=pre_linkage, method="average",
                             metric="euclidean", axis=0)
        assert dendro.linkage is pre_linkage
        assert len(dendro.reordered_index) == 10


class TestHierarchicalCluster:
    """Tests for the hierarchical_cluster wrapper function."""

    def test_basic_call(self):
        """Test basic hierarchical_cluster call."""
        np.random.seed(42)
        df = pd.DataFrame(np.random.rand(10, 4))
        result = hierarchical_cluster(data=df, method="average",
                                      metric="euclidean", axis=0)
        assert hasattr(result, "reordered_index")
        assert hasattr(result, "data")
        assert hasattr(result, "linkage")

    def test_default_axis_is_1(self):
        """Test that default axis is 1 (columns)."""
        np.random.seed(42)
        df = pd.DataFrame(np.random.rand(5, 10))
        result = hierarchical_cluster(data=df)
        # axis=1 clusters columns, so reordered_index should have 10 elements
        assert len(result.reordered_index) == 10

    def test_different_metrics(self):
        """Test various distance metrics."""
        np.random.seed(42)
        df = pd.DataFrame(np.random.rand(10, 4))
        for metric in ["euclidean", "cityblock", "cosine"]:
            result = hierarchical_cluster(data=df, metric=metric, axis=0)
            assert len(result.reordered_index) == 10

    def test_requires_scipy(self):
        """Verify the function works with scipy (it's installed in test env)."""
        np.random.seed(42)
        df = pd.DataFrame(np.random.rand(5, 3))
        result = hierarchical_cluster(data=df, axis=0)
        assert result is not None
