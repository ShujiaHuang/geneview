"""Tests for the geneview.mtdna mitochondrial visualization module.

Covers the rCRS reference backbone, the MitoQuest-style readers
(VCF / copynum TSV / BAM coverage) and the five plotting functions.

Author: Shujia Huang
"""
import textwrap

import numpy as np
import pandas as pd
import pytest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview import mtdna
from geneview.mtdna._reference import (
    MT_LENGTH,
    get_mt_genes,
    gene_at,
    genes_in_range,
    is_mt_contig,
)
from geneview.mtdna._utils import resolve_feature_colors, coerce_variant_frame
from geneview.plotstyle import use_style


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------
def _write_vcf(tmp_path):
    """Write a tiny two-sample MitoQuest-style VCF and return its path."""
    text = textwrap.dedent("""\
        ##fileformat=VCFv4.2
        ##contig=<ID=chrM,length=16569>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
        ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Heteroplasmy fraction">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
        chrM\t263\t.\tA\tG\t.\tPASS\t.\tGT:DP:AF\t1:1200:0.98\t1:900:0.95
        chrM\t3243\t.\tA\tG\t.\tPASS\t.\tGT:DP:AF\t0/1:1500:0.12\t0:1400:0.0
        chrM\t8860\t.\tA\tG,T\t.\tPASS\t.\tGT:DP:AF\t1/2:1000:0.6,0.3\t1:800:0.99
    """)
    path = tmp_path / "cohort.mt.vcf"
    path.write_text(text)
    return str(path)


def _write_copynum(tmp_path, name="S1.cn.tsv", cn=500.0):
    """Write a tiny MitoQuest copynum TSV (chr1 + chrM rows)."""
    text = textwrap.dedent("""\
        #Chromosome\tFragments\tChrom_Length\tGC_Content\tFragment_Normalized_Ratio\tCopyNum\tCopyNum-CI95-Lower\tCopyNum-CI95-Upper\tEffective_Length\tRegions_Used
        chr1\t100\t249000000\t0.41\t1.0\t2.0\t1.8\t2.2\t240000000\tauto
        chrM\t5000\t16569\t0.44\t250.0\t{cn}\t{lo}\t{hi}\t16000\tall
    """).format(cn=cn, lo=cn * 0.95, hi=cn * 1.05)
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def _make_bam(tmp_path, name="s1.bam"):
    """Build a tiny indexed BAM over a chrM contig; return its path."""
    pysam = pytest.importorskip("pysam")
    header = {"HD": {"VN": "1.6"},
              "SQ": [{"SN": "chrM", "LN": MT_LENGTH}]}
    path = tmp_path / name
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i in range(50):
            seg = pysam.AlignedSegment()
            seg.query_name = "r%d" % i
            seg.query_sequence = "A" * 100
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = i * 300  # 0 .. 14700, sorted, no wrap
            seg.mapping_quality = 60
            seg.cigartuples = [(0, 100)]
            seg.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(seg)
    pysam.index(str(path))
    return str(path)


def _synthetic_variants():
    """A small tidy variant frame spanning several feature types."""
    return pd.DataFrame({
        "sample": ["S1", "S1", "S2", "S2"],
        "chrom": ["chrM"] * 4,
        "pos": [263, 3243, 8860, 14747],
        "ref": ["A", "A", "A", "C"],
        "alt": ["G", "G", "T", "T"],
        "vaf": [0.98, 0.12, 0.30, 0.75],
        "depth": [1200, 1500, 1000, 900],
        "gt": ["1", "0/1", "1", "0/1"],
        "status": ["HOM", "HET", "HOM", "HET"],
        "var_type": ["SNV", "SNV", "SNV", "SNV"],
        "variant_id": ["."] * 4,
    })


# ---------------------------------------------------------------------------
# Reference backbone
# ---------------------------------------------------------------------------
class TestReference:

    def test_mt_length(self):
        assert MT_LENGTH == 16569

    def test_gene_count_by_type(self):
        genes = get_mt_genes()
        by = {}
        for g in genes:
            by[g["feature_type"]] = by.get(g["feature_type"], 0) + 1
        assert by["protein_coding"] == 13
        assert by["tRNA"] == 22
        assert by["rRNA"] == 2
        # D-loop is stored as two arcs to handle the origin wrap-around.
        assert by["control_region"] == 2

    def test_get_mt_genes_dataframe(self):
        df = get_mt_genes(as_dataframe=True)
        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == ["name", "start", "end", "strand", "feature_type"]
        assert len(df) == 39

    def test_gene_at_hits(self):
        assert gene_at(3308)["name"] == "MT-ND1"
        assert gene_at(263)["feature_type"] == "control_region"
        assert gene_at(8860)["feature_type"] == "protein_coding"

    def test_gene_at_miss_returns_none(self):
        # Position 0 is outside the 1-based genome.
        assert gene_at(0) is None

    def test_genes_in_range(self):
        hits = genes_in_range(3200, 3400)
        names = {g["name"] for g in hits}
        assert "MT-TL1" in names and "MT-ND1" in names

    def test_is_mt_contig(self):
        for name in ("chrM", "MT", "chrMT", "NC_012920.1", "M"):
            assert is_mt_contig(name)
        assert not is_mt_contig("chr1")
        assert not is_mt_contig(None)


# ---------------------------------------------------------------------------
# read_mito_vcf
# ---------------------------------------------------------------------------
class TestReadMitoVcf:

    def test_basic_shape_and_columns(self, tmp_path):
        pytest.importorskip("pysam")
        df = mtdna.read_mito_vcf(_write_vcf(tmp_path))
        assert list(df.columns) == mtdna.MITO_VCF_COLUMNS
        # S1@263, S2@263, S1@3243, S1@8860(x2 alleles), S2@8860 = 6 rows.
        assert len(df) == 6

    def test_het_hom_status(self, tmp_path):
        pytest.importorskip("pysam")
        df = mtdna.read_mito_vcf(_write_vcf(tmp_path))
        het = df[df["status"] == "HET"]
        # The 0/1 call and the two 1/2 alleles are heteroplasmic.
        assert set(het["pos"]) == {3243, 8860}

    def test_multiallelic_expands_to_two_rows(self, tmp_path):
        pytest.importorskip("pysam")
        df = mtdna.read_mito_vcf(_write_vcf(tmp_path))
        s1_8860 = df[(df["sample"] == "S1") & (df["pos"] == 8860)]
        assert len(s1_8860) == 2
        assert set(s1_8860["alt"]) == {"G", "T"}
        # Per-ALT VAF unpacked correctly.
        assert set(np.round(s1_8860["vaf"], 2)) == {0.6, 0.3}

    def test_ref_genotype_dropped(self, tmp_path):
        pytest.importorskip("pysam")
        df = mtdna.read_mito_vcf(_write_vcf(tmp_path))
        # S2 is homozygous-ref at 3243 -> no row.
        assert len(df[(df["sample"] == "S2") & (df["pos"] == 3243)]) == 0

    def test_sample_filter(self, tmp_path):
        pytest.importorskip("pysam")
        df = mtdna.read_mito_vcf(_write_vcf(tmp_path), samples=["S1"])
        assert set(df["sample"]) == {"S1"}

    def test_min_vaf_filter(self, tmp_path):
        pytest.importorskip("pysam")
        df = mtdna.read_mito_vcf(_write_vcf(tmp_path), min_vaf=0.5)
        assert (df["vaf"] >= 0.5).all()
        # The 0.12 (S1@3243) and 0.3 (S1@8860 T) alleles are removed; the
        # 0.98/0.95/0.60/0.99 alleles survive.
        assert len(df) == 4

    def test_missing_file_raises(self, tmp_path):
        pytest.importorskip("pysam")
        with pytest.raises(FileNotFoundError):
            mtdna.read_mito_vcf(str(tmp_path / "nope.vcf"))


# ---------------------------------------------------------------------------
# read_mito_copynumber
# ---------------------------------------------------------------------------
class TestReadMitoCopynumber:

    def test_single_file(self, tmp_path):
        df = mtdna.read_mito_copynumber(_write_copynum(tmp_path))
        assert list(df.columns) == ["sample", "chrom", "copy_number", "ci_low", "ci_high"]
        assert len(df) == 1  # only the chrM row is kept
        assert df.iloc[0]["copy_number"] == 500.0
        assert df.iloc[0]["ci_low"] == pytest.approx(475.0)

    def test_multi_file_sorted_desc(self, tmp_path):
        f1 = _write_copynum(tmp_path, name="S1.cn.tsv", cn=100.0)
        f2 = _write_copynum(tmp_path, name="S2.cn.tsv", cn=800.0)
        df = mtdna.read_mito_copynumber([f1, f2])
        assert len(df) == 2
        # Descending by copy number.
        assert list(df["copy_number"]) == [800.0, 100.0]

    def test_sample_names_from_filename(self, tmp_path):
        df = mtdna.read_mito_copynumber(_write_copynum(tmp_path, name="SampleA.cn.tsv"))
        assert df.iloc[0]["sample"] == "SampleA"

    def test_explicit_sample_names(self, tmp_path):
        df = mtdna.read_mito_copynumber(_write_copynum(tmp_path), sample_names=["X"])
        assert df.iloc[0]["sample"] == "X"

    def test_sample_names_length_mismatch(self, tmp_path):
        with pytest.raises(ValueError):
            mtdna.read_mito_copynumber(_write_copynum(tmp_path), sample_names=["a", "b"])


# ---------------------------------------------------------------------------
# read_mito_coverage
# ---------------------------------------------------------------------------
class TestReadMitoCoverage:

    def test_binned_coverage(self, tmp_path):
        pytest.importorskip("pysam")
        bam = _make_bam(tmp_path)
        df = mtdna.read_mito_coverage([bam], bins=100)
        assert list(df.columns) == ["sample", "chrom", "start", "end", "pos", "depth"]
        assert len(df) >= 100
        assert (df["chrom"] == "chrM").all()
        assert df["depth"].sum() > 0

    def test_sample_label_from_filename(self, tmp_path):
        pytest.importorskip("pysam")
        bam = _make_bam(tmp_path)
        df = mtdna.read_mito_coverage([bam], bins=50)
        assert set(df["sample"]) == {"s1"}

    def test_missing_contig_raises(self, tmp_path):
        pysam = pytest.importorskip("pysam")
        header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
        path = tmp_path / "nuc.bam"
        with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
            pass
        pysam.index(str(path))
        with pytest.raises(ValueError):
            mtdna.read_mito_coverage([str(path)], bins=10)


# ---------------------------------------------------------------------------
# _utils
# ---------------------------------------------------------------------------
class TestUtils:

    def test_resolve_feature_colors_defaults(self):
        colors = resolve_feature_colors()
        assert set(colors) >= {"protein_coding", "tRNA", "rRNA", "control_region"}

    def test_resolve_feature_colors_override(self):
        colors = resolve_feature_colors({"tRNA": "#000000"})
        assert colors["tRNA"] == "#000000"

    def test_coerce_variant_frame_adds_feature_type(self):
        df = coerce_variant_frame(pd.DataFrame({"pos": [3308], "vaf": [0.5]}))
        assert df.iloc[0]["feature_type"] == "protein_coding"

    def test_coerce_variant_frame_requires_pos(self):
        with pytest.raises(ValueError):
            coerce_variant_frame(pd.DataFrame({"vaf": [0.5]}))


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
class TestPlots:

    def test_mito_genome_map_returns_axes(self):
        ax = mtdna.mito_genome_map(_synthetic_variants())
        assert ax.name == "polar"
        plt.close(ax.figure)

    def test_mito_genome_map_bare_ring(self):
        ax = mtdna.mito_genome_map(None)
        assert ax.name == "polar"
        plt.close(ax.figure)

    def test_mito_genome_map_color_by_status(self):
        ax = mtdna.mito_genome_map(_synthetic_variants(), color_by="status")
        assert ax is not None
        plt.close(ax.figure)

    def test_heteroplasmy_scatter(self):
        ax = mtdna.heteroplasmy_scatter(_synthetic_variants())
        assert ax.get_xlim()[1] == MT_LENGTH
        plt.close(ax.figure)

    def test_heteroplasmy_scatter_hue_sample(self):
        ax = mtdna.heteroplasmy_scatter(_synthetic_variants(), hue="sample")
        assert ax is not None
        plt.close(ax.figure)

    def test_heteroplasmy_heatmap(self):
        ax = mtdna.heteroplasmy_heatmap(_synthetic_variants())
        assert ax is not None
        plt.close(ax.figure)

    def test_mito_coverage_plot(self):
        cov = pd.DataFrame({
            "sample": ["s1"] * 100,
            "pos": np.linspace(1, MT_LENGTH, 100),
            "depth": np.random.RandomState(0).poisson(500, 100),
        })
        ax = mtdna.mito_coverage_plot(cov)
        assert ax.get_xlim()[1] == MT_LENGTH
        plt.close(ax.figure)

    def test_mito_copynumber_plot(self):
        cn = pd.DataFrame({
            "sample": ["A", "B", "C"],
            "copy_number": [500.0, 300.0, 700.0],
            "ci_low": [475.0, 285.0, 665.0],
            "ci_high": [525.0, 315.0, 735.0],
        })
        ax = mtdna.mito_copynumber_plot(cn)
        assert ax is not None
        plt.close(ax.figure)

    def test_plots_with_journal_style(self):
        with use_style("nature"):
            ax = mtdna.heteroplasmy_scatter(_synthetic_variants())
        assert ax is not None
        plt.close(ax.figure)
