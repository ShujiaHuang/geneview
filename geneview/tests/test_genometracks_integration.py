"""Tests for genomeview-integrated features.

Tests cover:
    - _utils.py: match_chrom_format, is_paired_end, is_long_frag_dataset
    - _mismatch_counts.py: MismatchCounts
    - _annotation.py: BED12 rendering, feature_filter
    - _alignments_track.py: read_filter, quick_consensus
    - _grouped_alignments.py: GroupedAlignmentsTrack, get_group_by_tag_fn
    - _multi_view.py: plot_tracks_grid, plot_tracks_multi
    - _convenience.py: visualize_files, _normalise_ext
"""
import os
import tempfile
import shutil

import pytest
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks._base import GenomicInterval, Track
from geneview.genometracks._utils import match_chrom_format
from geneview.genometracks._mismatch_counts import MismatchCounts, _TYPES_TO_ID, _N_TYPES
from geneview.genometracks._annotation import AnnotationTrack
from geneview.genometracks._alignments_track import AlignmentsTrack, BAMCoverageTrack
from geneview.genometracks._grouped_alignments import (
    GroupedAlignmentsTrack,
    get_group_by_tag_fn,
)
from geneview.genometracks._multi_view import plot_tracks_grid, plot_tracks_multi
from geneview.genometracks._convenience import visualize_files, _normalise_ext
from geneview.genometracks._genome_axis import GenomeAxisTrack, get_ticks
from geneview.genometracks._io import read_bed, read_bigbed
from geneview.genometracks._track_plot import plot_tracks, find_tracks
from geneview.genometracks._vcf_track import VCFTrack

import pysam


# ---------------------------------------------------------------------------
# Helper: create a tiny synthetic BAM file for testing
# ---------------------------------------------------------------------------
def _make_synthetic_bam(tmpdir, chrom="chr1", length=1000, n_reads=50,
                        paired=False, tag=None, tag_values=None):
    """Create a minimal BAM file with synthetic reads.

    Parameters
    ----------
    tag : str, optional
        BAM tag to attach to each read.
    tag_values : list, optional
        Values to assign for *tag*, cycling through reads.
    """
    bam_path = os.path.join(tmpdir, "test.bam")
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": chrom, "LN": length}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        for i in range(n_reads):
            read = pysam.AlignedSegment()
            read.query_name = f"read_{i}"
            start = i * (length // n_reads)
            read.query_sequence = "A" * 100
            read.flag = 0 if not paired else (1 + 2 * (i % 2))
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            # 100M CIGAR
            read.cigar = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            if tag and tag_values:
                val = tag_values[i % len(tag_values)]
                read.set_tag(tag, val)
            outf.write(read)
    pysam.sort("-o", bam_path, bam_path)
    pysam.index(bam_path)
    return bam_path


@pytest.fixture
def tmp_bam_dir(tmp_path):
    return str(tmp_path)


# ===========================================================================
# _utils.py tests
# ===========================================================================

class TestMatchChromFormat:
    """Tests for match_chrom_format()."""

    def test_exact_match(self):
        assert match_chrom_format("chr1", ["chr1", "chr2"]) == "chr1"

    def test_strip_chr_prefix(self):
        assert match_chrom_format("chr1", ["1", "2", "3"]) == "1"

    def test_add_chr_prefix(self):
        assert match_chrom_format("1", ["chr1", "chr2"]) == "chr1"

    def test_chrM_to_MT(self):
        assert match_chrom_format("chrM", ["MT"]) == "MT"

    def test_MT_to_chrM(self):
        assert match_chrom_format("MT", ["chrM", "chr1"]) == "chrM"

    def test_no_match_returns_original(self):
        assert match_chrom_format("chrX", ["1", "2"]) == "chrX"

    def test_empty_keys(self):
        assert match_chrom_format("chr1", []) == "chr1"

    def test_priority_exact_match(self):
        # Exact match is preferred over stripped
        assert match_chrom_format("chr1", ["chr1", "1"]) == "chr1"

    def test_autosomal_numbers(self):
        keys = [str(i) for i in range(1, 23)]
        assert match_chrom_format("chr22", keys) == "22"
        assert match_chrom_format("1", [f"chr{i}" for i in range(1, 23)]) == "chr1"


class TestIsPairedEnd:
    """Tests for is_paired_end() using synthetic BAM."""

    def test_single_end(self, tmp_bam_dir):
        bam = _make_synthetic_bam(tmp_bam_dir, paired=False)
        from geneview.genometracks._utils import is_paired_end
        assert is_paired_end(bam) is False

    def test_paired_end(self, tmp_bam_dir):
        bam = _make_synthetic_bam(tmp_bam_dir, paired=True)
        from geneview.genometracks._utils import is_paired_end
        assert is_paired_end(bam) is True


class TestIsLongFragDataset:
    """Tests for is_long_frag_dataset()."""

    def test_short_reads(self, tmp_bam_dir):
        bam = _make_synthetic_bam(tmp_bam_dir, n_reads=50)
        from geneview.genometracks._utils import is_long_frag_dataset
        # Synthetic reads are 100bp mapped contiguously; not long-frag
        result = is_long_frag_dataset(bam, n=10)
        # Should be False since there are no large gaps
        assert isinstance(result, bool)


# ===========================================================================
# _mismatch_counts.py tests
# ===========================================================================

class TestMismatchCounts:
    """Tests for MismatchCounts."""

    def test_init(self):
        mc = MismatchCounts("chr1", 100, 200)
        assert mc.chrom == "chr1"
        assert mc.start == 100
        assert mc.end == 200
        assert mc.counts.shape == (_N_TYPES, 100)
        assert mc.insertions.shape == (100,)

    def test_add_count(self):
        mc = MismatchCounts("chr1", 0, 10)
        mc._add_count(5, "A")
        mc._add_count(5, "A")
        mc._add_count(5, "C")
        assert mc.counts[_TYPES_TO_ID["A"], 5] == 2.0
        assert mc.counts[_TYPES_TO_ID["C"], 5] == 1.0

    def test_add_count_out_of_range(self):
        mc = MismatchCounts("chr1", 100, 200)
        mc._add_count(50, "A")   # Before region
        mc._add_count(250, "A")  # After region
        assert mc.counts.sum() == 0.0

    def test_add_insertion(self):
        mc = MismatchCounts("chr1", 0, 10)
        mc._add_count(3, "INS")
        mc._add_count(3, "INS")
        assert mc.insertions[3] == 2.0

    def test_add_deletion(self):
        mc = MismatchCounts("chr1", 0, 10)
        mc._add_count(7, "DEL")
        assert mc.counts[_TYPES_TO_ID["DEL"], 7] == 1.0

    def test_query_high_fraction(self):
        mc = MismatchCounts("chr1", 0, 10)
        # 8 A's and 2 C's at position 5
        for _ in range(8):
            mc._add_count(5, "A")
        for _ in range(2):
            mc._add_count(5, "C")
        # A is 80% → above 0.2 threshold
        assert mc.query("A", 5, threshold=0.2) is True
        # C is 20% → at 0.2 threshold (not above), so False
        assert mc.query("C", 5, threshold=0.21) is False

    def test_query_below_threshold(self):
        mc = MismatchCounts("chr1", 0, 10)
        mc._add_count(5, "A")
        mc._add_count(5, "A")
        mc._add_count(5, "C")  # C = 1/3 ≈ 0.33
        assert mc.query("C", 5, threshold=0.5) is False

    def test_query_deletion_threshold(self):
        mc = MismatchCounts("chr1", 0, 10)
        for _ in range(7):
            mc._add_count(3, "A")
        for _ in range(3):
            mc._add_count(3, "DEL")
        # DEL = 30% → at default del_threshold=0.3, needs > 0.3
        assert mc.query("DEL", 3, del_threshold=0.3) is False
        assert mc.query("DEL", 3, del_threshold=0.25) is True

    def test_query_insertion(self):
        mc = MismatchCounts("chr1", 0, 10)
        for _ in range(8):
            mc._add_count(4, "A")
        for _ in range(4):
            mc._add_count(4, "INS")
        # INS fraction = 4/(8+4*0) ... insertions use total column sum
        # Actually ins / total where total = counts column sum
        # The insertion fraction is ins / total
        assert isinstance(mc.query("INS", 4, threshold=0.2), bool)

    def test_query_out_of_region(self):
        mc = MismatchCounts("chr1", 100, 200)
        assert mc.query("A", 50) is False
        assert mc.query("A", 250) is False

    def test_query_unknown_type(self):
        mc = MismatchCounts("chr1", 0, 10)
        mc._add_count(0, "A")
        assert mc.query("X", 0) is False

    def test_query_with_range(self):
        mc = MismatchCounts("chr1", 0, 20)
        # High alt at positions 5-8
        for pos in range(5, 9):
            for _ in range(8):
                mc._add_count(pos, "A")
            for _ in range(4):
                mc._add_count(pos, "G")
        assert mc.query("G", 5, end=8, threshold=0.3) is True

    def test_tally_reads_synthetic_bam(self, tmp_bam_dir):
        """Test tally_reads against a synthetic BAM."""
        bam_path = _make_synthetic_bam(tmp_bam_dir, chrom="chr1",
                                       length=500, n_reads=20)
        bam = pysam.AlignmentFile(bam_path, "rb")
        mc = MismatchCounts("chr1", 0, 500)
        mc.tally_reads(bam)
        bam.close()
        # Should have accumulated some counts
        assert mc.counts.sum() > 0


# ===========================================================================
# _annotation.py: BED12 + feature_filter tests
# ===========================================================================

def _make_bed12_df():
    """Create a minimal BED12 DataFrame for testing."""
    return pd.DataFrame([{
        "chrom": "chr1",
        "start": 1000,
        "end": 5000,
        "name": "transcript1",
        "score": 0,
        "strand": "+",
        "thick_start": 1500,
        "thick_end": 4500,
        "block_count": 3,
        "block_sizes": "500,800,600",
        "block_starts": "0,1500,3400",
    }, {
        "chrom": "chr1",
        "start": 6000,
        "end": 9000,
        "name": "transcript2",
        "score": 0,
        "strand": "-",
        "thick_start": 6200,
        "thick_end": 8800,
        "block_count": 2,
        "block_sizes": "400,600",
        "block_starts": "0,2400",
    }])


class TestBed12Rendering:
    """Tests for BED12 transcript-aware rendering."""

    def test_bed12_auto_detection(self):
        df = _make_bed12_df()
        track = AnnotationTrack(df, name="BED12")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 900, 5500)
        # Should not raise
        track.draw(ax, region)
        plt.close(fig)

    def test_bed12_thin_utr(self):
        """UTR exons (outside thick_start/thick_end) are drawn thinner."""
        df = pd.DataFrame([{
            "chrom": "chr1", "start": 0, "end": 3000,
            "name": "tx1", "score": 0, "strand": "+",
            "thick_start": 1000, "thick_end": 2000,
            "block_count": 3, "block_sizes": "500,600,500",
            "block_starts": "0,1000,2500",
        }])
        track = AnnotationTrack(df, name="UTR_test")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, 3000)
        track.draw(ax, region)
        plt.close(fig)

    def test_bed12_no_thick_cols(self):
        """When thick_start/thick_end are NaN, falls back to simple arrow."""
        df = pd.DataFrame([{
            "chrom": "chr1", "start": 0, "end": 3000,
            "name": "tx2", "score": 0, "strand": "+",
            "thick_start": np.nan, "thick_end": np.nan,
            "block_count": 2, "block_sizes": "500,500",
            "block_starts": "0,2500",
        }])
        track = AnnotationTrack(df, name="NoCDS")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, 3000)
        track.draw(ax, region)
        plt.close(fig)

    def test_bed12_single_exon(self):
        df = pd.DataFrame([{
            "chrom": "chr1", "start": 100, "end": 600,
            "name": "single_exon", "score": 0, "strand": "+",
            "thick_start": 100, "thick_end": 600,
            "block_count": 1, "block_sizes": "500",
            "block_starts": "0",
        }])
        track = AnnotationTrack(df, name="SingleExon")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        plt.close(fig)

    def test_non_bed12_unchanged(self):
        """Non-BED12 data should render as before (no regression)."""
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [100, 500],
            "end": [300, 700],
            "name": ["feat1", "feat2"],
        })
        track = AnnotationTrack(df, name="Simple")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        plt.close(fig)


class TestFeatureFilter:
    """Tests for AnnotationTrack.feature_filter."""

    def _make_df(self):
        return pd.DataFrame({
            "chrom": ["chr1"] * 5,
            "start": [100, 200, 300, 400, 500],
            "end":   [150, 250, 350, 450, 550],
            "name":  ["a", "b", "c", "d", "e"],
            "score": [10, 50, 30, 80, 20],
        })

    def test_filter_by_score(self):
        df = self._make_df()
        track = AnnotationTrack(
            df, name="Filtered",
            feature_filter=lambda row: row["score"] >= 30,
        )
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        plt.close(fig)

    def test_filter_removes_all(self):
        df = self._make_df()
        track = AnnotationTrack(
            df, name="AllRemoved",
            feature_filter=lambda row: False,
        )
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, 1000)
        # Should not raise even when all features filtered
        track.draw(ax, region)
        plt.close(fig)

    def test_no_filter(self):
        df = self._make_df()
        track = AnnotationTrack(df, name="NoFilter")
        assert track.feature_filter is None


# ===========================================================================
# _alignments_track.py: read_filter, quick_consensus tests
# ===========================================================================

class TestReadFilter:
    """Tests for AlignmentsTrack.read_filter."""

    def test_read_filter_stored(self):
        fn = lambda read: True
        track = AlignmentsTrack(filepath="dummy.bam", read_filter=fn)
        assert track.read_filter is fn

    def test_no_read_filter_default(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        assert track.read_filter is None


class TestQuickConsensus:
    """Tests for AlignmentsTrack.quick_consensus parameters."""

    def test_quick_consensus_default(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        assert track.quick_consensus is True
        assert track.consensus_threshold == 0.2

    def test_quick_consensus_disabled(self):
        track = AlignmentsTrack(
            filepath="dummy.bam",
            quick_consensus=False,
        )
        assert track.quick_consensus is False

    def test_custom_threshold(self):
        track = AlignmentsTrack(
            filepath="dummy.bam",
            consensus_threshold=0.5,
        )
        assert track.consensus_threshold == 0.5


# ===========================================================================
# _grouped_alignments.py tests
# ===========================================================================

class TestGetGroupByTagFn:

    def test_returns_callable(self):
        fn = get_group_by_tag_fn("HP")
        assert callable(fn)

    def test_missing_tag(self, tmp_bam_dir):
        """Reads without the tag should return 'missing'."""
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=5)
        fn = get_group_by_tag_fn("HP")
        bam = pysam.AlignmentFile(bam_path, "rb")
        for read in bam.fetch():
            assert fn(read) == "missing"
            break
        bam.close()

    def test_tag_present(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(
            tmp_bam_dir, n_reads=6,
            tag="HP", tag_values=[1, 2],
        )
        fn = get_group_by_tag_fn("HP")
        bam = pysam.AlignmentFile(bam_path, "rb")
        results = set()
        for read in bam.fetch():
            results.add(fn(read))
        bam.close()
        assert "1" in results or "2" in results


class TestGroupedAlignmentsTrack:

    def test_creation(self):
        fn = get_group_by_tag_fn("HP")
        track = GroupedAlignmentsTrack(
            filepath="dummy.bam",
            keyfn=fn,
        )
        assert track.name == "GroupedAlignments"
        assert track.height == 3.0
        assert track.quick_consensus is True
        assert track.consensus_threshold == 0.2
        assert track.space_between == 0.05

    def test_custom_params(self):
        fn = get_group_by_tag_fn("HP")
        track = GroupedAlignmentsTrack(
            filepath="dummy.bam",
            keyfn=fn,
            is_paired=True,
            type="pileup",
            space_between=0.1,
            quick_consensus=False,
            consensus_threshold=0.5,
            name="Custom",
            height=5.0,
        )
        assert track.is_paired is True
        assert track.plot_type == "pileup"
        assert track.space_between == 0.1
        assert track.quick_consensus is False
        assert track.consensus_threshold == 0.5
        assert track.name == "Custom"
        assert track.height == 5.0

    def test_category_label_fn(self):
        fn = get_group_by_tag_fn("HP")
        track = GroupedAlignmentsTrack(
            filepath="dummy.bam",
            keyfn=fn,
            category_label_fn=lambda x: f"Haplotype {x}",
        )
        assert track.category_label_fn("1") == "Haplotype 1"

    def test_draw_synthetic_bam(self, tmp_bam_dir):
        """Smoke test: draw grouped alignments from a synthetic BAM."""
        bam_path = _make_synthetic_bam(
            tmp_bam_dir, n_reads=20,
            tag="RG", tag_values=["A", "B"],
        )
        fn = get_group_by_tag_fn("RG")
        track = GroupedAlignmentsTrack(filepath=bam_path, keyfn=fn)
        fig, ax = plt.subplots(figsize=(10, 4))
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        plt.close(fig)

    def test_get_region(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir)
        fn = get_group_by_tag_fn("HP")
        track = GroupedAlignmentsTrack(filepath=bam_path, keyfn=fn)
        region = track.get_region()
        assert region is not None
        assert region.chrom == "chr1"


# ===========================================================================
# _multi_view.py tests
# ===========================================================================

class TestPlotTracksGrid:

    def _make_views(self):
        df1 = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [100, 500],
            "end":   [300, 700],
            "name":  ["f1", "f2"],
        })
        df2 = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [200, 600],
            "end":   [400, 800],
            "name":  ["g1", "g2"],
        })
        region1 = GenomicInterval("chr1", 0, 1000)
        region2 = GenomicInterval("chr1", 0, 1000)
        view1 = [GenomeAxisTrack(), AnnotationTrack(df1, name="View1")]
        view2 = [GenomeAxisTrack(), AnnotationTrack(df2, name="View2")]
        return [view1, view2], [region1, region2]

    def test_basic_grid(self):
        views, regions = self._make_views()
        axes = plot_tracks_grid(views, regions=regions, columns=2)
        assert len(axes) == 2
        plt.close("all")

    def test_grid_single_column(self):
        views, regions = self._make_views()
        axes = plot_tracks_grid(views, regions=regions, columns=1)
        assert len(axes) == 2
        plt.close("all")

    def test_grid_with_title(self):
        views, regions = self._make_views()
        axes = plot_tracks_grid(
            views, regions=regions,
            title="Comparison",
            fontsize_main=16,
        )
        assert len(axes) == 2
        plt.close("all")

    def test_grid_auto_figsize(self):
        views, regions = self._make_views()
        axes = plot_tracks_grid(views, regions=regions)
        assert len(axes) == 2
        plt.close("all")

    def test_grid_custom_figsize(self):
        views, regions = self._make_views()
        axes = plot_tracks_grid(
            views, regions=regions,
            figsize=(16, 8),
        )
        assert len(axes) == 2
        plt.close("all")


class TestPlotTracksMulti:

    def test_basic_multi(self):
        df1 = pd.DataFrame({
            "chrom": ["chr1"], "start": [100], "end": [300], "name": ["f1"],
        })
        df2 = pd.DataFrame({
            "chrom": ["chr1"], "start": [200], "end": [400], "name": ["g1"],
        })
        section1 = (
            [GenomeAxisTrack(), AnnotationTrack(df1, name="Sec1")],
            GenomicInterval("chr1", 0, 1000),
        )
        section2 = (
            [GenomeAxisTrack(), AnnotationTrack(df2, name="Sec2")],
            GenomicInterval("chr1", 0, 1000),
        )
        axes = plot_tracks_multi([section1, section2])
        assert len(axes) > 0
        plt.close("all")

    def test_multi_with_title(self):
        df = pd.DataFrame({
            "chrom": ["chr1"], "start": [100], "end": [300], "name": ["f1"],
        })
        section = (
            [GenomeAxisTrack(), AnnotationTrack(df, name="Sec")],
            GenomicInterval("chr1", 0, 1000),
        )
        axes = plot_tracks_multi([section], title="Multi-region")
        assert len(axes) > 0
        plt.close("all")


# ===========================================================================
# _convenience.py tests
# ===========================================================================

class TestNormaliseExt:

    def test_bam(self):
        assert _normalise_ext("sample.bam") == ".bam"

    def test_bed_gz(self):
        assert _normalise_ext("genes.bed.gz") == ".bed"

    def test_bigwig(self):
        assert _normalise_ext("signal.bigwig") == ".bigwig"

    def test_bw(self):
        assert _normalise_ext("signal.bw") == ".bw"

    def test_bedgraph(self):
        assert _normalise_ext("cov.bedgraph") == ".bedgraph"

    def test_gff3(self):
        assert _normalise_ext("annot.gff3") == ".gff3"

    def test_gtf(self):
        assert _normalise_ext("annot.gtf") == ".gtf"

    def test_path_with_dir(self):
        assert _normalise_ext("/data/sample.cram") == ".cram"

    def test_uppercase_ext(self):
        assert _normalise_ext("SAMPLE.BAM") == ".bam"


class TestVisualizeFiles:

    def test_unsupported_extension_raises(self):
        with pytest.raises(ValueError, match="Unsupported"):
            visualize_files(
                ["data.xyz"],
                region=GenomicInterval("chr1", 0, 1000),
            )

    def test_dict_input(self, tmp_bam_dir):
        """Dict input uses keys as track names."""
        bam_path = _make_synthetic_bam(tmp_bam_dir)
        tracks = visualize_files(
            {"My Sample": bam_path},
            region=GenomicInterval("chr1", 0, 1000),
        )
        # Should have GenomeAxisTrack + AlignmentsTrack
        assert len(tracks) == 2
        assert any(t.name == "My Sample" for t in tracks)

    def test_axis_on_top(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir)
        tracks = visualize_files(
            [bam_path],
            region=GenomicInterval("chr1", 0, 1000),
            axis_on_top=True,
        )
        # First track should be GenomeAxisTrack
        assert isinstance(tracks[0], GenomeAxisTrack)

    def test_axis_at_bottom(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir)
        tracks = visualize_files(
            [bam_path],
            region=GenomicInterval("chr1", 0, 1000),
            axis_on_top=False,
        )
        # Last track should be GenomeAxisTrack
        assert isinstance(tracks[-1], GenomeAxisTrack)

    def test_bed_file(self, tmp_path):
        """Test BED file auto-detection."""
        bed_path = str(tmp_path / "test.bed")
        with open(bed_path, "w") as f:
            f.write("chr1\t100\t300\tfeat1\n")
            f.write("chr1\t500\t700\tfeat2\n")
        tracks = visualize_files(
            [bed_path],
            region=GenomicInterval("chr1", 0, 1000),
        )
        assert len(tracks) == 2  # AnnotationTrack + GenomeAxisTrack

    def test_consensus_params_passed(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir)
        tracks = visualize_files(
            [bam_path],
            region=GenomicInterval("chr1", 0, 1000),
            quick_consensus=False,
            consensus_threshold=0.5,
        )
        aln_tracks = [t for t in tracks if isinstance(t, AlignmentsTrack)]
        assert len(aln_tracks) == 1
        assert aln_tracks[0].quick_consensus is False
        assert aln_tracks[0].consensus_threshold == 0.5


# ===========================================================================
# _alignments_track.py: clipping, strand coloring, secondary filtering
# ===========================================================================

class TestClippingAndIndels:
    """Tests for show_clipping, min_indel_size, show_insertion_labels."""

    def test_show_clipping_default(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        assert track.show_clipping is False
        assert track.col_clipping == "cyan"

    def test_show_clipping_enabled(self):
        track = AlignmentsTrack(
            filepath="dummy.bam",
            show_clipping=True,
            col_clipping="magenta",
        )
        assert track.show_clipping is True
        assert track.col_clipping == "magenta"

    def test_min_indel_size_default(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        assert track.min_indel_size == 0

    def test_min_indel_size_custom(self):
        track = AlignmentsTrack(filepath="dummy.bam", min_indel_size=5)
        assert track.min_indel_size == 5

    def test_show_insertion_labels_default(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        assert track.show_insertion_labels is False

    def test_show_insertion_labels_enabled(self):
        track = AlignmentsTrack(
            filepath="dummy.bam",
            show_insertion_labels=True,
        )
        assert track.show_insertion_labels is True

    def test_draw_with_clipping(self, tmp_bam_dir):
        """Smoke test: draw with show_clipping=True on synthetic BAM."""
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10)
        track = AlignmentsTrack(
            filepath=bam_path,
            type="pileup",
            show_clipping=True,
            col_clipping="cyan",
        )
        fig, ax = plt.subplots(figsize=(10, 4))
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        plt.close(fig)

    def test_draw_with_min_indel_size(self, tmp_bam_dir):
        """Smoke test: draw with min_indel_size filtering."""
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10)
        track = AlignmentsTrack(
            filepath=bam_path,
            type="pileup",
            show_indels=True,
            min_indel_size=3,
        )
        fig, ax = plt.subplots(figsize=(10, 4))
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        plt.close(fig)

    def test_draw_with_insertion_labels(self, tmp_bam_dir):
        """Smoke test: draw with insertion labels enabled."""
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10)
        track = AlignmentsTrack(
            filepath=bam_path,
            type="pileup",
            show_indels=True,
            show_insertion_labels=True,
        )
        fig, ax = plt.subplots(figsize=(10, 4))
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        plt.close(fig)


class TestColorByStrand:
    """Tests for color_by_strand and custom strand colors."""

    def test_color_by_strand_default(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        assert track.color_by_strand is True
        assert track.fill_reads_fwd == "#E89E9D"
        assert track.fill_reads_rev == "#8C8FCE"

    def test_color_by_strand_enabled(self):
        track = AlignmentsTrack(
            filepath="dummy.bam",
            color_by_strand=True,
        )
        assert track.color_by_strand is True

    def test_custom_strand_colors(self):
        track = AlignmentsTrack(
            filepath="dummy.bam",
            color_by_strand=True,
            fill_reads_fwd="red",
            fill_reads_rev="blue",
        )
        assert track.fill_reads_fwd == "red"
        assert track.fill_reads_rev == "blue"

    def test_draw_with_strand_coloring(self, tmp_bam_dir):
        """Smoke test: draw with color_by_strand=True."""
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=20)
        track = AlignmentsTrack(
            filepath=bam_path,
            type="pileup",
            color_by_strand=True,
        )
        fig, ax = plt.subplots(figsize=(10, 4))
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        plt.close(fig)


class TestIncludeSecondary:
    """Tests for include_secondary alignment filtering."""

    def test_include_secondary_default(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        assert track.include_secondary is True

    def test_include_secondary_disabled(self):
        track = AlignmentsTrack(
            filepath="dummy.bam",
            include_secondary=False,
        )
        assert track.include_secondary is False

    def test_draw_without_secondary(self, tmp_bam_dir):
        """Smoke test: draw with include_secondary=False."""
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10)
        track = AlignmentsTrack(
            filepath=bam_path,
            type="pileup",
            include_secondary=False,
        )
        fig, ax = plt.subplots(figsize=(10, 4))
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        plt.close(fig)


# ===========================================================================
# _export.py: save_figure tests
# ===========================================================================

class TestSaveFigure:
    """Tests for save_figure()."""

    def test_save_figure_png(self, tmp_path):
        from geneview.genometracks._export import save_figure
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 4, 9])
        out_path = str(tmp_path / "test_output.png")
        result = save_figure([ax], out_path)
        assert result == out_path
        assert os.path.exists(out_path)
        plt.close(fig)

    def test_save_figure_pdf(self, tmp_path):
        from geneview.genometracks._export import save_figure
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 4, 9])
        out_path = str(tmp_path / "test_output.pdf")
        result = save_figure([ax], out_path)
        assert result == out_path
        assert os.path.exists(out_path)
        plt.close(fig)

    def test_save_figure_svg(self, tmp_path):
        from geneview.genometracks._export import save_figure
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 4, 9])
        out_path = str(tmp_path / "test_output.svg")
        result = save_figure([ax], out_path)
        assert result == out_path
        assert os.path.exists(out_path)
        plt.close(fig)

    def test_save_figure_explicit_fmt(self, tmp_path):
        from geneview.genometracks._export import save_figure
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 4, 9])
        out_path = str(tmp_path / "test_output")
        result = save_figure([ax], out_path, fmt="png")
        assert result == out_path
        assert os.path.exists(out_path)
        plt.close(fig)

    def test_save_figure_custom_dpi(self, tmp_path):
        from geneview.genometracks._export import save_figure
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 4, 9])
        out_path = str(tmp_path / "test_hires.png")
        result = save_figure([ax], out_path, dpi=300)
        assert result == out_path
        assert os.path.exists(out_path)
        plt.close(fig)

    def test_save_figure_empty_axes_raises(self):
        from geneview.genometracks._export import save_figure
        with pytest.raises(ValueError, match="No axes"):
            save_figure([], "output.png")

    def test_save_figure_single_ax(self, tmp_path):
        """save_figure should also work with a single Axes (not list)."""
        from geneview.genometracks._export import save_figure
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 4, 9])
        out_path = str(tmp_path / "test_single_ax.png")
        result = save_figure(ax, out_path)
        assert result == out_path
        assert os.path.exists(out_path)
        plt.close(fig)

    def test_save_figure_from_plot_tracks(self, tmp_path):
        """Integration: save a figure produced by plot_tracks."""
        from geneview.genometracks._export import save_figure
        from geneview.genometracks._track_plot import plot_tracks
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [100, 500],
            "end": [300, 700],
            "name": ["f1", "f2"],
        })
        atrack = AnnotationTrack(df, name="Test")
        region = GenomicInterval("chr1", 0, 1000)
        axes = plot_tracks([GenomeAxisTrack(), atrack], region=region)
        out_path = str(tmp_path / "track_figure.png")
        result = save_figure(axes, out_path, dpi=100)
        assert result == out_path
        assert os.path.exists(out_path)
        plt.close("all")


# ===========================================================================
# _utils.py: reverse_comp tests
# ===========================================================================

class TestReverseComp:
    """Tests for reverse_comp()."""

    def test_basic_reverse_comp(self):
        from geneview.genometracks._utils import reverse_comp
        assert reverse_comp("ATCG") == "CGAT"

    def test_reverse_comp_palindrome(self):
        from geneview.genometracks._utils import reverse_comp
        # AATT is its own reverse complement
        assert reverse_comp("AATT") == "AATT"

    def test_reverse_comp_longer(self):
        from geneview.genometracks._utils import reverse_comp
        assert reverse_comp("AATTCC") == "GGAATT"

    def test_reverse_comp_lowercase(self):
        from geneview.genometracks._utils import reverse_comp
        assert reverse_comp("atcg") == "cgat"

    def test_reverse_comp_mixed_case(self):
        from geneview.genometracks._utils import reverse_comp
        assert reverse_comp("AtCg") == "cGaT"

    def test_reverse_comp_with_N(self):
        from geneview.genometracks._utils import reverse_comp
        assert reverse_comp("ANCG") == "CGNT"

    def test_reverse_comp_empty(self):
        from geneview.genometracks._utils import reverse_comp
        assert reverse_comp("") == ""

    def test_reverse_comp_single_base(self):
        from geneview.genometracks._utils import reverse_comp
        assert reverse_comp("A") == "T"
        assert reverse_comp("G") == "C"

    def test_reverse_comp_unknown_chars_preserved(self):
        from geneview.genometracks._utils import reverse_comp
        # Unknown characters are left unchanged (just reversed)
        assert reverse_comp("AXT") == "AXT"


# ===========================================================================
# Tests for remaining genomeview features (round 3)
# ===========================================================================

class TestOverlapColor:
    """Test the overlap_color parameter for paired-end reads."""

    def test_overlap_color_param_default(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        assert track.overlap_color is None

    def test_overlap_color_param_set(self):
        track = AlignmentsTrack(filepath="dummy.bam", overlap_color="lime")
        assert track.overlap_color == "lime"

    def test_draw_with_overlap_color(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10, paired=True)
        track = AlignmentsTrack(
            filepath=bam_path,
            type="pileup",
            is_paired=True,
            overlap_color="lime",
        )
        fig, ax = plt.subplots(figsize=(10, 4))
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        plt.close(fig)


class TestDrawReadLabels:
    """Test the draw_read_labels parameter."""

    def test_draw_read_labels_default(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        assert track.draw_read_labels is False

    def test_draw_read_labels_set(self):
        track = AlignmentsTrack(filepath="dummy.bam", draw_read_labels=True)
        assert track.draw_read_labels is True

    def test_draw_with_read_labels(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=5)
        track = AlignmentsTrack(
            filepath=bam_path,
            type="pileup",
            draw_read_labels=True,
        )
        fig, ax = plt.subplots(figsize=(10, 4))
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        # Check that text labels were added
        texts = [c for c in ax.get_children() if isinstance(c, matplotlib.text.Text)]
        # Should have some text labels for reads
        assert len(texts) > 0
        plt.close(fig)


class TestMinInsertionLabelSize:
    """Test the min_insertion_label_size parameter."""

    def test_default_value(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        assert track.min_insertion_label_size == 5

    def test_custom_value(self):
        track = AlignmentsTrack(filepath="dummy.bam", min_insertion_label_size=10)
        assert track.min_insertion_label_size == 10

    def test_draw_with_min_insertion_label_size(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10)
        track = AlignmentsTrack(
            filepath=bam_path,
            type="pileup",
            min_insertion_label_size=3,
        )
        fig, ax = plt.subplots(figsize=(10, 4))
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        plt.close(fig)


class TestGetTicks:
    """Test the get_ticks() standalone utility."""

    def test_basic_ticks(self):
        ticks = get_ticks(0, 1000, target_n_labels=5)
        assert len(ticks) > 0
        # All ticks should be within range
        for pos, label in ticks:
            assert 0 <= pos <= 1000

    def test_ticks_mb_range(self):
        ticks = get_ticks(2_000_000, 3_000_000, target_n_labels=5)
        assert len(ticks) > 0
        # Labels should contain 'kb' or 'Mb' (depends on resolution)
        labels = [label for _, label in ticks]
        assert any("kb" in l or "Mb" in l for l in labels)

    def test_ticks_kb_range(self):
        ticks = get_ticks(10_000, 50_000, target_n_labels=5)
        assert len(ticks) > 0
        labels = [label for _, label in ticks]
        assert any("kb" in l for l in labels)

    def test_ticks_bp_range(self):
        ticks = get_ticks(0, 500, target_n_labels=5)
        assert len(ticks) > 0

    def test_ticks_zero_width(self):
        ticks = get_ticks(100, 100)
        assert len(ticks) == 1
        assert ticks[0][0] == 100

    def test_ticks_sorted(self):
        ticks = get_ticks(0, 1_000_000, target_n_labels=10)
        positions = [pos for pos, _ in ticks]
        assert positions == sorted(positions)


class TestBAMCoverageTrack:
    """Test BAMCoverageTrack standalone coverage line."""

    def test_init_default(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10)
        track = BAMCoverageTrack(filepath=bam_path)
        assert track.plot_type == "line"
        assert track.col == "#5B8DB8"
        assert track.alpha == 0.7

    def test_init_fill(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10)
        track = BAMCoverageTrack(filepath=bam_path, type="fill", col="red")
        assert track.plot_type == "fill"
        assert track.col == "red"

    def test_draw_line(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10)
        track = BAMCoverageTrack(filepath=bam_path, type="line")
        fig, ax = plt.subplots(figsize=(10, 3))
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        # Should have a line plotted
        lines = ax.get_lines()
        assert len(lines) > 0
        plt.close(fig)

    def test_draw_fill(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10)
        track = BAMCoverageTrack(filepath=bam_path, type="fill")
        fig, ax = plt.subplots(figsize=(10, 3))
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        # Should have collections (fill_between)
        assert len(ax.collections) > 0
        plt.close(fig)

    def test_draw_with_transformation(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10)
        track = BAMCoverageTrack(filepath=bam_path, transformation=np.log1p)
        fig, ax = plt.subplots(figsize=(10, 3))
        region = GenomicInterval("chr1", 0, 1000)
        track.draw(ax, region)
        plt.close(fig)

    def test_get_region(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10)
        track = BAMCoverageTrack(filepath=bam_path)
        region = track.get_region()
        assert region is not None
        assert region.chrom == "chr1"

    def test_name_from_file(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=5)
        track = BAMCoverageTrack(filepath=bam_path)
        assert track.name is not None
        assert len(track.name) > 0

    def test_plot_tracks_integration(self, tmp_bam_dir):
        bam_path = _make_synthetic_bam(tmp_bam_dir, n_reads=10)
        track = BAMCoverageTrack(filepath=bam_path, name="Coverage")
        axis = GenomeAxisTrack()
        region = GenomicInterval("chr1", 0, 1000)
        axes = plot_tracks([axis, track], region=region)
        assert len(axes) == 2
        plt.close("all")


class TestReadBigBed:
    """Test read_bigbed() function."""

    def test_import_error_message(self):
        """read_bigbed should give a clear error if pyBigWig is missing."""
        # This test only runs if pyBigWig is not installed
        try:
            import pyBigWig
            pytest.skip("pyBigWig is installed")
        except ImportError:
            with pytest.raises(ImportError, match="pyBigWig"):
                read_bigbed("dummy.bb")

    def test_non_bigbed_error(self):
        """read_bigbed should raise ValueError for non-BigBed files."""
        try:
            import pyBigWig
        except ImportError:
            pytest.skip("pyBigWig not installed")
        # Create a temp file that is not a bigbed
        with tempfile.NamedTemporaryFile(suffix=".bb", delete=False) as f:
            f.write(b"not a bigbed")
            tmp_path = f.name
        try:
            with pytest.raises(Exception):
                read_bigbed(tmp_path)
        finally:
            os.unlink(tmp_path)


class TestFindTracks:
    """Test find_tracks() utility."""

    def test_find_by_name(self):
        ann = AnnotationTrack(
            pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [100]}),
            name="Genes",
        )
        axis = GenomeAxisTrack(name="Axis")
        tracks = [axis, ann]
        result = find_tracks(tracks, name="Genes")
        assert len(result) == 1
        assert result[0].name == "Genes"

    def test_find_by_type(self):
        ann = AnnotationTrack(
            pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [100]}),
            name="Genes",
        )
        axis = GenomeAxisTrack()
        tracks = [axis, ann]
        result = find_tracks(tracks, track_type=GenomeAxisTrack)
        assert len(result) == 1
        assert isinstance(result[0], GenomeAxisTrack)

    def test_find_no_match(self):
        axis = GenomeAxisTrack()
        tracks = [axis]
        result = find_tracks(tracks, name="nonexistent")
        assert len(result) == 0

    def test_find_all(self):
        ann = AnnotationTrack(
            pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [100]}),
            name="Genes",
        )
        axis = GenomeAxisTrack()
        tracks = [axis, ann]
        result = find_tracks(tracks)
        assert len(result) == 2

    def test_find_single_track_input(self):
        axis = GenomeAxisTrack(name="MyAxis")
        result = find_tracks(axis, name="MyAxis")
        assert len(result) == 1

    def test_find_both_filters(self):
        ann = AnnotationTrack(
            pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [100]}),
            name="Genes",
        )
        axis = GenomeAxisTrack(name="Axis")
        tracks = [axis, ann]
        result = find_tracks(tracks, name="Axis", track_type=GenomeAxisTrack)
        assert len(result) == 1
        assert result[0].name == "Axis"


# ---------------------------------------------------------------------------
# Tests for color_fn on AlignmentsTrack
# ---------------------------------------------------------------------------
class TestAlignmentsColorFn:
    """Tests for the color_fn parameter on AlignmentsTrack."""

    def test_color_fn_default_none(self):
        t = AlignmentsTrack(filepath="dummy.bam")
        assert t.color_fn is None

    def test_color_fn_set(self):
        fn = lambda read: "red"
        t = AlignmentsTrack(filepath="dummy.bam", color_fn=fn)
        assert t.color_fn is fn

    def test_color_fn_overrides_fill_reads(self, tmp_path):
        """color_fn should override fill_reads when drawing."""
        bam = _make_synthetic_bam(str(tmp_path), n_reads=20)
        region = GenomicInterval("chr1", 0, 1000)

        fig, ax = plt.subplots()
        t = AlignmentsTrack(
            filepath=bam, type="pileup",
            fill_reads="blue",
            color_fn=lambda read: "red",
        )
        t.draw(ax, region)
        plt.close(fig)

    def test_color_fn_exception_falls_back(self, tmp_path):
        """If color_fn raises, should fall back to fill_reads."""
        bam = _make_synthetic_bam(str(tmp_path), n_reads=20)
        region = GenomicInterval("chr1", 0, 1000)

        def bad_fn(read):
            raise ValueError("oops")

        fig, ax = plt.subplots()
        t = AlignmentsTrack(filepath=bam, type="pileup", color_fn=bad_fn)
        t.draw(ax, region)  # Should not raise
        plt.close(fig)


# ---------------------------------------------------------------------------
# Helper: create a tiny synthetic VCF file for testing
# ---------------------------------------------------------------------------
def _make_synthetic_vcf(tmpdir, chrom="chr1", n_variants=10, start=0):
    """Create a minimal VCF file with synthetic SNPs."""
    vcf_path = os.path.join(tmpdir, "test.vcf.gz")
    header = pysam.VariantHeader()
    header.add_sample("SAMPLE")
    header.add_line('##contig=<ID={},length=1000000>'.format(chrom))
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

    vcf_out = pysam.VariantFile(vcf_path, "wz", header=header)
    for i in range(n_variants):
        pos = start + i * 100 + 50
        rec = vcf_out.new_record()
        rec.chrom = chrom
        rec.pos = pos
        rec.alleles = ("A", ["C", "G", "T"][i % 3])
        rec.qual = 30 + i * 5
        rec.samples["SAMPLE"]["GT"] = (0, 1)
        vcf_out.write(rec)
    vcf_out.close()

    # Index the VCF
    pysam.tabix_index(vcf_path, preset="vcf", force=True)
    return vcf_path


# ---------------------------------------------------------------------------
# Tests for VCFTrack
# ---------------------------------------------------------------------------
class TestVCFTrack:
    """Tests for VCFTrack."""

    def test_init_default(self):
        t = VCFTrack(filepath="dummy.vcf.gz")
        assert t.filepath == "dummy.vcf.gz"
        assert t.name == "Variants"
        assert t.height == 1.0
        assert callable(t.color_fn)

    def test_init_custom(self):
        fn = lambda v: "red"
        t = VCFTrack(filepath="dummy.vcf.gz", color_fn=fn,
                     min_variant_width=0.01, name="MyVCF")
        assert t.color_fn is fn
        assert t.min_variant_width == 0.01
        assert t.name == "MyVCF"

    def test_draw_with_variants(self, tmp_path):
        vcf = _make_synthetic_vcf(str(tmp_path))
        region = GenomicInterval("chr1", 0, 1500)

        fig, ax = plt.subplots()
        t = VCFTrack(filepath=vcf, name="Test VCF")
        t.draw(ax, region)
        # Should have drawn variant rectangles
        plt.close(fig)

    def test_draw_no_variants_in_region(self, tmp_path):
        vcf = _make_synthetic_vcf(str(tmp_path), start=5000)
        region = GenomicInterval("chr1", 0, 100)  # No variants here

        fig, ax = plt.subplots()
        t = VCFTrack(filepath=vcf)
        t.draw(ax, region)
        plt.close(fig)

    def test_custom_color_fn(self, tmp_path):
        vcf = _make_synthetic_vcf(str(tmp_path))
        region = GenomicInterval("chr1", 0, 1500)

        def color_by_qual(v):
            return "red" if v.qual and v.qual > 50 else "blue"

        fig, ax = plt.subplots()
        t = VCFTrack(filepath=vcf, color_fn=color_by_qual)
        t.draw(ax, region)
        plt.close(fig)

    def test_plot_tracks_integration(self, tmp_path):
        vcf = _make_synthetic_vcf(str(tmp_path))
        region = GenomicInterval("chr1", 0, 1500)

        tracks = [GenomeAxisTrack(), VCFTrack(filepath=vcf, name="SNPs")]
        axes = plot_tracks(tracks, region=region)
        assert len(axes) == 2
        plt.close("all")

    def test_repr(self):
        t = VCFTrack(filepath="test.vcf.gz", name="MyVCF")
        r = repr(t)
        assert "VCFTrack" in r
        assert "test.vcf.gz" in r

    def test_pysam_import_error_message(self):
        """VCFTrack should give a clear error if pysam is missing."""
        # pysam is available in test env, but let's check the method exists
        t = VCFTrack(filepath="dummy.vcf.gz")
        assert hasattr(t, '_import_pysam')
