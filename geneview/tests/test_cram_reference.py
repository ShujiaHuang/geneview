"""Tests for reference-aware BAM/CRAM opening.

`geneview tracks` must decode CRAM files using the user-supplied ``--reference``
FASTA (passed to pysam as ``reference_filename``) instead of relying on the
reference path embedded in the CRAM header.  These tests lock in that the
central ``open_alignment_file`` helper picks the right mode and forwards the
reference, and that every track / helper that reads alignments threads its
``reference`` through to that helper.
"""
import pytest
import matplotlib
matplotlib.use("Agg")

from geneview.genometracks import _io
from geneview.genometracks import _utils
from geneview.genometracks._alignments_track import AlignmentsTrack, BAMCoverageTrack
from geneview.genometracks._grouped_alignments import GroupedAlignmentsTrack


class _FakeAln:
    """Minimal stand-in for a pysam.AlignmentFile handle."""

    def __init__(self, paired=False):
        self.references = ["chr1"]
        self.lengths = [10000]
        self._paired = paired

    def fetch(self, *args, **kwargs):
        read = type("R", (), {"is_paired": self._paired, "query_length": 50})()
        return [read]

    def close(self):
        pass


def _spy(record, fake=None):
    """Return an ``open_alignment_file`` replacement that records its call."""

    def _open(filepath, reference=None, mode=None):
        record.append({"filepath": filepath, "reference": reference, "mode": mode})
        return fake if fake is not None else _FakeAln()

    return _open


class TestOpenAlignmentFile:
    """The helper must pick mode by extension and forward the reference."""

    def _patch_pysam(self, monkeypatch, captured):
        pysam = pytest.importorskip("pysam")

        def fake_af(filepath, mode, **kwargs):
            captured["filepath"] = filepath
            captured["mode"] = mode
            captured["kwargs"] = kwargs
            return object()

        monkeypatch.setattr(pysam, "AlignmentFile", fake_af)

    def test_cram_uses_rc_mode_and_reference_filename(self, monkeypatch):
        captured = {}
        self._patch_pysam(monkeypatch, captured)
        _io.open_alignment_file("sample.cram", reference="ref.fa")
        assert captured["mode"] == "rc"
        assert captured["kwargs"]["reference_filename"] == "ref.fa"

    def test_bam_uses_rb_and_omits_reference(self, monkeypatch):
        captured = {}
        self._patch_pysam(monkeypatch, captured)
        _io.open_alignment_file("sample.bam")
        assert captured["mode"] == "rb"
        assert "reference_filename" not in captured["kwargs"]

    def test_uppercase_cram_extension_detected(self, monkeypatch):
        captured = {}
        self._patch_pysam(monkeypatch, captured)
        _io.open_alignment_file("SAMPLE.CRAM")
        assert captured["mode"] == "rc"

    def test_explicit_mode_overrides_extension(self, monkeypatch):
        captured = {}
        self._patch_pysam(monkeypatch, captured)
        _io.open_alignment_file("sample.cram", mode="rb")
        assert captured["mode"] == "rb"


class TestReferencePlumbing:
    """Each alignment reader must forward its reference to the open helper."""

    def test_alignments_track_forwards_reference(self, monkeypatch):
        import geneview.genometracks._alignments_track as mod
        rec = []
        monkeypatch.setattr(mod, "open_alignment_file", _spy(rec))
        AlignmentsTrack(filepath="s.cram", reference="ref.fa").get_region()
        assert rec and rec[0]["reference"] == "ref.fa"
        assert rec[0]["filepath"] == "s.cram"

    def test_bamcoverage_track_forwards_reference(self, monkeypatch):
        import geneview.genometracks._alignments_track as mod
        rec = []
        monkeypatch.setattr(mod, "open_alignment_file", _spy(rec))
        BAMCoverageTrack(filepath="s.cram", reference="ref.fa").get_region()
        assert rec and rec[0]["reference"] == "ref.fa"

    def test_grouped_alignments_forwards_reference(self, monkeypatch):
        import geneview.genometracks._grouped_alignments as mod
        rec = []
        monkeypatch.setattr(mod, "open_alignment_file", _spy(rec))
        GroupedAlignmentsTrack(filepath="s.cram", keyfn=lambda r: "g",
                               reference="ref.fa").get_region()
        assert rec and rec[0]["reference"] == "ref.fa"

    def test_is_paired_end_forwards_reference(self, monkeypatch):
        rec = []
        monkeypatch.setattr(_io, "open_alignment_file",
                            _spy(rec, fake=_FakeAln(paired=True)))
        assert _utils.is_paired_end("s.cram", reference="ref.fa") is True
        assert rec and rec[0]["reference"] == "ref.fa"

    def test_is_long_frag_dataset_forwards_reference(self, monkeypatch):
        rec = []
        monkeypatch.setattr(_io, "open_alignment_file",
                            _spy(rec, fake=_FakeAln(paired=False)))
        _utils.is_long_frag_dataset("s.cram", reference="ref.fa")
        assert rec and rec[0]["reference"] == "ref.fa"
