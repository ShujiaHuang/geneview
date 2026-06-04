"""
UcscTrack - Stub for UCSC Table Browser queries.

This is a placeholder that mirrors Gviz's UcscTrack API but raises
NotImplementedError on draw, since it requires rtracklayer or a live
UCSC connection.

Use the appropriate ``read_*`` function with a downloaded file instead.
"""

from typing import Optional, Dict, Any, Union

from ._base import Track, GenomicInterval


class UcscTrack(Track):
    """Stub for a track that queries the UCSC Table Browser.

    This class accepts the same constructor parameters as the Gviz
    ``UcscTrack`` but does **not** implement live UCSC queries.
    Calling ``draw()`` raises ``NotImplementedError``.

    To achieve equivalent functionality, download the relevant file from
    the UCSC Table Browser and use the appropriate reader::

        from geneview.genometracks import read_bed, AnnotationTrack

        data = read_bed("ucsc_track.bed")
        track = AnnotationTrack(data)

    Parameters
    ----------
    genome : str, optional
        Genome assembly (e.g. ``"hg38"``).
    chromosome : str, optional
        Chromosome to query.
    track : str, optional
        UCSC track name.
    table : str, optional
        UCSC table name.
    from_ : int, optional
        Start position (named ``from_`` to avoid Python keyword conflict).
    to : int, optional
        End position.
    track_type : str, optional
        Track type hint (e.g. ``"bed"``, ``"bigWig"``).
    name : str
        Track name.
    height : float
        Relative track height.
    display_params : dict, optional
        Display parameters.
    **kwargs
        Additional field mapping keyword arguments.
    """

    def __init__(
        self,
        genome: Optional[str] = None,
        chromosome: Optional[str] = None,
        track: Optional[str] = None,
        table: Optional[str] = None,
        from_: Optional[int] = None,
        to: Optional[int] = None,
        track_type: Optional[str] = None,
        name: str = "UCSC",
        height: float = 1.0,
        display_params: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        super().__init__(name=name, height=height, display_params=display_params)
        self.genome = genome
        self.chromosome = chromosome
        self.track = track
        self.table = table
        self.from_ = from_
        self.to = to
        self.track_type = track_type
        self._field_mapping = kwargs

    def draw(self, ax, region: GenomicInterval) -> None:
        """Raise NotImplementedError — UCSC queries are not implemented.

        Use the appropriate read_* function with a downloaded file instead.
        """
        raise NotImplementedError(
            "UcscTrack requires rtracklayer/UCSC connection. "
            "Use the appropriate read_* function with a downloaded file "
            "instead. Example:\n"
            "    from geneview.genometracks import read_bed, AnnotationTrack\n"
            "    data = read_bed('ucsc_track.bed')\n"
            "    track = AnnotationTrack(data)"
        )

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the requested region if chromosome/from_/to were provided."""
        if self.chromosome and self.from_ is not None and self.to is not None:
            return GenomicInterval(self.chromosome, self.from_, self.to)
        return None
