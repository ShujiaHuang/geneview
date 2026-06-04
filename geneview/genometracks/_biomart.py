"""
BiomartGeneRegionTrack - Stub for Biomart-based gene region queries.

This is a placeholder that mirrors Gviz's BiomartGeneRegionTrack API but
raises NotImplementedError on draw, since it requires a live Biomart
connection.

Use GeneRegionTrack with a local GTF/GFF file instead.
"""

from typing import Optional, Dict, Any, Union, List

from ._gene_region import GeneRegionTrack
from ._base import GenomicInterval


class BiomartGeneRegionTrack(GeneRegionTrack):
    """Stub for a gene region track that queries a Biomart database.

    This class accepts the same constructor parameters as the Gviz
    ``BiomartGeneRegionTrack`` but does **not** implement live Biomart
    queries.  Calling ``draw()`` raises ``NotImplementedError``.

    To achieve equivalent functionality, download the relevant GTF/GFF
    file from Ensembl Biomart and use :class:`GeneRegionTrack` directly::

        from geneview.genometracks import GeneRegionTrack, read_gff

        data = read_gff("genes.gtf")
        track = GeneRegionTrack(data, collapse_transcripts="longest")

    Parameters
    ----------
    genome : str, optional
        Genome assembly (e.g. ``"hg38"``).
    chromosome : str, optional
        Chromosome to query.
    start : int, optional
        Start position.
    end : int, optional
        End position.
    symbol : str or list of str, optional
        Gene symbol(s) to query.
    gene : str or list of str, optional
        Ensembl gene ID(s).
    transcript : str or list of str, optional
        Ensembl transcript ID(s).
    entrez : str or list of str, optional
        Entrez gene ID(s).
    filter : dict, optional
        Additional Biomart filters.
    biomart : str, optional
        Biomart dataset name.
    name : str
        Track name.
    height : float
        Relative track height.
    display_params : dict, optional
        Display parameters.
    """

    def __init__(
        self,
        genome: Optional[str] = None,
        chromosome: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        symbol: Optional[Union[str, List[str]]] = None,
        gene: Optional[Union[str, List[str]]] = None,
        transcript: Optional[Union[str, List[str]]] = None,
        entrez: Optional[Union[str, List[str]]] = None,
        filter: Optional[Dict[str, Any]] = None,
        biomart: Optional[str] = None,
        name: str = "BiomartGeneRegion",
        height: float = 1.5,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(data=None, name=name, height=height,
                         display_params=display_params)
        self.genome = genome
        self.chromosome = chromosome
        self._start = start
        self._end = end
        self.symbol = symbol
        self.gene = gene
        self.transcript = transcript
        self.entrez = entrez
        self.filter = filter
        self.biomart = biomart

    def draw(self, ax, region: GenomicInterval) -> None:
        """Raise NotImplementedError — Biomart queries are not implemented.

        Use GeneRegionTrack with a local GTF/GFF file instead.
        """
        raise NotImplementedError(
            "BiomartGeneRegionTrack requires a Biomart connection. "
            "Use GeneRegionTrack with a local GTF/GFF file instead. "
            "Example:\n"
            "    from geneview.genometracks import GeneRegionTrack, read_gff\n"
            "    data = read_gff('genes.gtf')\n"
            "    track = GeneRegionTrack(data)"
        )

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the requested region if chromosome/start/end were provided."""
        if self.chromosome and self._start is not None and self._end is not None:
            return GenomicInterval(self.chromosome, self._start, self._end)
        return None
