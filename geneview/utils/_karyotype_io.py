"""
Shared reader/normalizer for cytoband (karyotype) data.

Both :func:`geneview.karyotype.karyoplot` and
:class:`geneview.genometracks.IdeogramTrack` need to turn assorted cytoband
inputs (a UCSC/Gviz karyotype file, a DataFrame, or an array of rows) into a
single canonical schema.  Centralizing that here keeps the two rendering paths
in sync.
"""
import pandas as pd

# Canonical output schema (UCSC cytoBandIdeo style).
CYTOBAND_COLUMNS = ["chrom", "chromStart", "chromEnd", "name", "gieStain"]

# Accepted alias -> canonical column name (compared case-insensitively).
_COLUMN_ALIASES = {
    "chrom": "chrom", "chromosome": "chrom", "chr": "chrom",
    "chromstart": "chromStart", "start": "chromStart",
    "chromend": "chromEnd", "end": "chromEnd",
    "name": "name", "band": "name", "band_name": "name",
    "giestain": "gieStain", "gie_stain": "gieStain",
    "stain": "gieStain", "type": "gieStain",
}


def _peek_first_char(path):
    """Return the first character of ``path`` (local file or URL)."""
    try:
        with open(path) as fh:
            return fh.read(1)
    except OSError:
        # Not a local file (e.g. an http/S3 URL); peek via urllib.
        from urllib.request import urlopen
        with urlopen(path) as resp:
            return resp.read(1).decode("utf-8", "ignore")


def _read_cytoband_file(path):
    """Read a cytoband table from a file path or URL into a DataFrame."""
    if _peek_first_char(path) == "#":
        # UCSC/Gviz style: the header line starts with '#', e.g.
        # ``#chrom  chromStart  chromEnd  name  gieStain``.
        return pd.read_table(
            path, comment="#", header=None, names=CYTOBAND_COLUMNS,
        )
    # Plain TSV with a regular header row.
    return pd.read_table(path)


def _normalize_columns(df):
    """Map recognized/positional columns onto the canonical schema."""
    if set(CYTOBAND_COLUMNS).issubset(df.columns):
        return df

    col_map = {
        c: _COLUMN_ALIASES[str(c).lower().strip()]
        for c in df.columns
        if str(c).lower().strip() in _COLUMN_ALIASES
    }
    df = df.rename(columns=col_map)

    # Fall back to positional assignment when names are unknown but the column
    # count matches (e.g. list-of-lists input).
    if not set(CYTOBAND_COLUMNS).issubset(df.columns) and \
            df.shape[1] == len(CYTOBAND_COLUMNS):
        df = df.copy()
        df.columns = CYTOBAND_COLUMNS
    return df


def read_cytoband(data):
    """Read cytoband/karyotype data into a canonical DataFrame.

    Parameters
    ----------
    data : str, pandas.DataFrame, or array-like
        Either a path/URL to a tab-separated cytoband file (UCSC ``#chrom``
        comment header or a plain header row are both accepted), a DataFrame,
        or an array-like of rows in ``chrom, chromStart, chromEnd, name,
        gieStain`` order.

    Returns
    -------
    pandas.DataFrame
        A DataFrame whose leading columns are ``chrom, chromStart, chromEnd,
        name, gieStain`` (with ``chromStart``/``chromEnd`` cast to ``int``).
        Any extra input columns are preserved after the canonical ones.
    """
    if isinstance(data, str):
        df = _read_cytoband_file(data)
    elif isinstance(data, pd.DataFrame):
        df = data.copy()
    else:
        df = pd.DataFrame(list(data))

    df = _normalize_columns(df)

    missing = set(CYTOBAND_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required cytoband columns: {sorted(missing)}")

    df["chromStart"] = df["chromStart"].astype(int)
    df["chromEnd"] = df["chromEnd"].astype(int)

    extra = [c for c in df.columns if c not in CYTOBAND_COLUMNS]
    return df[CYTOBAND_COLUMNS + extra]
