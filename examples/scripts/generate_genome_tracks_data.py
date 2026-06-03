"""Generate synthetic example data for geneview genome tracks demonstrations.

This script creates small, self-contained genomic data files that can be used
with the genome tracks example scripts and documentation. Data is modeled after
the Gviz R package's test datasets (chr7 region ~26.5-26.8 Mb).

Author: geneview contributors
"""
import os
import shutil
import numpy as np


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "genome_tracks")
GVIZ_EXT = os.path.join(SCRIPT_DIR, "..", "..", "learn_from", "Gviz", "inst", "extdata")


def ensure_dirs():
    os.makedirs(DATA_DIR, exist_ok=True)
    gviz_out = os.path.join(DATA_DIR, "gviz_samples")
    os.makedirs(gviz_out, exist_ok=True)
    return gviz_out


# ---------------------------------------------------------------------------
# 1. CpG islands (BED6) -- mimics Gviz cpgIslands on chr7
# ---------------------------------------------------------------------------
def generate_cpg_islands(output_dir):
    """Generate CpG island BED file on chr7:26,500,000-26,800,000."""
    path = os.path.join(output_dir, "cpg_islands.bed")
    features = [
        ("chr7", 26500100, 26500800, "CpG_1",  0, "+"),
        ("chr7", 26520000, 26520500, "CpG_2",  0, "+"),
        ("chr7", 26560200, 26560900, "CpG_3",  0, "-"),
        ("chr7", 26580000, 26580350, "CpG_4",  0, "+"),
        ("chr7", 26610500, 26611200, "CpG_5",  0, "-"),
        ("chr7", 26660000, 26660600, "CpG_6",  0, "+"),
        ("chr7", 26700800, 26701500, "CpG_7",  0, "+"),
        ("chr7", 26730000, 26730400, "CpG_8",  0, "-"),
        ("chr7", 26770200, 26770700, "CpG_9",  0, "+"),
    ]
    with open(path, "w") as f:
        for chrom, s, e, name, score, strand in features:
            f.write(f"{chrom}\t{s}\t{e}\t{name}\t{score}\t{strand}\n")
    print(f"[INFO] CpG islands written to {path} ({len(features)} features)")
    return path


# ---------------------------------------------------------------------------
# 2. Gene models (GTF) -- 3 genes, multiple transcripts, on chr7
# ---------------------------------------------------------------------------
def generate_gene_models(output_dir):
    """Generate gene model GTF file with 3 genes on chr7."""
    path = os.path.join(output_dir, "gene_models.gtf")
    lines = []

    def _add(chrom, feat, start, end, strand, gene_id, tx_id, exon_num, gene_name):
        attrs = (
            f'gene_id "{gene_id}"; transcript_id "{tx_id}"; '
            f'exon_number "{exon_num}"; gene_name "{gene_name}";'
        )
        lines.append(f"{chrom}\tgeneview\t{feat}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}")

    # Gene 1: BRCA-like gene on + strand, 2 transcripts
    gene1 = ("chr7", "+", "ENSG00001", "BRCA_L1")
    # Transcript 1a (3 exons)
    _add(*gene1[:2], 26500500, 26500800, "+", gene1[2], "ENST00001a", 1, gene1[3])
    _add(*gene1[:2], 26500500, 26500650, "+", gene1[2], "ENST00001a", 1, gene1[3])  # 5UTR
    _add(*gene1[:2], 26500651, 26500800, "+", gene1[2], "ENST00001a", 1, gene1[3])  # CDS
    _add(*gene1[:2], 26510000, 26510400, "+", gene1[2], "ENST00001a", 2, gene1[3])
    _add(*gene1[:2], 26510000, 26510400, "+", gene1[2], "ENST00001a", 2, gene1[3])  # CDS
    _add(*gene1[:2], 26520000, 26520600, "+", gene1[2], "ENST00001a", 3, gene1[3])
    _add(*gene1[:2], 26520000, 26520400, "+", gene1[2], "ENST00001a", 3, gene1[3])  # CDS
    _add(*gene1[:2], 26520401, 26520600, "+", gene1[2], "ENST00001a", 3, gene1[3])  # 3UTR
    # Transcript 1b (2 exons, shorter)
    _add(*gene1[:2], 26500500, 26500800, "+", gene1[2], "ENST00001b", 1, gene1[3])
    _add(*gene1[:2], 26500651, 26500800, "+", gene1[2], "ENST00001b", 1, gene1[3])  # CDS
    _add(*gene1[:2], 26510000, 26510400, "+", gene1[2], "ENST00001b", 2, gene1[3])
    _add(*gene1[:2], 26510000, 26510300, "+", gene1[2], "ENST00001b", 2, gene1[3])  # CDS
    _add(*gene1[:2], 26510301, 26510400, "+", gene1[2], "ENST00001b", 2, gene1[3])  # 3UTR

    # Gene 2: Kinase-like gene on - strand, 2 transcripts
    gene2 = ("chr7", "-", "ENSG00002", "KIN_L1")
    # Transcript 2a (4 exons)
    _add(*gene2[:2], 26600000, 26600500, "-", gene2[2], "ENST00002a", 1, gene2[3])
    _add(*gene2[:2], 26600000, 26600200, "-", gene2[2], "ENST00002a", 1, gene2[3])  # 5UTR
    _add(*gene2[:2], 26600201, 26600500, "-", gene2[2], "ENST00002a", 1, gene2[3])  # CDS
    _add(*gene2[:2], 26620000, 26620600, "-", gene2[2], "ENST00002a", 2, gene2[3])
    _add(*gene2[:2], 26620000, 26620600, "-", gene2[2], "ENST00002a", 2, gene2[3])  # CDS
    _add(*gene2[:2], 26640000, 26640350, "-", gene2[2], "ENST00002a", 3, gene2[3])
    _add(*gene2[:2], 26640000, 26640350, "-", gene2[2], "ENST00002a", 3, gene2[3])  # CDS
    _add(*gene2[:2], 26660000, 26660800, "-", gene2[2], "ENST00002a", 4, gene2[3])
    _add(*gene2[:2], 26660000, 26660500, "-", gene2[2], "ENST00002a", 4, gene2[3])  # CDS
    _add(*gene2[:2], 26660501, 26660800, "-", gene2[2], "ENST00002a", 4, gene2[3])  # 3UTR
    # Transcript 2b (3 exons)
    _add(*gene2[:2], 26600000, 26600500, "-", gene2[2], "ENST00002b", 1, gene2[3])
    _add(*gene2[:2], 26600201, 26600500, "-", gene2[2], "ENST00002b", 1, gene2[3])  # CDS
    _add(*gene2[:2], 26620000, 26620600, "-", gene2[2], "ENST00002b", 2, gene2[3])
    _add(*gene2[:2], 26620000, 26620600, "-", gene2[2], "ENST00002b", 2, gene2[3])  # CDS
    _add(*gene2[:2], 26640000, 26640350, "-", gene2[2], "ENST00002b", 3, gene2[3])
    _add(*gene2[:2], 26640000, 26640250, "-", gene2[2], "ENST00002b", 3, gene2[3])  # CDS
    _add(*gene2[:2], 26640251, 26640350, "-", gene2[2], "ENST00002b", 3, gene2[3])  # 3UTR

    # Gene 3: Small single-exon gene on + strand
    gene3 = ("chr7", "+", "ENSG00003", "SNRNP_L1")
    _add(*gene3[:2], 26700000, 26700600, "+", gene3[2], "ENST00003a", 1, gene3[3])
    _add(*gene3[:2], 26700000, 26700150, "+", gene3[2], "ENST00003a", 1, gene3[3])  # 5UTR
    _add(*gene3[:2], 26700151, 26700450, "+", gene3[2], "ENST00003a", 1, gene3[3])  # CDS
    _add(*gene3[:2], 26700451, 26700600, "+", gene3[2], "ENST00003a", 1, gene3[3])  # 3UTR

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[INFO] Gene models written to {path} ({len(lines)} features)")
    return path


# ---------------------------------------------------------------------------
# 3. Coverage data (bedGraph) -- numeric signal on chr7
# ---------------------------------------------------------------------------
def generate_coverage(output_dir, seed=42):
    """Generate bedGraph coverage data on chr7:26,500,000-26,800,000."""
    rng = np.random.RandomState(seed)
    path = os.path.join(output_dir, "coverage.bedgraph")
    start = 26_500_000
    end = 26_800_000
    n_bins = 100
    bin_size = (end - start) // n_bins
    values = rng.normal(5, 2, n_bins).cumsum() / 10  # random walk

    with open(path, "w") as f:
        for i in range(n_bins):
            s = start + i * bin_size
            e = s + bin_size
            v = values[i]
            f.write(f"chr7\t{s}\t{e}\t{v:.4f}\n")
    print(f"[INFO] Coverage data written to {path} ({n_bins} bins)")
    return path


# ---------------------------------------------------------------------------
# 4. Annotation features (BED6) -- simple annotation on chr7
# ---------------------------------------------------------------------------
def generate_annotations(output_dir):
    """Generate annotation BED file with strand info on chr7."""
    path = os.path.join(output_dir, "annotations.bed")
    features = [
        ("chr7", 26505000, 26508000, "Promoter_A",  0, "+"),
        ("chr7", 26530000, 26535000, "Enhancer_1",  0, "+"),
        ("chr7", 26570000, 26575000, "CTCF_site",   0, "-"),
        ("chr7", 26610000, 26615000, "Promoter_B",  0, "-"),
        ("chr7", 26650000, 26655000, "Enhancer_2",  0, "+"),
        ("chr7", 26700000, 26705000, "Promoter_C",  0, "+"),
        ("chr7", 26740000, 26745000, "DNase_HSS",   0, "-"),
        ("chr7", 26780000, 26785000, "Enhancer_3",  0, "+"),
    ]
    with open(path, "w") as f:
        for chrom, s, e, name, score, strand in features:
            f.write(f"{chrom}\t{s}\t{e}\t{name}\t{score}\t{strand}\n")
    print(f"[INFO] Annotations written to {path} ({len(features)} features)")
    return path


# ---------------------------------------------------------------------------
# 5. Multi-sample data (bedGraph-like with extra columns) for DataTrack
# ---------------------------------------------------------------------------
def generate_multi_sample(output_dir, seed=42):
    """Generate multi-sample numeric data as a TSV for DataTrack grouping."""
    rng = np.random.RandomState(seed)
    path = os.path.join(output_dir, "multi_sample.tsv")
    start = 26_500_000
    end = 26_800_000
    n_bins = 80
    bin_size = (end - start) // n_bins

    header = "chrom\tstart\tend\tsample_ctrl_1\tsample_ctrl_2\tsample_ctrl_3\tsample_treat_1\tsample_treat_2\tsample_treat_3"
    with open(path, "w") as f:
        f.write(header + "\n")
        for i in range(n_bins):
            s = start + i * bin_size
            e = s + bin_size
            base = np.sin(i / n_bins * 4 * np.pi) * 5 + 10
            c1 = base + rng.normal(0, 0.5)
            c2 = base + rng.normal(0, 0.5)
            c3 = base + rng.normal(0, 0.5)
            t1 = base + 3 + rng.normal(0, 0.8)   # treatment shifts signal up
            t2 = base + 3 + rng.normal(0, 0.8)
            t3 = base + 3 + rng.normal(0, 0.8)
            f.write(f"chr7\t{s}\t{e}\t{c1:.3f}\t{c2:.3f}\t{c3:.3f}\t{t1:.3f}\t{t2:.3f}\t{t3:.3f}\n")
    print(f"[INFO] Multi-sample data written to {path} ({n_bins} bins x 6 samples)")
    return path


# ---------------------------------------------------------------------------
# 6. Copy Gviz sample test files
# ---------------------------------------------------------------------------
def copy_gviz_samples(gviz_out_dir):
    """Copy real Gviz test files for users who want to try loading real data."""
    filenames = ["test.bed", "test.bedGraph", "test.gtf", "test.gff3"]
    copied = 0
    for fname in filenames:
        src = os.path.join(GVIZ_EXT, fname)
        dst = os.path.join(gviz_out_dir, fname)
        if os.path.exists(src):
            shutil.copy2(src, dst)
            copied += 1
            print(f"[INFO] Copied {fname} -> {dst}")
        else:
            print(f"[WARN] Gviz sample file not found: {src}")
    print(f"[INFO] {copied}/{len(filenames)} Gviz sample files copied")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    gviz_out = ensure_dirs()
    generate_cpg_islands(DATA_DIR)
    generate_gene_models(DATA_DIR)
    generate_coverage(DATA_DIR)
    generate_annotations(DATA_DIR)
    generate_multi_sample(DATA_DIR)
    copy_gviz_samples(gviz_out)
    print("\n[DONE] All genome tracks example data generated successfully.")
