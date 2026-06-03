"""Generate synthetic example data for geneview CLI testing and demonstrations.

This script creates small, self-contained test datasets that can be used with
the geneview CLI subcommands without requiring network access.

Author: Shujia Huang
"""
import os
import numpy as np
import pandas as pd


def generate_gwas_data(output_path, n_chroms=5, snps_per_chrom=200, seed=42):
    """Generate a synthetic GWAS association file (tab-delimited).

    Parameters
    ----------
    output_path : str
        Path to the output file.
    n_chroms : int
        Number of chromosomes to simulate.
    snps_per_chrom : int
        Number of SNPs per chromosome.
    seed : int
        Random seed for reproducibility.
    """
    rng = np.random.RandomState(seed)
    rows = []
    for chrom in range(1, n_chroms + 1):
        chrom_name = "chr%d" % chrom
        positions = np.sort(rng.randint(1, 5_000_000, size=snps_per_chrom))
        # Most p-values are uniform, but inject a few significant ones
        pvalues = rng.uniform(0.01, 1.0, size=snps_per_chrom)
        sign_indices = rng.choice(snps_per_chrom, size=5, replace=False)
        pvalues[sign_indices] = rng.uniform(1e-12, 1e-6, size=5)

        for i, (pos, pv) in enumerate(zip(positions, pvalues)):
            rows.append({
                "#CHROM": chrom_name,
                "POS": int(pos),
                "P": pv,
                "ID": "rs%d_%d" % (chrom, i),
            })
    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep="\t", index=False)
    print("[INFO] GWAS data written to %s (%d rows)" % (output_path, len(df)))


def generate_venn_files(output_dir, n_files=3, n_items=200, overlap=50, seed=42):
    """Generate synthetic gene list files for Venn diagram testing.

    Parameters
    ----------
    output_dir : str
        Directory to write the output files.
    n_files : int
        Number of gene list files to create (2-6).
    n_items : int
        Number of items in each list.
    overlap : int
        Number of shared items between consecutive lists.
    seed : int
        Random seed for reproducibility.
    """
    rng = np.random.RandomState(seed)
    os.makedirs(output_dir, exist_ok=True)

    all_genes = ["GENE_%04d" % i for i in range(n_items * n_files)]
    rng.shuffle(all_genes)

    shared_core = all_genes[:overlap]
    offset = overlap
    file_paths = []
    for i in range(n_files):
        unique_genes = all_genes[offset:offset + n_items - overlap]
        offset += n_items - overlap
        genes = shared_core + list(unique_genes)
        rng.shuffle(genes)

        path = os.path.join(output_dir, "gene_list_%d.txt" % (i + 1))
        with open(path, "w") as f:
            f.write("\n".join(genes) + "\n")
        file_paths.append(path)
        print("[INFO] Venn file written to %s (%d genes)" % (path, len(genes)))

    return file_paths


def generate_admixture_data(output_dir, n_groups=4, samples_per_group=30,
                            k=3, seed=42):
    """Generate synthetic ADMIXTURE output files.

    Parameters
    ----------
    output_dir : str
        Directory to write the output files.
    n_groups : int
        Number of population groups.
    samples_per_group : int
        Number of samples per group.
    k : int
        Number of ancestral components.
    seed : int
        Random seed for reproducibility.
    """
    rng = np.random.RandomState(seed)
    os.makedirs(output_dir, exist_ok=True)

    q_rows = []
    pop_labels = []
    group_names = ["POP%d" % (i + 1) for i in range(n_groups)]

    for g_idx, g in enumerate(group_names):
        for _ in range(samples_per_group):
            # Each group has a dominant component
            proportions = rng.dirichlet(np.ones(k) * 0.5)
            proportions[g_idx % k] += 0.5  # boost one component
            proportions /= proportions.sum()  # re-normalize
            q_rows.append(proportions)
            pop_labels.append(g)

    # Write Q file
    q_path = os.path.join(output_dir, "synthetic.Q")
    with open(q_path, "w") as f:
        for row in q_rows:
            f.write(" ".join("%.6f" % v for v in row) + "\n")
    print("[INFO] Admixture Q file written to %s" % q_path)

    # Write population info file
    pop_path = os.path.join(output_dir, "synthetic_pop.txt")
    with open(pop_path, "w") as f:
        f.write("\n".join(pop_labels) + "\n")
    print("[INFO] Population info written to %s" % pop_path)

    return q_path, pop_path


if __name__ == "__main__":
    data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
    os.makedirs(data_dir, exist_ok=True)

    generate_gwas_data(os.path.join(data_dir, "gwas_example.assoc"))
    generate_venn_files(os.path.join(data_dir, "venn_sets"), n_files=3)
    generate_admixture_data(os.path.join(data_dir, "admixture_example"))
