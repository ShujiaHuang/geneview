#!/bin/bash
# Example CLI usage for geneview subcommands.
#
# Before running, make sure geneview is installed:
#   pip install -e .
#
# Author: Shujia Huang

set -e

# --- Setup: Generate test data if needed ---
python examples/scripts/generate_test_data.py

DATA_DIR="examples/data"
OUT_DIR="examples/figures/cli"
mkdir -p "$OUT_DIR"

# --- Manhattan plot ---
echo "=== Creating Manhattan plot ==="
geneview manhattan \
    -i "$DATA_DIR/gwas_example.assoc" \
    -o "$OUT_DIR/manhattan.png" \
    --title "Synthetic GWAS" \
    --sign-marker-p 1e-6 \
    --annotate-topsnp \
    --figsize 12 4

# --- Q-Q plot ---
echo "=== Creating Q-Q plot ==="
geneview qq \
    -i "$DATA_DIR/gwas_example.assoc" \
    -o "$OUT_DIR/qq.png" \
    --title "Synthetic GWAS " \
    --figsize 6 6

# --- Venn diagram (3 sets) ---
echo "=== Creating Venn diagram ==="
geneview venn \
    -i "$DATA_DIR/venn_sets/gene_list_1.txt" \
       "$DATA_DIR/venn_sets/gene_list_2.txt" \
       "$DATA_DIR/venn_sets/gene_list_3.txt" \
    -o "$OUT_DIR/venn3.png" \
    --names "Dataset A" "Dataset B" "Dataset C" \
    --palette "viridis" \
    --legend-use-petal-color

# --- Admixture plot ---
echo "=== Creating Admixture plot ==="
geneview admixture \
    -i "$DATA_DIR/admixture_example/synthetic.Q" \
    -p "$DATA_DIR/admixture_example/synthetic_pop.txt" \
    -o "$OUT_DIR/admixture.png" \
    --palette "Set1" \
    --edgewidth 2.0 \
    --set-xticklabel-top

# --- Genome tracks (basic) ---
echo "=== Creating Genome tracks (basic) ==="
geneview tracks \
    --region chr7:26490000-26720000 \
    --ideogram \
    -a "$DATA_DIR/genome_tracks/cpg_islands.bed" \
    -g "$DATA_DIR/genome_tracks/gene_models.gtf" \
    -d "$DATA_DIR/genome_tracks/coverage.bedgraph" \
    -o "$OUT_DIR/tracks_basic.png" \
    --figsize 14 8

# --- Genome tracks (with highlights and custom appearance) ---
echo "=== Creating Genome tracks (with highlights) ==="
geneview tracks \
    --region chr7:26M-27M \
    --ideogram \
    -a "$DATA_DIR/genome_tracks/annotations.bed" --annotation-shape ellipse \
    -g "$DATA_DIR/genome_tracks/gene_models.gtf" --collapse-transcripts longest \
    -d "$DATA_DIR/genome_tracks/coverage.bedgraph" --data-type line --data-color blue \
    --highlight "$DATA_DIR/genome_tracks/annotations.bed" --highlight-fill yellow \
    -o "$OUT_DIR/tracks_custom.png" \
    --figsize 14 10

echo ""
echo "=== All figures saved to $OUT_DIR/ ==="
ls -la "$OUT_DIR/"
