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

echo ""
echo "=== All figures saved to $OUT_DIR/ ==="
ls -la "$OUT_DIR/"
