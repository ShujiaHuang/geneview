## Geneview Tutorials

This directory contains Jupyter notebook tutorials for **geneview** — a Python package for genomics data visualization.

### Available Tutorials

| Tutorial | Description | Notebook |
|----------|-------------|----------|
| **GWAS Plots** | Manhattan plots and Q-Q plots for genome-wide association studies | [gwas_plot.ipynb](./gwas_plot.ipynb) |
| **Admixture** | Population structure visualization from ADMIXTURE output | [admixture.ipynb](./admixture.ipynb) |
| **Venn Diagrams** | Set intersection diagrams for 2–6 datasets | [venn.ipynb](./venn.ipynb) |
| **Color Palettes** | Color schemes and palette customization for genomics figures | [palettes.ipynb](./palettes.ipynb) |
| **Genome Tracks** | Gviz-style genome track visualization (axis, annotations, gene models, data tracks, highlights) | [genome_tracks.ipynb](./genome_tracks.ipynb) |

### Additional Documentation

| Document | Description |
|----------|-------------|
| [User Guide](../user_guide.md) | Comprehensive guide covering all geneview modules (GWAS, Venn, Admixture, Karyotype, Genome Tracks, CLI, and more) |
| [Genome Tracks Guide](../genome_tracks_guide.md) | Detailed Gviz-style guide for the `geneview.genometracks` module — track types, display parameters, file I/O, and complete examples |

### Example Scripts

Runnable Python scripts demonstrating genome tracks features are available in [`examples/scripts/`](../../examples/scripts/):

| Script | Description |
|--------|-------------|
| `genome_tracks_basic.py` | Basic GenomeAxisTrack + AnnotationTrack from BED |
| `genome_tracks_gene_region.py` | Gene models from GTF with collapsing modes |
| `genome_tracks_data.py` | DataTrack plot types: line, histogram, polygon, heatmap, points, mountain, gradient |
| `genome_tracks_highlight.py` | Cross-track highlight regions |
| `genome_tracks_comprehensive.py` | Full showcase with all track types combined |

Run any script with: `python examples/scripts/<script_name>.py`

Synthetic example data can be regenerated with: `python examples/scripts/generate_genome_tracks_data.py`

### Set Up the Environment for Jupyter Notebook

```bash
$ pip install jupyter        # Skip if already installed
$ pip install ipykernel
$ pip install geneview       # Install geneview and its dependencies
$ python -m ipykernel install --user --name venv --display-name "Python3(geneview)"
```

### Launch a Tutorial

```bash
$ jupyter notebook
```

Then navigate to any `.ipynb` file above and select the **Python3(geneview)** kernel to get started.

> **Note:** The `genome_tracks.ipynb` notebook does not yet exist as a Jupyter file. Use the [genome tracks example scripts](../../examples/scripts/) and the [Genome Tracks Guide](../genome_tracks_guide.md) in the meantime.

### Dataset

The tutorials use example datasets hosted in [geneview-data](https://github.com/ShujiaHuang/geneview-data).
You can load them programmatically:

```python
import geneview as gv

# List available datasets
names = gv.get_dataset_names()

# Load a dataset by name
df = gv.load_dataset("gwas")
```
