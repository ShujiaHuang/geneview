## Geneview Tutorials

This directory contains Jupyter notebook tutorials for **geneview** — a Python package for genomics data visualization.

### Available Tutorials

| Tutorial | Description | Notebook |
|----------|-------------|----------|
| **GWAS Plots** | Manhattan plots and Q-Q plots for genome-wide association studies | [gwas_plot.ipynb](./gwas_plot.ipynb) |
| **Admixture** | Population structure visualization from ADMIXTURE output | [admixture.ipynb](./admixture.ipynb) |
| **Venn Diagrams** | Set intersection diagrams for 2–6 datasets | [venn.ipynb](./venn.ipynb) |
| **Color Palettes** | Color schemes and palette customization for genomics figures | [palettes.ipynb](./palettes.ipynb) |

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
