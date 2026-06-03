# Installation Guide

<cite>
**Referenced Files in This Document**
- [setup.py](file://setup.py)
- [requirements.txt](file://requirements.txt)
- [README.md](file://README.md)
- [.github/workflows/main.yml](file://.github/workflows/main.yml)
- [ci/utils.txt](file://ci/utils.txt)
- [docs/tutorial/README.md](file://docs/tutorial/README.md)
- [geneview/cli/main.py](file://geneview/cli/main.py)
- [release.py](file://release.py)
</cite>

## Update Summary
**Changes Made**
- Updated Python version requirements from 3.7+ to 3.9+ to reflect current CI workflow and setup.py classifiers
- Removed mention of Python 3.8 support as it's no longer tested in CI
- Updated system requirements section to reflect new Python version constraints
- Updated verification steps to align with new supported Python versions

## Table of Contents
1. [Introduction](#introduction)
2. [System Requirements](#system-requirements)
3. [Installation Methods](#installation-methods)
4. [Core Dependencies](#core-dependencies)
5. [Optional Dependencies](#optional-dependencies)
6. [Virtual Environment Setup](#virtual-environment-setup)
7. [Verification Steps](#verification-steps)
8. [Operating System Specific Guidance](#operating-system-specific-guidance)
9. [Conda Users Guide](#conda-users-guide)
10. [Command-Line Interface (CLI) Installation](#command-line-interface-cli-installation)
11. [Troubleshooting Common Issues](#troubleshooting-common-issues)
12. [Performance Considerations](#performance-considerations)
13. [Conclusion](#conclusion)

## Introduction
GeneView is a Python package designed for genomics data visualization. It provides high-level abstractions for creating attractive and informative genomics graphics, built on top of the PyData stack including NumPy, pandas, matplotlib, scipy, and seaborn. **Version 0.3.0** introduces a powerful command-line interface (CLI) alongside the traditional Python API, allowing users to create publication-quality plots directly from the terminal without writing any Python code.

## System Requirements
- **Python Version**: Python 3.9, 3.10, or 3.12 (Python 2.x and Python 3.8 are not supported)
- **Operating Systems**: POSIX-compliant systems including Unix, Linux, and macOS
- **Hardware**: Standard desktop or server hardware sufficient for typical genomics data visualization tasks

**Updated** Updated Python version requirements to reflect current CI workflow and setup.py classifiers, removing support for Python 3.7 and 3.8.

**Section sources**
- [setup.py:57-68](file://setup.py#L57-L68)
- [setup.py:17](file://setup.py#L17)
- [.github/workflows/main.yml:22](file://.github/workflows/main.yml#L22)

## Installation Methods

### Method 1: Using pip (Recommended)
The simplest and most straightforward installation method:

```bash
pip install geneview
```

This command installs GeneView along with all core dependencies automatically. The installation includes both the Python library and the command-line interface.

**Section sources**
- [README.md:30-36](file://README.md#L30-L36)

### Method 2: From Source Code
For development or customization purposes:

```bash
git clone https://github.com/ShujiaHuang/geneview.git
cd geneview
pip install .
```

### Method 3: Development Installation
For contributing to the project or using unreleased features:

```bash
pip install -e .
```

### Method 4: Installing with Optional Dependencies
To include optional dependencies for enhanced functionality:

```bash
pip install geneview[all]
```

**Section sources**
- [.github/workflows/main.yml:51](file://.github/workflows/main.yml#L51)

## Core Dependencies
GeneView requires the following essential dependencies:

| Package | Purpose | Version Requirements |
|---------|---------|---------------------|
| **NumPy** | Numerical computing foundation | Essential for array operations |
| **SciPy** | Scientific computing and statistics | Required for advanced statistical functions |
| **pandas** | Data manipulation and analysis | Essential for DataFrame operations |
| **matplotlib** | 2D plotting library | Core visualization engine |
| **seaborn** | Statistical data visualization | Enhanced plotting aesthetics |

These dependencies are automatically installed when you install GeneView via pip.

**Section sources**
- [setup.py:44-50](file://setup.py#L44-L50)
- [requirements.txt:1-5](file://requirements.txt#L1-L5)
- [README.md:463-471](file://README.md#L463-L471)

## Optional Dependencies
While not required for basic functionality, these optional packages enhance specific features:

### statsmodels
- **Purpose**: Advanced statistical modeling and econometrics functions
- **When to install**: Required for certain statistical analyses and advanced plotting features
- **Installation**: `pip install statsmodels`

### fastcluster
- **Purpose**: High-performance hierarchical clustering algorithms
- **When to install**: Recommended for large-scale clustering operations to improve performance
- **Installation**: `pip install fastcluster`

**Section sources**
- [README.md:463-471](file://README.md#L463-L471)
- [.github/workflows/main.yml:51](file://.github/workflows/main.yml#L51)

## Virtual Environment Setup
It is highly recommended to use virtual environments to manage GeneView dependencies separately from your system Python installation.

### Using venv (Python 3.3+)
```bash
# Create virtual environment
python -m venv geneview_env

# Activate environment (Linux/macOS)
source geneview_env/bin/activate

# Activate environment (Windows)
geneview_env\Scripts\activate

# Install GeneView
pip install geneview
```

### Using conda
```bash
# Create conda environment
conda create -n geneview_env python=3.9

# Activate environment
conda activate geneview_env

# Install GeneView
pip install geneview
```

**Section sources**
- [docs/tutorial/README.md:3-6](file://docs/tutorial/README.md#L3-L6)

## Verification Steps
After installation, verify that GeneView is working correctly:

### Basic Import Test
```python
import geneview as gv
print(gv.__version__)
```

### Simple Function Test
```python
import geneview as gv
import matplotlib.pyplot as plt

# Load sample dataset and create a basic plot
df = gv.load_dataset("gwas")
ax = gv.manhattanplot(df)
plt.show()
```

### Dependency Verification
```python
import sys
print(f"Python version: {sys.version}")

# Verify core dependencies
import numpy, pandas, matplotlib, scipy, seaborn
print(f"NumPy: {numpy.__version__}")
print(f"pandas: {pandas.__version__}")
print(f"matplotlib: {matplotlib.__version__}")
print(f"SciPy: {scipy.__version__}")
print(f"seaborn: {seaborn.__version__}")
```

### CLI Verification
```bash
# Check if CLI is available
geneview --help

# Check version
geneview --version
```

**Section sources**
- [geneview/cli/main.py:50-53](file://geneview/cli/main.py#L50-L53)

## Operating System Specific Guidance

### Linux (Ubuntu/Debian)
Ensure you have the required system packages:
```bash
sudo apt-get update
sudo apt-get install python3-dev python3-pip libfreetype6-dev
pip install geneview
```

### macOS
For macOS users, especially those with Apple Silicon:
```bash
# Install Xcode command line tools if not present
xcode-select --install

# Install GeneView
pip install geneview
```

### Windows
Windows users can install GeneView using either:
- **Standard Python**: `pip install geneview`
- **WSL (Windows Subsystem for Linux)**: Recommended for optimal compatibility
- **Docker**: Containerized deployment option

**Section sources**
- [setup.py:61-68](file://setup.py#L61-L68)

## Conda Users Guide
Conda users have several options for installing GeneView:

### Option 1: Pure pip within conda environment
```bash
conda create -n geneview python=3.9
conda activate geneview
pip install geneview
```

### Option 2: Bioconda channel (if available)
```bash
conda install -c bioconda geneview
```

### Option 3: Environment isolation
```bash
# Create dedicated environment
conda create -n geneview_env python=3.9

# Activate environment
conda activate geneview_env

# Install dependencies
conda install numpy pandas matplotlib scipy seaborn

# Install GeneView
pip install geneview
```

**Section sources**
- [docs/tutorial/README.md:3-6](file://docs/tutorial/README.md#L3-L6)

## Command-Line Interface (CLI) Installation

### Installing the CLI
After installing GeneView via pip, the command-line interface becomes immediately available. The CLI provides four main subcommands for different types of genomics visualizations:

```bash
# View all available subcommands
geneview --help

# Output includes:
# subcommands:
#   manhattan    Create a Manhattan plot from GWAS association results.
#   qq           Create a Q-Q plot from GWAS association results.
#   venn         Create a Venn diagram from 2-6 input files.
#   admixture    Create an Admixture plot from ADMIXTURE .Q output.
```

### CLI Usage Examples

#### Manhattan Plot
```bash
# Basic Manhattan plot from PLINK2.x association results
geneview manhattan -i gwas_results.assoc -o manhattan.png

# Advanced Manhattan plot with custom options
geneview manhattan -i gwas_results.assoc -o manhattan.png \
    --title "My GWAS" \
    --sign-marker-p 1e-6 \
    --annotate-topsnp
```

#### Q-Q Plot
```bash
# Basic Q-Q plot
geneview qq -i gwas_results.assoc -o qq.png

# Customized Q-Q plot
geneview qq -i gwas_results.assoc -o qq.png \
    --title "GWAS QQ Plot" \
    --marker "o" --figsize 6 6
```

#### Venn Diagram
```bash
# Basic Venn diagram for 2 datasets
geneview venn -i genes_A.txt genes_B.txt -o venn2.png

# Advanced Venn diagram with 3 datasets
geneview venn -i DEG_list1.txt DEG_list2.txt DEG_list3.txt \
    --names "Study A" "Study B" "Study C" \
    --palette plasma \
    --legend-use-petal-color \
    -o venn3.png
```

#### Admixture Plot
```bash
# Basic Admixture plot
geneview admixture -i output.3.Q -p population.txt -o admixture.png

# Advanced Admixture plot
geneview admixture -i output.5.Q -p population.txt \
    --palette Set1 --edgewidth 2.0 \
    --group-order POP1 POP2 POP3 POP4 POP5 \
    --set-xticklabel-top \
    -o admixture_K5.png
```

### CLI Subcommand Help
Use `geneview <subcommand> --help` for detailed options of each command:

```bash
# Get help for specific subcommand
geneview manhattan --help
geneview qq --help
geneview venn --help
geneview admixture --help
```

**Section sources**
- [README.md:52-149](file://README.md#L52-L149)
- [geneview/cli/main.py:62-71](file://geneview/cli/main.py#L62-L71)

## Troubleshooting Common Issues

### Issue 1: Permission Errors during Installation
**Problem**: `Permission denied` errors when installing system-wide
**Solution**: Use `--user` flag or virtual environments
```bash
pip install --user geneview
# OR
pip install geneview  # in virtual environment
```

### Issue 2: matplotlib Backend Problems
**Problem**: Display issues in headless environments
**Solution**: Configure matplotlib backend
```python
import matplotlib
matplotlib.use('Agg')  # For server environments
import geneview as gv
```

### Issue 3: Large Memory Usage with Clustering
**Problem**: Slow hierarchical clustering on large datasets
**Solution**: Install fastcluster for improved performance
```bash
pip install fastcluster
```

### Issue 4: Import Errors
**Problem**: `ImportError` for missing modules
**Solution**: Reinstall with all dependencies
```bash
pip uninstall geneview
pip install geneview
```

### Issue 5: Version Compatibility Issues
**Problem**: Conflicts with existing package versions
**Solution**: Use clean virtual environment
```bash
python -m venv fresh_env
source fresh_env/bin/activate  # Linux/macOS
# OR
fresh_env\Scripts\activate  # Windows
pip install geneview
```

### Issue 6: CLI Not Found
**Problem**: `command not found` when running `geneview`
**Solution**: Ensure pip installation completed successfully and PATH is configured
```bash
# Check if CLI is available
which geneview
pip show geneview

# Reinstall if necessary
pip install --force-reinstall geneview
```

### Issue 7: CLI Version Mismatch
**Problem**: CLI shows different version than expected
**Solution**: Update to latest version
```bash
pip install --upgrade geneview
```

**Section sources**
- [.github/workflows/main.yml:11-14](file://.github/workflows/main.yml#L11-L14)
- [setup.py:17](file://setup.py#L17)
- [geneview/cli/main.py:50-53](file://geneview/cli/main.py#L50-L53)

## Performance Considerations

### Memory Management
- For large genomic datasets, consider installing `fastcluster` for improved clustering performance
- Use appropriate data types (float32 vs float64) when working with large arrays

### Caching Strategies
- GeneView caches datasets internally to reduce loading times
- Clear cache if you encounter memory issues: `gv.clear_cache()`

### Parallel Processing
- Some operations support multiprocessing for improved performance
- Configure appropriate thread counts based on your system resources

### CLI Performance Optimization
- Use appropriate file formats (PLINK2.x association files for Manhattan plots)
- Consider using compressed input files for faster processing
- Batch multiple CLI operations when possible

## Conclusion
GeneView provides a robust foundation for genomics data visualization in Python, now enhanced with a powerful command-line interface. The installation process is straightforward, with automatic dependency resolution handled by pip. Version 0.3.0 introduces seamless CLI integration alongside the traditional Python API, offering users flexibility in how they create publication-quality genomics visualizations. 

**Updated** The current version supports Python 3.9, 3.10, and 3.12, with Python 3.7 and 3.8 no longer supported. For optimal results, use virtual environments with supported Python versions, ensure compatible Python versions, and consider optional dependencies for enhanced functionality. The comprehensive CLI provides quick access to Manhattan plots, Q-Q plots, Venn diagrams, and Admixture plots without requiring any Python programming knowledge.