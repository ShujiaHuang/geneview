geneview: A python package for genomics data visualization
==========================================================

**geneview** is a library for making attractive and informative genomics graphics in Python.
It is built on top of matplotlib and tightly integrated with the PyData stack, including
support for numpy and pandas data structures. And now it is actively developed.

Some of the features that geneview offers are:

- Functions for visualizing gerneral genomics plots.
- High-level abstractions for structuring grids of plots that let you easily build complex visualizations.


Installation
------------

To install the released version, just do

```bash
pip install geneview
```

You may instead want to use the development version from Github, by
running

```bash
pip install git+git://github.com/ShujiaHuang/geneview.git#egg=geneview
```

Tutorials
---------

-   **Manhattan** and **Q-Q** plot from the PLINK2.x GWAS result:
    [tutorial](./docs/tutorial/gwas_plot.ipynb)
-   Generate **Admixture** plot from the raw admixture output result:
    [tutorial](./docs/tutorial/admixture.ipynb)
-   **Venn diagrams for 2, 3, 4, 5, 6 sets**:
    [tutorial](./docs/tutorial/venn.ipynb)

Dependencies
------------

Geneview supports Python 3.7+ and no longer supports Python 2.

Installation requires [numpy](http://www.numpy.org/), [scipy](http://www.scipy.org/), [pandas](http://pandas.pydata.org/), and [matplotlib](http://matplotlib.org/). Some functions will use [statsmodels](http://statsmodels.sourceforge.net/).

We need the data structures: `DataFrame` and `Series` in **pandas**. 
It's easy and worth to learn, click 
[here](http://pda.readthedocs.org/en/latest/chp5.html) to see more detail 
tutorial for these two data type.

License
-------

Released under a BSD (3-clause) license
