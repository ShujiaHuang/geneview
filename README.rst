geneview: A python package for genomics data visualization
==========================================================

**geneview** is a Python visualization library base on matplotlib. It provides a
high-level interface for drawing attractive genomics graphics. And now it is 
actively developed.

Installation
------------

To install the released version, just do

.. code-block::

   pip install geneview

You may instead want to use the development version from Github, by running

.. code-block::

   pip install git+git://github.com/ShujiaHuang/geneview.git#egg=geneview

Tutorials
---------

* **Manhattan** and **Q-Q** plot: `tutorial <./docs/tutorial/gwas_plot.ipynb>`_
* **Venn diagrams for 2, 3, 4, 5, 6 sets**: `tutorial <./docs/tutorial/venn.ipynb>`_

Dependencies
------------

* Python 3.7+

Mandatory
^^^^^^^^^

* `numpy <http://www.numpy.org/>`_
* `scipy <http://www.scipy.org/>`_
* `pandas <http://pandas.pydata.org/>`_
* `matplotlib <http://matplotlib.org/>`_

Recommended
^^^^^^^^^^^

* `statsmodels <http://statsmodels.sourceforge.net/>`_

We need the data structures: ``DataFrame`` and ``Series`` in **pandas**. It's easy 
and worth to learn, click `here <http://pda.readthedocs.org/en/latest/chp5.html>`_ 
to see more detail tutorial for these two data type.

License
-------

Released under a BSD (3-clause) license


