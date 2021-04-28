
Geneview: A python package for genomics data visualization
==========================================================

Geneview is a Python visualization library base on matplotlib. It provides a 
high-level interface for drawing attractive genomics graphics. And now it is 
actively developed.

Documentation
-------------

Examples
--------

..

   We use a dataset ``GOYA_preview`` which is in 
   `geneview-data <https://github.com/ShujiaHuang/geneview-data>`_ repository, as 
   the input for the plots below. Here is the format preview of ``GOYA_preview``\ :

   .. image:: https://github.com/ShujiaHuang/geneview/blob/master/examples/data/goya_format.png
      :target: https://github.com/ShujiaHuang/geneview/blob/master/examples/data/goya_format.png
      :alt: GOYA format


A basic example for **Manhattan** plot.

.. code-block:: python

   import matplotlib.pyplot as plt
   import geneview as gv

   df = gv.util.load_dataset('GOYA_preview')  # df is DataFrame of pandas
   xtick = ['chr'+c for c in map(str, range(1, 15) + ['16', '18', '20', '22'])]
   gv.gwas.manhattanplot(df[['chrID','position','pvalue']],  
                            xlabel="Chromosome", 
                            ylabel="-Log10(P-value)", 
                            xticklabel_kws={'rotation': 'vertical'},
                            xtick_label_set = set(xtick))
   plt.show()


.. image:: https://github.com/ShujiaHuang/geneview/blob/master/examples/manhattan.png
   :target: https://github.com/ShujiaHuang/geneview/blob/master/examples/manhattan.png
   :alt: manhattanplot


A basic example for **QQ** plot.

.. code-block:: python

   import matplotlib.pyplot as plt
   import geneview as gv

   df = gv.util.load_dataset('GOYA_preview')  # df is DataFrame of pandas
   gv.gwas.qqplot(df['pvalue'], color="#00bb33",
                  xlabel="Expected p-value(-log10)",
                  ylabel="Observed p-value(-log10)")
   plt.show()


.. image:: https://github.com/ShujiaHuang/geneview/blob/master/examples/qq.png
   :target: https://github.com/ShujiaHuang/geneview/blob/master/examples/qq.png
   :alt: qqplot


Citing
------

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

Installation
------------

To install the released version, just do

.. code-block::

   pip install geneview

You may instead want to use the development version from Github, by running

.. code-block::

   pip install git+git://github.com/ShujiaHuang/geneview.git#egg=geneview

Testing
-------

Development
-----------

License
-------

Released under a BSD (3-clause) license
