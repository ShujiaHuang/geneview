Geneview: A python package for genomics data visualization
==========================================================

Geneview is a Python visualization library base on matplotlib. It provides a 
high-level interface for drawing attractive genomics graphics. And now it is 
actively developed.

Documentation
-------------

Examples
--------

A basical example for Manhattan plot.

```
import matplotlib.pyplot as plt
import geneview as gv

df = gv.util.load_dataset('GOYA_preview')  # df is DataFrame of pandas
xtick = map(str, range(1, 15) + ['16','18', '20','22'])
gv.gwas.manhattanplot(df[['chrID','position','pvalue']],  
                         xlabel="Chromosome", 
                         ylabel="-Log10(P-value)", 
                         xtick_label_set = set(xtick))
plt.show()
```

![manhattanplot](https://github.com/ShujiaHuang/geneview/blob/master/examples/manhattan.png)

A basical example for QQ plot.

```
import matplotlib.pyplot as plt
import geneview as gv

df = gv.util.load_dataset('GOYA_preview')  # df is DataFrame of pandas
gv.gwas.qqplot(df['pvalue'], color="#00bb33",
               xlabel="Expected p-value(-log10)",
               ylabel="Observed p-value(-log10)")
plt.show()
```
![qqplot](https://github.com/ShujiaHuang/geneview/blob/master/examples/qq.png)

> We use a dataset `GOYA_preview` which is in 
[geneview-data](https://github.com/ShujiaHuang/geneview-data) repository, as 
the input for the plots above. Here is the format preview of `GOYA_preview`:
![GOYA format](https://github.com/ShujiaHuang/geneview/blob/master/examples/data/goya_format.png)


Citing
------

Dependencies
------------

- Python 2.7 or Python 3.4+

### Mandatory

- [numpy](http://www.numpy.org/)

- [scipy](http://www.scipy.org/)

- [pandas](http://pandas.pydata.org/)

- [matplotlib](http://matplotlib.org/)

We need the data structures: `DataFrame` and `Seris` in **pandas**. It's easy 
and worth to learn, click [here](http://pda.readthedocs.org/en/latest/chp5.html) 
to see more detail tutorial for these two data type.

Installation
------------

To install the released version, just do
```
pip install seaborn
```

You may instead want to use the development version from Github, by running

```
pip install git+git://github.com/ShujiaHuang/geneview.git#egg=geneview
```

Testing
-------

Development
-----------

License
-------

Released under a BSD (3-clause) license
