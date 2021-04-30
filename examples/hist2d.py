import numpy as np
import matplotlib.pyplot as plt

import geneview as gv

x = np.random.randn(100000)
y = np.random.randn(100000) + 5
gv.hist2d(x, y)
plt.show()
