from itertools import chain, islice
from string import ascii_uppercase
from numpy.random import choice

import matplotlib.pyplot as plt
from geneview import venn


_, top_axs = plt.subplots(ncols=3, nrows=1, figsize=(18, 5))
_, bot_axs = plt.subplots(ncols=2, nrows=1, figsize=(18, 8))

cmaps = ["cool", list("rgb"), "plasma", "viridis", "Set1"]
letters = iter(ascii_uppercase)

for n_sets, cmap, ax in zip(range(2, 7), cmaps, chain(top_axs, bot_axs)):
    dataset_dict = {
        name: set(choice(1000, 700, replace=False))
        for name in islice(letters, n_sets)
    }
    venn(dataset_dict,
         fmt="{percentage:.1f}%", # "{size}", "{logic}"
         palette=cmap,
         fontsize=12,
         legend_use_petal_color=True,
         ax=ax)

plt.tight_layout()
plt.show()

