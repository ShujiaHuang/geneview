import matplotlib.pyplot as plt
import geneview as gv

df = gv.utils.load_dataset("gwas")
gv.qqplot(df.loc[:, "P"]) 

plt.tight_layout()
plt.savefig("images/qq.png")
plt.show()
