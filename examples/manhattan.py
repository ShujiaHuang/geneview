import matplotlib.pyplot as plt
import geneview as gv

xtick = ["chr"+c for c in list(map(str, range(1, 15))) + ["16", "18", "20", "22"]]
df = gv.utils.load_dataset("gwas")
gv.manhattanplot(data=df.loc[:,["#CHROM","POS","P"]], 
                 hline_kws={"linestyle": "--", "lw": 1.3},
                 xlabel="Chromosome", 
                 ylabel="-Log10(P-value)",
                 xticklabel_kws={"rotation": 45},
                 xtick_label_set = set(xtick))
plt.tight_layout()
plt.savefig("images/manhattan.png")
plt.show()
