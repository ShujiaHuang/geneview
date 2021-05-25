import matplotlib.pyplot as plt
import geneview as gv

admixture_output_fn = gv.utils.load_dataset("admixture_output.Q")
population_group_fn = gv.utils.load_dataset("admixture_population.info")

# define the order for population to plot
pop_group_1kg = ["KHV", "CDX", "CHS", "CHB", "JPT", "BEB", "STU", "ITU", "GIH", "PJL", "FIN", 
                 "CEU", "GBR", "IBS", "TSI", "PEL", "PUR", "MXL", "CLM", "ASW", "ACB", "GWD", 
                 "MSL", "YRI", "ESN", "LWK"]

f, ax = plt.subplots(1, 1, figsize=(14, 2), facecolor="w", constrained_layout=True)
gv.admixtureplot(data=admixture_output_fn, 
                 population_info=population_group_fn,
                 group_order=pop_group_1kg,
                 shuffle_popsample_kws={"frac": 0.5},
                 palette="Set1", 
                 xticklabel_kws={"rotation": 45},
                 ylabel_kws={"rotation": 0, "ha": "right"},
                 ax=ax)
plt.show()




