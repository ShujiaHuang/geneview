"""
Cache test datasets before running test suites to avoid
race conditions to due tests parallelization
"""
import geneview as gv

datasets = (
    "iris",
    "planets",
    "gwas",
    "admixture_output.Q",
    "admixture_population.info",
    "karyotype_human_hg38.txt",
    "karyotype_human_hg19.txt",
)
list(map(gv.utils.load_dataset, datasets))
