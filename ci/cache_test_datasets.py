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
    "admixture_population.info"
)
list(map(gv.utils.load_dataset, datasets))
