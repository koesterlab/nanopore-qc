import pickle
import common


full_tree = common.load_classification_tree(snakemake.input[0], min_percentage=0.0)
cmap = common.TreeCMap(full_tree)
with open(snakemake.output[0], "wb") as out:
    pickle.dump(cmap, out)

