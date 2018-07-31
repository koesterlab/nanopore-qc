import pickle
from networkx.drawing.nx_agraph import write_dot
import common

cmap = pickle.load(open(snakemake.input.colormap, "rb"))

tree = common.load_classification_tree(snakemake.input.classification, min_percentage=1.0)
cmap.color_tree(tree)

write_dot(tree, snakemake.output[0])

