import matplotlib
matplotlib.use("Agg")
import networkx as nx
from networkx.drawing.nx_agraph import write_dot
import common


full_tree = common.load_classification_tree(snakemake.input[0], min_percentage=0.0)
cmap = common.TreeCMap(full_tree)

tree = common.load_classification_tree(snakemake.input[0], min_percentage=1.0)
cmap.color_tree(tree)

write_dot(tree, snakemake.output[0])

