import matplotlib
matplotlib.use("Agg")
from networkx.drawing.nx_agraph import write_dot
import common


full_tree = common.load_classification_tree(snakemake.input[0], min_percentage=0.0)
cmap = common.TreeCMap(full_tree)
colors = {node: tuple(cmap.rgb(node)) for node in full_tree}

tree = common.load_classification_tree(snakemake.input[0], min_percentage=1.0)
nx.set_node_attributes(tree, colors, name="color")

write_dot(tree, snakemake.output[0])

