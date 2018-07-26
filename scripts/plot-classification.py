import re
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import yaml
from networkx.drawing.nx_agraph import write_dot


classification = pd.read_table(snakemake.input[0], header=None)

classification.columns = ["percentage", "n_fragments_covered", "n_fragments_assigned", "code", "taxid", "name"]
print(classification)

regexp = re.compile("(?P<indent> *)(?P<name>.+)")
def fmt_node(row):
    m = regexp.match(row.name)
    return '{}"{} ({}%)":'.format(m.group("indent"), m.group("name"), row.percentage)

tree = "\n".join(fmt_node(row) for row in classification.itertuples() if row.percentage >= 1.0)
tree = yaml.load(tree)


def traverse(subtree, G, parent=None):
    for node, children in subtree.items():
        if parent is not None:
            print(node)
            G.add_node(node)
            G.add_edge(parent, node)
            
        if children:
            traverse(children, G, parent=node)


G = nx.DiGraph()
traverse(tree, G)


write_dot(G, snakemake.output[0])

