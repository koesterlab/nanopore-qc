import matplotlib
matplotlib.use("agg")
import re
import networkx as nx
import pandas as pd
import yaml
import numpy as np
from networkx.drawing.nx_agraph import graphviz_layout, write_dot
import math


def load_classification_tree(f, min_percentage=0.0, node_fmt="{name} ({perc}%)"):
    classification = pd.read_table(f, header=None)

    classification.columns = ["percentage", "n_fragments_covered", "n_fragments_assigned", "code", "taxid", "name"]

    regexp = re.compile("(?P<indent> *)(?P<name>.+)")
    def fmt_node(row):
        m = regexp.match(row.name)
        return '{}"{}":'.format(m.group("indent"), node_fmt.format(name=m.group("name"), perc=row.percentage, taxid=row.taxid))

    tree = "\n".join(fmt_node(row) for row in classification.itertuples() if row.percentage >= min_percentage)
    tree = yaml.load(tree)


    def traverse(subtree, G, parent=None):
        for node, children in subtree.items():
            if parent is not None:
                G.add_node(node)
                G.add_edge(parent, node)
            
            if children:
                traverse(children, G, parent=node)


    G = nx.DiGraph()
    traverse(tree, G)

    return G




class TreeCMap:
    """for coloring taxids we layout the classification tree circular, and map the angle to HSL colorspace"""

    def __init__(self, tree):
        self.tree = tree
        root = next(n for n, d in tree.in_degree() if d == 0)
        self.pos = graphviz_layout(tree, prog='twopi', args='-Groot="{}"'.format(root))
        nx.set_node_attributes(tree, {node: '{},{}!'.format(*pos) for node, pos in self.pos.items()}, name="pos")
        x = [p[0] for p in self.pos.values()]
        y = [p[1] for p in self.pos.values()]
        xmin, xmax = min(x), max(x)
        ymin, ymax = min(y), max(y)

        self.center = np.array(((xmax - xmin) / 2, (ymax - ymin) / 2))
        self.maxdist = np.linalg.norm(np.array([xmax, ymax]) - self.center)

    def _normed_angle_and_mag(self, p, min_sat=0.0):
        p = np.array(p) - self.center
        ang = np.arctan2(*p[::-1])
        mag = np.linalg.norm(p) / self.maxdist
        return np.rad2deg(ang % (2 * np.pi)) / 360, min_sat + (1.0 - min_sat) * mag

    def rgb(self, taxid):
        if pd.isnull(taxid):
            return (1, 1, 1)
        try:
            angle, mag = self._normed_angle_and_mag(self.pos[taxid])
            return matplotlib.colors.hsv_to_rgb([angle, mag, 0.9])
        except KeyError:
            return (0, 0 , 0)

    def color_tree(self, tree):
        colors = {node: "#{:02x}{:02x}{:02x}".format(*map(lambda v: int(v * 255), self.rgb(node).tolist())) for node in tree}
        nx.set_node_attributes(tree, colors, name="color")

    def write_dot(self, path):
        self.color_tree(self.tree)
        write_dot(self.tree, path)
