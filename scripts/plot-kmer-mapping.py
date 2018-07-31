import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import common
from networkx.drawing.nx_agraph import graphviz_layout


tree = common.load_classification_tree(snakemake.input.classification, min_percentage=0.0, node_fmt="{taxid}")


d = pd.read_table(snakemake.input.mapping, names=["status", "seqid", "taxid", "seqlen", "kmer_mapping"])

mappings = []
for i, per_read_mapping in enumerate(d.kmer_mapping.str.split(" ")):
    j = 0
    def mapping_stretches(row):
        global j
        taxid, count = row.split(":")
        count = int(count)
        yield from ([taxid, x] for x in range(j, j + count))
        j += count

    m = pd.DataFrame([item for row in per_read_mapping for item in mapping_stretches(row)])
    m.columns = ["taxid", "kmer"]
    m["read"] = i
    m = m.set_index(["read", "kmer"], drop=False)
    major_taxid = m.taxid.value_counts().index[0]
    def encode(taxid):
        if taxid == major_taxid:
            return 2
        elif taxid == 0:
            return 1
        else:
            return 0
    m["is_major"] = m.taxid.apply(encode)
    mappings.append(m)
    if i == 100:
        break

mappings = pd.concat(mappings)

#import pdb; pdb.set_trace()
#sns.scatterplot(x="kmer", y="read", hue="taxid", data=mappings, size=2, edgecolors="face", marker=".", rasterized=True, legend=False, palette="Set2")
#sns.stripplot(x="kmer", y="read", hue="taxid", data=mappings, palette="tab20", rasterized=True, jitter=0.0)
mappings = mappings.taxid.unstack()
mappings = mappings.iloc[:, :100]

# color image in HSV representation of the circular layout of the taxid tree
cmap = common.TreeCMap(tree)
img = mappings.applymap(cmap.rgb)
shape = list(img.shape) + [3]
img = np.vstack(img.values.tolist()).reshape(shape)

plt.imshow(img, interpolation="nearest", aspect="auto")

plt.xlabel("k-mer")
plt.ylabel("read")

plt.savefig(snakemake.output.svg)

cmap.write_dot(snakemake.output.dot)
