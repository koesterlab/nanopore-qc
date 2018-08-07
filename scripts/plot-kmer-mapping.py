import pickle
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import common


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
    mappings.append(m)
    if i == 100:
        break

if mappings:
    mappings = pd.concat(mappings)
    mappings = mappings.taxid.unstack()
    mappings = mappings.iloc[:, :100]

    # color image in HSV representation of the circular layout of the taxid tree
    cmap = pickle.load(open(snakemake.input.colormap, "rb"))
    img = mappings.applymap(cmap.rgb)

    shape = list(img.shape) + [3]
    img = np.vstack(img.values.tolist()).reshape(shape)

    #plt.figure(figsize=(30, 10))

    plt.imshow(img, interpolation="nearest", aspect="auto")

    plt.xlabel("k-mer")
    plt.ylabel("read")

plt.savefig(snakemake.output[0])
