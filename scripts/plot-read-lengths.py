import json

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import seaborn as sns
import matplotlib.pyplot as plt


fastqc = json.load(open(snakemake.input[0]))

lendist = pd.DataFrame.from_records(fastqc["Sequence Length Distribution"]["contents"])
lendist["Length"] = lendist["Length"].apply(lambda l: int(l.split("-")[0]) if isinstance(l, str) else l)
lendist = lendist[lendist["Count"] > 1]
if not lendist.empty:
    sns.barplot(x="Length", y="Count", data=lendist, color="#649600")

sns.despine()
plt.yscale("log")
plt.xlabel("min read length")
plt.ylabel("count")
plt.xticks(rotation="vertical")

plt.savefig(snakemake.output[0], bbox_inches="tight")

