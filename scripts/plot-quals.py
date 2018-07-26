import json

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

fastqc = json.load(open(snakemake.input[0]))

qualdist = pd.DataFrame.from_records(fastqc["Per base sequence quality"]["contents"])
qualdist.Base = qualdist.Base.apply(lambda b: int(b.split("-")[0]))
qualdist = qualdist.replace("NaN", np.NaN)
single_value = qualdist["Lower Quartile"].isnull()
qualdist = qualdist[~single_value]

if not qualdist.empty:
    plt.fill_between(x=qualdist.Base, y2=qualdist["Lower Quartile"], y1=qualdist["Upper Quartile"], color="#cee3a3")
    plt.plot(qualdist.Base, qualdist.Mean, "-", color="#649600")

sns.despine()
plt.xlabel("read base")
plt.ylabel("quality (PHRED scale)")
plt.xticks(rotation="vertical")

plt.savefig(snakemake.output[0], bbox_inches="tight")

