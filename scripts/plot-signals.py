import json

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import h5py


def get_signals():
    for f in snakemake.input:
        f = h5py.File(f, "r")
        for read, group in f["Raw"]["Reads"].items():
            yield group["Signal"]


for signal in get_signals():
    l, u = snakemake.params.steps
    signal = signal[l:u]
    plt.plot(np.arange(signal.shape[0]), signal, "-", color="#649600")


sns.despine()
plt.xlabel("")
plt.ylabel("raw signal")
plt.xticks([])

plt.savefig(snakemake.output[0], bbox_inches="tight")

