
import pandas as pd
from snakemake.io import apply_wildcards

configfile: "config.yaml"


samples = []
units = []
fqs = []
print("Detecting available FASTQ files...")
for pattern in config["fastq-patterns"]:
    w = glob_wildcards(pattern)
    samples.extend(w.sample)
    units.extend(w.unit)
    fqs.extend(apply_wildcards(pattern, {"sample": sample, "unit": unit}) for sample, unit in zip(w.sample, w.unit))

samples = pd.DataFrame({"sample": samples, "unit": units, "fq": fqs}).set_index(["sample", "unit"], drop=False)

rule all:
    input:
        expand("plots/{item.sample}-{item.unit}.{plot}.svg", item=samples.itertuples(), plot=["read-lengths", "quals"])


include: "rules/qc.smk"
