
import pandas as pd
from snakemake.io import apply_wildcards

configfile: "config.yaml"


samples = []
units = []
fqs = []
fast5s = []
print("Detecting available FASTQ files...")
for pattern in config["fastq-patterns"]:
    w = glob_wildcards(pattern)
    samples.extend(w.sample)
    units.extend(w.unit)
    items = list(zip(w.sample, w.unit))
    fqs.extend(apply_wildcards(pattern, {"sample": sample, "unit": unit}) for sample, unit in items)
    fast5s.extend(config["fast5-pattern"].format(sample=sample, unit=unit) for sample, unit in items)


samples = pd.DataFrame({"sample": samples, "unit": units, "fq": fqs, "fast5_pattern": fast5s}).set_index(["sample", "unit"], drop=False)


report: "report/workflow.rst"


rule all:
    input:
        expand("plots/{item.sample}-{item.unit}.{plot}.svg", item=samples.itertuples(), plot=["read-lengths", "quals"]),
        expand("plots/{sample}.signals.svg", sample=samples["sample"].unique()),
        expand("kraken/{item.sample}/{item.unit}.tsv", item=samples.itertuples())


include: "rules/qc.smk"
include: "rules/kraken.smk"
