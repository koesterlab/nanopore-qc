
import pandas as pd
from snakemake.io import apply_wildcards

configfile: "config.yaml"


logger.info("Detecting available FASTQ files...")
samples = pd.DataFrame(columns=["sample", "fast5_pattern"])
units = pd.DataFrame(columns=["sample", "unit", "fq"])
for pattern in config["fastq-patterns"]:
    w = glob_wildcards(pattern)
    for sample in set(w.sample):
        samples = samples.append({"sample": sample, "fast5_pattern": config["fast5-pattern"].format(sample=sample)}, ignore_index=True)
    for sample, unit in zip(w.sample, w.unit):
        units = units.append({"sample": sample, "unit": unit, "fq": apply_wildcards(pattern, {"sample": sample, "unit": unit})}, ignore_index=True)
samples = samples.set_index("sample", drop=False)
units = units.set_index(["sample", "unit"], drop=False)

report: "report/workflow.rst"


targets_qc = expand("plots/{sample}.{plot}.svg", sample=samples["sample"], plot=["read-lengths", "quals", "signals", "kmer-mapping"])

targets_classification = expand("plots/{sample}.classification.svg", sample=samples["sample"])


rule all:
    input:
        targets_qc, targets_classification


rule qc:
    input:
        targets_qc


include: "rules/qc.smk"
include: "rules/kraken.smk"
