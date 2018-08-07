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
barcodes = expand("BC{barcode:02d}", barcode=range(1, 13))
