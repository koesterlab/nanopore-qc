
include: "rules/common.smk"

report: "report/workflow.rst"


targets = (expand("plots/{sample}.{plot}.svg",
                    sample=samples["sample"],
                    plot=["read-lengths", "quals", "signals"]) +
           expand("plots/{sample}-{barcode}.{plot}.svg",
                    sample=samples["sample"],
                    barcode=barcodes,
                    plot=["kmer-mapping", "classification"]))


rule all:
    input:
        targets


include: "rules/utils.smk"
include: "rules/preprocess.smk"
include: "rules/qc.smk"
include: "rules/kraken.smk"
include: "rules/eval.smk"
