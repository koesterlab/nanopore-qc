
include: "rules/common.smk"

report: "report/workflow.rst"


targets_qc = expand("plots/{sample}.{plot}.svg",
                    sample=samples["sample"],
                    plot=["read-lengths", "quals", "signals", "kmer-mapping"])
targets_classification = expand("plots/{sample}.classification.svg",
                                sample=samples["sample"])


rule all:
    input:
        targets_qc, targets_classification


rule qc:
    input:
        targets_qc


include: "rules/utils.smk"
include: "rules/qc.smk"
include: "rules/kraken.smk"
include: "rules/eval.smk"
