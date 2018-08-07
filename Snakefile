
include: "rules/common.smk"

report: "report/workflow.rst"


def targets(samples):
    return (expand("plots/{sample}.{plot}.svg",
                    sample=samples,
                    plot=["read-lengths", "quals", "signals"]) +
           expand("plots/{sample}-{barcode}.{plot}.svg",
                    sample=samples,
                    barcode=barcodes,
                    plot=["kmer-mapping", "classification"]))


rule all:
    input:
        targets(samples["sample"])


include: "rules/utils.smk"
include: "rules/preprocess.smk"
include: "rules/qc.smk"
include: "rules/kraken.smk"
include: "rules/eval.smk"


rule per_sample_report:
    input:
        lambda w: targets([w.sample])
    output:
        "reports/{sample}.report.html"
    shell:
        "snakemake {input} --report {output}"
