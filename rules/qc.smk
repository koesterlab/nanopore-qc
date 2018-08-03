import glob


rule fastqc:
    input:
        "reads/{sample}/{barcode}.fastq"
    output:
        html="qc/fastqc/{sample}-{barcode}.html",
        zip="qc/fastqc/{sample}-{barcode}.zip"
    params: ""
    log:
        "logs/fastqc/{sample}-{barcode}.log"
    wrapper:
        "0.27.1/bio/fastqc"


rule crimson:
    input:
        "qc/fastqc/{sample}-{barcode}.zip"
    output:
        "qc/fastqc/{sample}-{barcode}.json"
    conda:
        "../envs/crimson.yaml"
    shell:
        "crimson fastqc {input} {output}"
