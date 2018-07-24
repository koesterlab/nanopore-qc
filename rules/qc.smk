rule fastqc:
    input:
        get_fastq
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}.zip"
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "0.27.1/bio/fastqc"


rule crimson:
    input:
        "qc/fastqc/{sample}.zip"
    output:
        "qc/fastqc/{sample}.json"
    conda:
        "../envs/crimson.yaml"
    shell:
        "crimson fastqc {input} {output}"


rule plot_read_lengths:
    input:
        "qc/fastqc/{sample}.json"
    output:
        "plots/read-lengths.svg"
    conda:
        "../envs/eval.yaml"
    script:
        "scripts/plot-read-lenghts.py"
