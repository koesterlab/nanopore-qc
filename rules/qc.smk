import glob


rule plot_raw:
    input:
        lambda w: glob.glob(samples.loc[w.sample, "fast5_pattern"])[config["signal"]["read"]]
    output:
        report("plots/{sample}.signals.svg", caption="../report/signal.rst", category="Raw signal")
    conda:
        "../envs/eval.yaml"
    params:
        **config["signal"]
    script:
        "../scripts/plot-signals.py"


rule merge_fastq:
    input:
        lambda w: units.loc[w.sample, "fq"]
    output:
        "reads/{sample}.fq"
    shell:
        "cat {input} > {output}"


rule fastqc:
    input:
        "reads/{sample}.fq"
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
        report("plots/{sample}.read-lengths.svg", caption="../report/read-lengths.rst", category="Read length")
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-read-lengths.py"


rule plot_quals:
    input:
        "qc/fastqc/{sample}.json"
    output:
        report("plots/{sample}.quals.svg", caption="../report/quals.rst", category="Base quality")
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-quals.py"


rule plot_kmer_mapping:
    input:
        mapping="kraken/{sample}.tsv",
        classification="tables/{sample}.classification.tsv"
    output:
        "plots/{sample}.kmer-mapping.svg"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-kmer-mapping.py"
