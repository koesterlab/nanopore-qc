rule fastqc:
    input:
        lambda w: samples.loc[(w.sample, w.unit), "fq"]
    output:
        html="qc/fastqc/{sample}/{unit}.html",
        zip="qc/fastqc/{sample}/{unit}.zip"
    params: ""
    log:
        "logs/fastqc/{sample}/{unit}.log"
    wrapper:
        "0.27.1/bio/fastqc"


rule crimson:
    input:
        "qc/fastqc/{sample}/{unit}.zip"
    output:
        "qc/fastqc/{sample}/{unit}.json"
    conda:
        "../envs/crimson.yaml"
    shell:
        "crimson fastqc {input} {output}"


rule plot_read_lengths:
    input:
        "qc/fastqc/{sample}/{unit}.json"
    output:
        report("plots/{sample}-{unit}.read-lengths.svg", caption="report/read-lengths.rst")
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-read-lengths.py"


rule plot_quals:
    input:
        "qc/fastqc/{sample}/{unit}.json"
    output:
        report("plots/{sample}-{unit}.quals.svg", caption="report/quals.rst")
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-quals.py"
