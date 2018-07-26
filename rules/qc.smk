import glob


rule plot_raw:
    input:
        lambda w: glob.glob(samples.loc[w.sample, "fast5_pattern"].iloc[0])[config["signal"]["read"]]
    output:
        report("plots/{sample}.signals.svg", caption="../report/signal.rst", category="Raw signal")
    conda:
        "../envs/eval.yaml"
    params:
        **config["signal"]
    script:
        "../scripts/plot-signals.py"


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
        report("plots/{sample}-{unit}.read-lengths.svg", caption="../report/read-lengths.rst", category="Read length")
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-read-lengths.py"


rule plot_quals:
    input:
        "qc/fastqc/{sample}/{unit}.json"
    output:
        report("plots/{sample}-{unit}.quals.svg", caption="../report/quals.rst", category="Base quality")
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-quals.py"
