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
        mapping="kraken/{sample}-{barcode}.tsv",
        colormap="colormap/{sample}-{barcode}.pickle"
    output:
        report("plots/{sample}-{barcode}.kmer-mapping.svg", caption="../report/kmer-mapping.rst", category="K-mer mapping")
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-kmer-mapping.py"


rule plot_classification_tree:
    input:
        "kraken/{sample}-{barcode}.classification.dot"
    output:
        report("plots/{sample}-{barcode}.classification.svg", caption="../report/classification-tree.rst", category="classification")
    conda:
        "../envs/eval.yaml"
    params:
        color=config["style"]["color"]
    shell:
        "dot -Tsvg {input} "
        "-Grankdir=TB -Nshape=box -Nstyle=rounded -Nfontname=sans "
        "-Nfontsize=10 -Npenwidth=2 -Epenwidth=2 -Ecolor=grey -Nbgcolor=white " # -Ncolor='{params.color}'"
        "> {output}"
