rule kraken:
    input:
        reads="reads/{sample}.fq",
        db=config["kraken"]["db"]
    output:
        "kraken/{sample}.tsv"
    conda:
        "../envs/kraken.yaml"
    threads: 64
    shell:
        "kraken --threads {threads} --db {input.db} {input.reads} > {output}"


rule kraken_report:
    input:
        tsv="kraken/{sample}.tsv",
        db=config["kraken"]["db"]
    output:
        report("tables/{sample}.classification.tsv", caption="../report/kraken.rst", category="classification")
    conda:
        "../envs/kraken.yaml"
    shell:
        "kraken-report --db {input.db} {input.tsv} > {output}"


rule extract_classification_tree:
    input:
        "tables/{sample}.classification.tsv"
    output:
        "kraken/{sample}.classification.dot"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-classification.py"


rule plot_classification_tree:
    input:
        "kraken/{sample}.classification.dot"
    output:
        report("plots/{sample}.classification.svg", caption="../report/classification-tree.rst", category="classification")
    conda:
        "../envs/eval.yaml"
    params:
        color=config["style"]["color"]
    shell:
        "dot -Tsvg {input} "
        "-Grankdir=TB -Nshape=box -Nstyle=rounded -Nfontname=sans "
        "-Nfontsize=10 -Npenwidth=2 -Epenwidth=2 -Ecolor=grey -Nbgcolor=white -Ncolor='{params.color}'"
        "> {output}"
