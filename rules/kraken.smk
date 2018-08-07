rule kraken:
    input:
        reads="reads/{sample}/{barcode}.fastq.gz",
        db=config["kraken"]["db"]
    output:
        "kraken/{sample}-{barcode}.tsv"
    conda:
        "../envs/kraken.yaml"
    log:
        "logs/kraken/{sample}-{barcode}.log"
    threads: 64
    shell:
        """
        if [ -s {input.reads} ]
        then
            kraken --threads {threads} --db {input.db} {input.reads} > {output} 2> {log}
        else
            touch {output}
        fi
        """


rule kraken_report:
    input:
        tsv="kraken/{sample}-{barcode}.tsv",
        db=config["kraken"]["db"]
    output:
        report("tables/{sample}-{barcode}.classification.tsv", caption="../report/kraken.rst", category="classification")
    conda:
        "../envs/kraken.yaml"
    log:
        "logs/kraken-report/{sample}-{barcode}.log"
    shell:
        "kraken-report --db {input.db} {input.tsv} > {output} 2> {log}"


rule extract_classification_tree:
    input:
        classification="tables/{sample}-{barcode}.classification.tsv",
        colormap="colormap/{sample}-{barcode}.pickle"
    output:
        "kraken/{sample}-{barcode}.classification.dot"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-classification.py"
