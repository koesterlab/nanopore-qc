rule kraken:
    input:
        lambda w: samples.loc[(w.sample, w.unit), "fq"]
    output:
        "kraken/{sample}/{unit}.tsv"
    conda:
        "../envs/kraken.yaml"
    params:
        **config["kraken"]
    shell:
        "kraken --db {params.db} {input} > {output}"
