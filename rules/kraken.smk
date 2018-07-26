rule kraken:
    input:
        reads="reads/{sample}.fq",
        db=config["kraken"]["db"]
    output:
        "kraken/{sample}.tsv"
    conda:
        "../envs/kraken.yaml"
    shell:
        "kraken --db {input.db} {input.reads} > {output}"
