rule merge_fastq:
    input:
        lambda w: units.loc[w.sample, "fq"]
    output:
        "reads/{sample}/all.fastq.gz"
    shell:
        "cat {input} | gzip > {output}"


rule porechop:
    input:
        "reads/{sample}/all.fastq.gz"
    output:
        expand("reads/{{sample}}/{barcode}.fastq.gz", barcode=barcodes),
        "reads/{sample}/none.fastq.gz"
    params:
        output_dir=lambda _, output: os.path.dirname(output[0])
    conda:
        "../envs/porechop.yaml"
    log:
        "logs/porechop/{sample}.log"
    threads: 8
    shell:
        "porechop -t {threads} -i {input} "
        "-b {params.output_dir} --format fastq.gz > {log} 2>&1; "
        "touch {output}"  # ensure that all output files are present
