configfile: "config.yaml"


rule all:
    input:
        "plots/read-lengths.svg"


include: "rules/qc.smk"
