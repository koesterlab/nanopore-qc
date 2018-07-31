
rule calc_tree_colormap:
    input:
        "tables/{sample}.classification.tsv"
    output:
        "colormap/{sample}.pickle"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/calc-colormap.py"
