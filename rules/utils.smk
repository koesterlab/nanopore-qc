

rule calc_tree_colormap:
    input:
        "tables/{sample}-{barcode}.classification.tsv"
    output:
        "colormap/{sample}-{barcode}.pickle"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/calc-colormap.py"
