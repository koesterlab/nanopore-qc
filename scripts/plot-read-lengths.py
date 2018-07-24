import json

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt


fastqc = json.load(open(snakemake.input[0]))

import pdb; pdb.set_trace()
