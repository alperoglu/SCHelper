import scrublet as scr
import scipy.io
import numpy as np
import os
import csv
import pandas
import sys

cur_sample = sys.argv[1]
input_dir = sys.argv[2]

counts_matrix = scipy.io.mmread(input_dir + 'matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + 'features.tsv', delimiter='\t', column=1))
df = pandas.read_csv(input_dir + 'barcodes.tsv', sep="\t", header=None)


scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=85,
                                                          n_prin_comps=30)

doublets = scrub.call_doublets(threshold=0.25)

tdf = df[doublets]

tdf.to_csv('doublets_scrublet_' + cur_sample + '.csv', index=False, header = False)
