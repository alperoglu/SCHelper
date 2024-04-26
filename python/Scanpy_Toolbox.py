import os
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import pairwise_distances
from scipy.sparse import coo_matrix
import leidenalg
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.special import betaln
import re
import warnings
import json

warnings.filterwarnings("ignore")

def plot_gene_coexpression(adata, genes, layer, dim_space):
    import scipy
    import scanpy as sc
    import numpy as np
    import pandas as pd
    
    """
    Plot the co-expression of a list of genes in a given dimension space.

    Parameters:
    adata (AnnData): The annotated data matrix.
    genes (list): The list of genes to plot.
    layer (str): The layer of the adata to use.
    dim_space (str): The dimension space to plot (e.g., 'pca', 'umap').

    Returns:
    None
    """
    # Get the indices of the genes
    gene_indices = [adata.var_names.get_loc(gene) for gene in genes]

    # check if the layer exists
    if layer not in adata.layers.keys():
        print(f"Layer {layer} not found in adata.")
        return

    if type(adata.layers[layer]) == np.ndarray:
        coex = np.all([adata.layers[layer][:, gene_index] > 0 for gene_index in gene_indices], axis=0)
    elif type(adata.layers[layer]) == scipy.sparse.csr.csr_matrix:
        coex = np.all([adata.layers[layer][:, gene_index].todense() > 0 for gene_index in gene_indices], axis=0)
    else:
        print(f"Unknown layer type: {type(adata.layers[layer])}")
        return
    
    coex_list = [item for sublist in coex.tolist() for item in sublist]
    coex_text = 'CoEx' + ''.join([f'_{gene}' for gene in genes])
    adata.obs[coex_text] = pd.Categorical(coex_list, categories=[True, False])

    # Plot the co-expression in the specified dimension space
    if dim_space == 'pca':
        sc.pl.pca(adata, color=coex_text, groups=[True])
    elif dim_space == 'umap':
        sc.pl.umap(adata, color=coex_text, groups=[True])
    else:
        print(f"Unknown dimension space: {dim_space}")

def evaluate_formula(formula, x):
    if formula == "":
        return True
    if not re.search(r"(x|\(|\)|\+|\!)", formula) and formula not in x.keys():
        return False
    if not re.search(r"(x|\(|\)|\+|\!)", formula):
        return x[formula]
    if not re.search(r"(x|\(|\)|\+)", formula) and re.search(r"!", formula) and formula.replace("!", "") not in x.keys():
        return False
    if not re.search(r"(x|\(|\)|\+)", formula) and re.search(r"!", formula):
        return not x[formula.replace("!", "")]
    if not re.search(r"(x|\(|\))", formula) and re.search(r"\+", formula):
        return any([evaluate_formula(f, x) for f in formula.split("+")])
    if re.search(r"(x|\(|\))", formula):
        return all([evaluate_formula(f, x) for f in re.split(r"x|\(|\)", formula)])
    
# entropy function for a array of percentages and return zero if percentage is zero
def entropy(p):
    tmp = 0
    for i in range(len(p)):
        if p[i] == 0:
           tmp = tmp + 0
        else:
            tmp = tmp + p[i] * np.log(p[i])
    return -tmp

def mutual_info_score(contingency=None):
    nzx, nzy = np.nonzero(contingency)
    nz_val = contingency[nzx, nzy]

    contingency_sum = contingency.sum()
    pi = np.ravel(contingency.sum(axis=1))
    pj = np.ravel(contingency.sum(axis=0))

    log_contingency_nm = np.log(nz_val)
    contingency_nm = nz_val / contingency_sum

    outer = pi.take(nzx) * pj.take(nzy)
    log_outer = -np.log(outer) + np.log(pi.sum()) + np.log(pj.sum())

    mi = (
        contingency_nm * (log_contingency_nm - np.log(contingency_sum))
        + contingency_nm * log_outer
    )

    mi = np.where(np.abs(mi) < np.finfo(mi.dtype).eps, 0.0, mi)

    return np.clip(mi.sum(), 0.0, None)

def norm_mi(contingency = None, average_method = 'arithmetic'):
    
    mi = mutual_info_score(contingency)
    
    h_1, h_2 = entropy(contingency.sum(axis=1)), entropy(contingency.sum(axis=0))

    if average_method == "min":
        normalizer = min(h_1, h_2)
    elif average_method == "geometric":
        normalizer = np.sqrt(h_1 * h_2)
    elif average_method == "arithmetic":
        normalizer = np.mean([h_1, h_2])
    elif average_method == "max":
        normalizer = max(h_1, h_2)
    elif average_method == "custom":
        normalizer = h_1
    else:
        normalizer = 1

    return mi / normalizer

# function to calculate the normalized mutual information between two categorical variables from adata.obs using norm_mi function
def calculate_nmi(adata, clustering_key, control_key):
        
        # convert the clustering key to string and control key to category 
        adata.obs[clustering_key] = adata.obs[clustering_key].astype(str)
        adata.obs[control_key] = adata.obs[control_key].astype("category")
    
        # create the confusion matrix 
        confusion_matrix = pd.crosstab(adata.obs[control_key], adata.obs[clustering_key])

        # normalize the confusion matrix by column so thay they sum up to 1
        confusion_matrix = confusion_matrix.div(confusion_matrix.values.sum())

        nmi = norm_mi(confusion_matrix.values, average_method="custom")
    
        return nmi

def odds_ratio_differential(refCounts, otherCounts):
    
    totalRefCounts = np.sum(refCounts)
    totalOtherCounts = np.sum(otherCounts)
    
    term1 = betaln(refCounts + 1, totalRefCounts - refCounts + 1)
    term2 = betaln(otherCounts + 1, totalOtherCounts - otherCounts + 1)
    term3 = betaln(refCounts + otherCounts + 1, totalRefCounts + totalOtherCounts - refCounts - otherCounts + 1)
    
    result = (term1 + term2) - term3
    result = result / np.log(10)
    
    return result

def odds_ratio_differential_efficient(refCounts, totalRefCounts, otherCounts, totalOtherCounts):
        
    term1 = betaln(refCounts + 1, totalRefCounts - refCounts + 1)
    term2 = betaln(otherCounts + 1, totalOtherCounts - otherCounts + 1)
    term3 = betaln(refCounts + otherCounts + 1, totalRefCounts + totalOtherCounts - refCounts - otherCounts + 1)
    
    result = (term1 + term2) - term3
    result = result / np.log(10)
    
    return result