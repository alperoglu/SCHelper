import warnings

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
    """
    Evaluate a logical formula based on the values in a dictionary.

    Parameters:
    formula (str): The logical formula to evaluate. The formula can contain the following elements:
                   - Variable names that are keys in the dictionary x.
                   - The logical operators '!', '+', and 'x' (for AND).
                   - Parentheses for grouping.
    x (dict): A dictionary that maps variable names to boolean values.

    Returns:
    bool: The result of evaluating the formula. If the formula is empty, the function returns True.
          If the formula contains a variable name that is not in the dictionary, the function returns False.
    """
    import re

    # If the formula is empty, return True
    if formula == "":
        return True

    # If the formula contains a variable name that is not in the dictionary, return False
    if not re.search(r"(x|\(|\)|\+|\!)", formula) and formula not in x.keys():
        return False

    # If the formula is a variable name, return its value from the dictionary
    if not re.search(r"(x|\(|\)|\+|\!)", formula):
        return x[formula]

    # If the formula is a negated variable name and the variable is not in the dictionary, return False
    if not re.search(r"(x|\(|\)|\+)", formula) and re.search(r"!", formula) and formula.replace("!", "") not in x.keys():
        return False

    # If the formula is a negated variable name, return the negated value from the dictionary
    if not re.search(r"(x|\(|\)|\+)", formula) and re.search(r"!", formula):
        return not x[formula.replace("!", "")]

    # If the formula is a disjunction (OR operation), evaluate each part of the disjunction and return True if any part is True
    if not re.search(r"(x|\(|\))", formula) and re.search(r"\+", formula):
        return any([evaluate_formula(f, x) for f in formula.split("+")])

    # If the formula is a conjunction (AND operation) or contains parentheses, evaluate each part of the conjunction or parentheses and return True if all parts are True
    if re.search(r"(x|\(|\))", formula):
        return all([evaluate_formula(f, x) for f in re.split(r"x|\(|\)", formula)])
    
def entropy(p):
    """
    Calculate the entropy of a list of percentages.

    Parameters:
    p (list): A list of percentages. Each percentage is a number between 0 and 1.

    Returns:
    float: The entropy of the list of percentages. The entropy is a measure of the uncertainty or randomness of the data.
           The entropy is 0 if all percentages are 0.
    """
    import numpy as np

    # Initialize the entropy to 0
    entropy = 0

    # Iterate over the list of percentages
    for i in range(len(p)):
        # If the percentage is 0, add 0 to the entropy
        if p[i] == 0:
           entropy += 0
        # Otherwise, add the percentage times its logarithm to the entropy
        else:
            entropy += p[i] * np.log(p[i])

    # Return the negative of the entropy
    return -entropy

def mutual_info_score(contingency=None):
    """
    Calculate the mutual information score of a contingency table.

    Parameters:
    contingency (numpy.ndarray): A 2D array that represents a contingency table.

    Returns:
    float: The mutual information score. The mutual information score is a measure of the mutual dependence between the two variables.
           It is equal to zero if and only if two random variables are independent, and higher values mean higher dependency.
    """
    import numpy as np

    # Get the indices of the non-zero elements in the contingency table
    nzx, nzy = np.nonzero(contingency)
    # Get the non-zero elements in the contingency table
    nz_val = contingency[nzx, nzy]

    # Calculate the sum of the contingency table
    contingency_sum = contingency.sum()
    # Calculate the sum of each row in the contingency table
    pi = np.ravel(contingency.sum(axis=1))
    # Calculate the sum of each column in the contingency table
    pj = np.ravel(contingency.sum(axis=0))

    # Calculate the logarithm of the non-zero elements in the contingency table
    log_contingency_nm = np.log(nz_val)
    # Normalize the non-zero elements in the contingency table
    contingency_nm = nz_val / contingency_sum

    # Calculate the outer product of the row sums and the column sums
    outer = pi.take(nzx) * pj.take(nzy)
    # Calculate the logarithm of the outer product
    log_outer = -np.log(outer) + np.log(pi.sum()) + np.log(pj.sum())

    # Calculate the mutual information score
    mi = (
        contingency_nm * (log_contingency_nm - np.log(contingency_sum))
        + contingency_nm * log_outer
    )

    # Replace very small values with zero
    mi = np.where(np.abs(mi) < np.finfo(mi.dtype).eps, 0.0, mi)

    # Return the mutual information score, clipped to be at least zero
    return np.clip(mi.sum(), 0.0, None)

def norm_mi(contingency=None, average_method='arithmetic'):
    """
    Calculate the normalized mutual information of a contingency table.

    Parameters:
    contingency (numpy.ndarray): A 2D array that represents a contingency table.
    average_method (str): The method to use to average the entropies of the two variables. 
                          Options are 'min', 'geometric', 'arithmetic', 'max', 'custom', and 'default'.

    Returns:
    float: The normalized mutual information. This is the mutual information score divided by a normalizer based on the average_method.
    """
    import numpy as np

    # Calculate the mutual information score of the contingency table
    mi = mutual_info_score(contingency)

    # Calculate the entropies of the two variables
    h_1, h_2 = entropy(contingency.sum(axis=1)), entropy(contingency.sum(axis=0))

    # Calculate the normalizer based on the average_method
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

    # Return the normalized mutual information
    return mi / normalizer


def calculate_nmi(adata, clustering_key, control_key):
    """
    Calculate the normalized mutual information between two categorical variables from adata.obs.

    Parameters:
    adata (AnnData): The annotated data matrix.
    clustering_key (str): The key of the first categorical variable in adata.obs.
    control_key (str): The key of the second categorical variable in adata.obs.

    Returns:
    float: The normalized mutual information between the two categorical variables.
    """
    import pandas as pd

    # Convert the clustering key to string and control key to category
    adata.obs[clustering_key] = adata.obs[clustering_key].astype(str)
    adata.obs[control_key] = adata.obs[control_key].astype("category")

    # Create the confusion matrix
    confusion_matrix = pd.crosstab(adata.obs[control_key], adata.obs[clustering_key])

    # Normalize the confusion matrix by column so that they sum up to 1
    confusion_matrix = confusion_matrix.div(confusion_matrix.values.sum())

    # Calculate the normalized mutual information
    nmi = norm_mi(confusion_matrix.values, average_method="custom")

    return nmi

def odds_ratio_differential(refCounts, otherCounts):
    """
    Calculate the log odds ratio differential between two sets of counts.

    Parameters:
    refCounts (numpy.ndarray): The counts for the reference set.
    otherCounts (numpy.ndarray): The counts for the other set.

    Returns:
    float: The log odds ratio differential. This is a measure of the difference in the odds ratios between the two sets of counts.
    """
    import numpy as np
    from scipy.special import betaln

    # Calculate the total counts for the reference set and the other set
    totalRefCounts = np.sum(refCounts)
    totalOtherCounts = np.sum(otherCounts)

    # Calculate the terms of the log odds ratio differential
    term1 = betaln(refCounts + 1, totalRefCounts - refCounts + 1)
    term2 = betaln(otherCounts + 1, totalOtherCounts - otherCounts + 1)
    term3 = betaln(refCounts + otherCounts + 1, totalRefCounts + totalOtherCounts - refCounts - otherCounts + 1)

    # Calculate the log odds ratio differential
    result = (term1 + term2) - term3
    result = result / np.log(10)

    return result


def odds_ratio_differential_efficient(refCounts, totalRefCounts, otherCounts, totalOtherCounts):
    """
    Calculate the log odds ratio differential between two sets of counts, given the total counts.

    Parameters:
    refCounts (numpy.ndarray): The counts for the reference set.
    totalRefCounts (int): The total counts for the reference set.
    otherCounts (numpy.ndarray): The counts for the other set.
    totalOtherCounts (int): The total counts for the other set.

    Returns:
    float: The log odds ratio differential. This is a measure of the difference in the odds ratios between the two sets of counts.
    """
    import numpy as np
    from scipy.special import betaln

    # Calculate the terms of the log odds ratio differential
    term1 = betaln(refCounts + 1, totalRefCounts - refCounts + 1)
    term2 = betaln(otherCounts + 1, totalOtherCounts - otherCounts + 1)
    term3 = betaln(refCounts + otherCounts + 1, totalRefCounts + totalOtherCounts - refCounts - otherCounts + 1)

    # Calculate the log odds ratio differential
    result = (term1 + term2) - term3
    result = result / np.log(10)

    return result