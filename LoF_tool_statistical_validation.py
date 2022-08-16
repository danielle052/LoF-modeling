import numpy as np
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt


# p-value correction for all phenotypes
def bh_correction_for_pvalue(p_value_dict):
    p_values = list(p_value_dict.values())
    p_values.sort(reverse=True)
    n = len(p_values)
    best_i = -1
    for i, pvalue in enumerate(p_values):
        if pvalue <= (((i + 1) * 0.05) / n):
            best_i = i
            break
    if best_i == -1:
        return None
    return p_values[i:]


def check_null_hypothesis(p_value_dict):
    pvalues_to_reject = bh_correction_for_pvalue(p_value_dict)
    phenotypes_to_reject = []
    if pvalues_to_reject == None:
        return None
    for phenotype in p_value_dict.keys():
        if p_value_dict[phenotype] in pvalues_to_reject:
            print(phenotype, p_value_dict[phenotype])
            phenotypes_to_reject.append(phenotype)
    return phenotypes_to_reject


def plot_PCA(normalized_data_sick_all_phen):
    pca = PCA(n_components=2)
    normalized_data_sick_all_phen = normalized_data_sick_all_phen.T
    print(normalized_data_sick_all_phen)
    projected = pca.fit_transform(normalized_data_sick_all_phen)
    fig, ax = plt.subplots()
    ax.set_xlabel('Principal Component 1', fontsize=15)
    ax.set_ylabel('Principal Component 2', fontsize=15)
    ax.scatter(projected[:, 0], projected[:, 1])
    plt.show()
