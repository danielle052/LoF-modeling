import matplotlib.pyplot as plt
import pandas as pd
from ukbb_parser.shared_utils.util import summarize
from scipy.stats import chi2_contingency
import numpy as np
import glob
import os
import math
import ntpath


def parse_data():
    genes = glob.glob(os.path.join(input_dir_path, "*.csv"))
    #     genes = ['/cs/labs/michall/danielle.levi/pwas_results/variants_per_gene_with_relative_pos/1663.csv']
    return genes


def sigmoid_hypothesis(genes):
    for gene in genes:
        variants_with_relative_position = pd.read_csv(gene)
        val = variants_with_relative_position[variants_with_relative_position['effect_score'] < 0.2][
            'relative_pos'].dropna().apply(lambda x: 1 / (1 + math.exp(-x)))
        variants_with_relative_position['effect_score'].update(val)
        file_name = ntpath.basename(gene)
        file_path = output_dir_path_sigmoid + file_name
        variants_with_relative_position.drop(['relative_pos'], axis=1).to_csv(file_path, index=False)


def relative_pos_hypothesis(genes):
    for gene in genes:
        variants_with_relative_position = pd.read_csv(gene)
        val = variants_with_relative_position[variants_with_relative_position['effect_score'] < 0.2][
            'relative_pos'].dropna()
        variants_with_relative_position['effect_score'].update(val)
        file_name = ntpath.basename(gene)
        file_path = output_dir_path_relative_pos + file_name
        variants_with_relative_position.drop(['relative_pos'], axis=1).to_csv(file_path, index=False)


def five_bins_hypothesis(genes, new_scores):
    for gene in genes:
        variants_with_relative_position = pd.read_csv(gene)
        lof_variants = variants_with_relative_position[variants_with_relative_position['effect_score'] < 0.2][
            'relative_pos'].dropna()
        
        lof_variants_first_bin = lof_variants[lof_variants < 0.2].apply(lambda x: new_scores[0])
        if len(lof_variants_first_bin):
            variants_with_relative_position['effect_score'].update(lof_variants_first_bin)

        lof_variants_second_bin = lof_variants[(lof_variants >= 0.2) & (lof_variants < 0.4)].apply(lambda x: new_scores[1])
        if len(lof_variants_second_bin):
            variants_with_relative_position['effect_score'].update(lof_variants_second_bin)

        lof_variants_third_bin = lof_variants[(lof_variants >= 0.4) & (lof_variants < 0.6)].apply(lambda x: new_scores[2])
        if len(lof_variants_third_bin):
            variants_with_relative_position['effect_score'].update(lof_variants_third_bin)

        lof_variants_fourth_bin = lof_variants[(lof_variants >= 0.6) & (lof_variants < 0.8)].apply(lambda x: new_scores[3])
        if len(lof_variants_fourth_bin):
            variants_with_relative_position['effect_score'].update(lof_variants_fourth_bin)

        lof_variants_fith_bin = lof_variants[(lof_variants >= 0.8) & (lof_variants <= 1)].apply(lambda x: new_scores[4])
        if len(lof_variants_fith_bin):
            variants_with_relative_position['effect_score'].update(lof_variants_fith_bin)
            
        file_name = ntpath.basename(gene)
        file_path = output_dir_path_linear_five_bins_half + file_name
        variants_with_relative_position.drop(['relative_pos'], axis=1).to_csv(file_path, index=False)
        
        
def two_bins_hypothesis(genes):
    for gene in genes:
        variants_with_relative_position = pd.read_csv(gene)
        lof_variants = variants_with_relative_position[variants_with_relative_position['effect_score'] < 0.2][
            'relative_pos'].dropna()
        
        lof_variants_first_bin = lof_variants[lof_variants <= 0.2].apply(lambda x: 0)
        if len(lof_variants_first_bin):
            variants_with_relative_position['effect_score'].update(lof_variants_first_bin)

        lof_variants_second_bin = lof_variants[(lof_variants > 0.2)].apply(lambda x: 1)
        if len(lof_variants_second_bin):
            variants_with_relative_position['effect_score'].update(lof_variants_second_bin)
            
        file_name = ntpath.basename(gene)
        file_path = output_dir_path_linear_five_bins_half + file_name
        variants_with_relative_position.drop(['relative_pos'], axis=1).to_csv(file_path, index=False)


genes = parse_data()
# sigmoid_hypothesis(genes)
# relative_pos_hypothesis(genes)
# five_bins_hypothesis(genes, [0, 0.1, 0.2, 0.3, 0.4])
# five_bins_hypothesis(genes, [0, 0.3, 0.5, 0.7, 0.9])
# two_bins_hypothesis(genes)
