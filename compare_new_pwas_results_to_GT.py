import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import import_ipynb
import Parse_open_target_data
from matplotlib_venn import venn3

pwas_control_paths = {'diabetes': pwas_diabetes_GT, 'asthma': pwas_asthma_GT, 'lung_cancer': pwas_lung_cancer_GT,
                      'hypertension': pwas_hypertension_GT}
sigmoid_hypothesis_paths = {'diabetes': sigmoid_diabetes_GT, 'asthma': sigmoid_asthma_GT,
                            'lung_cancer': sigmoid_lung_cancer_GT, 'hypertension': sigmoid_hypertension_GT}
relative_position_hypothesis_paths = {'diabetes': relative_position_diabetes_GT, 'asthma': relative_position_asthma_GT,
                                      'lung_cancer': relative_position_lung_cancer_GT,
                                      'hypertension': relative_position_hypertension_GT}
# pwas_control_paths = {'diabetes': pwas_diabetes_GT,'asthma': pwas_asthma_GT}
# sigmoid_hypothesis_paths = {'diabetes':sigmoid_diabetes_GT, 'asthma': sigmoid_asthma_GT}
# relative_position_hypothesis_paths = {'diabetes': relative_position_diabetes_GT,'asthma':relative_position_asthma_GT}

phenotypes = {}


# OT = open target.
class Phenotype():
    def __init__(self, phenotype_name):
        self.phenotype_name = phenotype_name
        self.OT_overall_score_GT = pd.DataFrame()
        self.OT_genetic_association_GT = pd.DataFrame()
        self.OT_litrature_score_GT = pd.DataFrame()
        self.pwas_control_results = pd.DataFrame()
        self.sigmoid_hypothesis = pd.DataFrame()
        self.relative_position_hypothesis = pd.DataFrame()

    def update_open_target_GT(self, OT_overall_score_df, OT_genetic_association_df, OT_litrature_df):
        self.OT_overall_score_GT = OT_overall_score_df
        self.OT_genetic_association_GT = OT_genetic_association_df
        self.OT_litrature_score_GT = OT_litrature_df

    def update_pwas_GT(self, pwas_control_df):
        self.pwas_control_results = pwas_control_df

    def update_hypothesis_data(self, sigmoid_df, relative_position_df):
        self.sigmoid_hypothesis = sigmoid_df
        self.relative_position_hypothesis = relative_position_df


# OT = open target.
def parse_GT_data():
    OT_overall_score_all_phens, OT_genetic_association_all_phens, OT_litrature_score_all_phens = Parse_open_target_data.create_GTs()
    for phenotype_name in OT_overall_score_all_phens.keys():
        # open target data
        cur_phenotype = Phenotype(phenotype_name)
        overall_score_GT = OT_overall_score_all_phens[phenotype_name]
        gwas_GT = OT_genetic_association_all_phens[phenotype_name]
        litrature_GT = OT_litrature_score_all_phens[phenotype_name]
        cur_phenotype.update_open_target_GT(overall_score_GT, gwas_GT, litrature_GT)

        # pwas data
        pwas_control_GT = pd.read_csv(pwas_control_paths[phenotype_name])
        cur_phenotype.update_pwas_GT(pwas_control_GT)

        # hypotheses data
        sigmoid_df = pd.read_csv(sigmoid_hypothesis_paths[phenotype_name])
        #         relative_position_df =  pd.read_csv(relative_position_hypothesis_paths[phenotype_name])
        relative_position_df = pd.DataFrame()  # delete
        cur_phenotype.update_hypothesis_data(sigmoid_df, relative_position_df)

        phenotypes[phenotype_name] = cur_phenotype


def plot_new_discovered_genes():
    for phenotype_name in phenotypes:
        phenotype = phenotypes[phenotype_name]
        sigmoid_genes = set(
            phenotype.sigmoid_hypothesis[phenotype.sigmoid_hypothesis['fdr_qval'] < 0.05]['symbol'].values)
        #         relative_position_genes = set(phenotype.relative_position_hypothesis[phenotype.relative_position_
        #         hypothesis['fdr_qval'] < 0.05]['symbol'].values)
        pwas_control_genes = set(
            phenotype.pwas_control_results[phenotype.pwas_control_results['fdr_qval'] < 0.05]['symbol'].values)
        litrature_OT_genes = set(
            phenotype.OT_litrature_score_GT[phenotype.OT_litrature_score_GT['literature_score'] > 0.45][
                'approvedSymbol'].values)
        venn3([sigmoid_genes, pwas_control_genes, litrature_OT_genes],
              ('Sigmoid hypothesis', 'PWAS control', 'Litrature'))
        plt.title(phenotype_name)
        plt.show()
        inter = sigmoid_genes.intersection(pwas_control_genes)
        print(pwas_control_genes.intersection(litrature_OT_genes))
        print(inter.intersection(litrature_OT_genes))
        print(sigmoid_genes.intersection(litrature_OT_genes))

        if phenotype_name == 'hypertension':
            print(phenotype.OT_overall_score_GT.shape[0])
        OT_overall_score_GT = set(phenotype.OT_overall_score_GT[phenotype.OT_overall_score_GT['total_score'] > 0.485][
                                      'approvedSymbol'].values)
        venn3([sigmoid_genes, pwas_control_genes, OT_overall_score_GT],
              ('Sigmoid hypothesis', 'Vanilla PWAS', 'Open target - overall score'))
        plt.title("Genes associated with hypertension:", x=0.37, y=1.05, fontweight="bold")
        plt.show()
        inter = sigmoid_genes.intersection(pwas_control_genes)

        print(inter.intersection(OT_overall_score_GT))
        print(sigmoid_genes.intersection(OT_overall_score_GT))


parse_GT_data()
plot_new_discovered_genes()
