import matplotlib.pyplot as plt
import pandas as pd
from ukbb_parser.shared_utils.util import summarize
from scipy.stats import chi2_contingency
import numpy as np
import import_ipynb
import LOF_tool_statistical_validation as stat


def get_data(path_lof_variants, path_missense_variants, path_lof_sick_population, path_missense_sick_population,
             path_lof_healthy_population, path_missense_healthy_population):
    # Create dataframes
    lof_variants = pd.read_csv(path_lof_variants)
    missense_variants = pd.read_csv(path_missense_variants)
    lof_sick_population = pd.read_csv(path_lof_sick_population)
    missense_sick_population = pd.read_csv(path_missense_sick_population)
    lof_healthy_population = pd.read_csv(path_lof_healthy_population)
    missense_healthy_population = pd.read_csv(path_missense_healthy_population)
    return lof_variants, missense_variants, lof_sick_population, missense_sick_population, lof_healthy_population, missense_healthy_population


# parse the data before plotting
def parse_data(population_df):
    #     population_df = population_df.replace(2.0, 1.0)
    # Count poeple with variants
    population_df.loc['count_with_variant'] = population_df[population_df > 0].sum(axis=0)

    # count the number of nan values in each column
    population_df.loc['count_nan'] = population_df.isna().sum(axis=0)

    # Count relative population with variant
    total_population = population_df.shape[0]

    population_df.loc['count_population'] = population_df.loc['count_with_variant'] / (
                2 * (total_population - population_df.loc['count_nan']))
    #     print(population_df['rs2229094'])

    # tuen count_population row into a column
    population_df = population_df.loc[['count_population'], :].T

    # create a snp column for later merge with variants table
    population_df['snp'] = population_df.index

    return population_df


def create_list_of_genes(variants_df):
    return variants_df['gene_name'].unique()


# filter the data frame acording to the relevant genes (and divide it into chuncks of relative positions?)
def filter_variants(gene_name, variants):
    variants = variants[variants['gene_name'] == gene_name]
    return variants


# this function prepare the data sets for visualization and call the functions that create the visualization
def create_visualization_per_gene(lof_variants, missense_variants, lof_sick_population, missense_sick_population,
                                  lof_healthy_population, missense_healthy_population, list_of_genes, pathology):
    df_agg_sick_per_phen = []  # todo change - create this data frame elsewhere
    df_agg_healty_per_phen = []  # todo change - create this data frame elsewhere
    for gene in list_of_genes:
        lof_variants_per_gene = filter_variants(gene, lof_variants)
        missense_variants_per_gene = filter_variants(gene, missense_variants)

        # merge variants with relative number of population
        lof_per_gene_with_sick_count = lof_variants_per_gene.merge(lof_sick_population, on='snp', how='inner')
        missense_per_gene_with_sick_count = missense_variants_per_gene.merge(missense_sick_population, on='snp',
                                                                             how='inner')
        lof_per_gene_with_healthy_count = lof_variants_per_gene.merge(lof_healthy_population, on='snp', how='inner')
        missense_per_gene_with_healthy_count = missense_variants_per_gene.merge(missense_healthy_population, on='snp',
                                                                                how='inner')

        # plot data per gene - scatter
        #         scatter_visualization_per_gene(gene, pathology, lof_per_gene_with_sick_count,
        #         missense_per_gene_with_sick_count, lof_per_gene_with_healthy_count, missense_per_gene_with_healthy_count)

        # plot data per gene, while gene is splitted to ten segments.
        df_agg_sick = create_aggregated_data_per_gene(lof_per_gene_with_sick_count, 10, "sick")
        df_agg_healty = create_aggregated_data_per_gene(lof_per_gene_with_healthy_count, 10, "healthy")

        # plot data per gene - bar
        bar_visualization_per_gene(gene, pathology, df_agg_sick, df_agg_healty)

        # todo change - create this data frame elsewhere
        sick_norm_sum = df_agg_sick["sick"].sum()
        healthy_norm_sum = df_agg_healty["healthy"].sum()
        if sick_norm_sum:
            df_agg_sick_per_phen.append(df_agg_sick / sick_norm_sum)
        if sick_norm_sum:
            df_agg_healty_per_phen.append(df_agg_healty / healthy_norm_sum)

    combined_df_sick_per_phen = pd.concat(df_agg_sick_per_phen, axis=1)
    combined_df_healthy_per_phen = pd.concat(df_agg_healty_per_phen, axis=1)
    return combined_df_sick_per_phen, combined_df_healthy_per_phen


def scatter_visualization_per_gene(gene, pathology, lof_per_gene_with_sick_count, missense_per_gene_with_sick_count,
                                   lof_per_gene_with_healthy_count, missense_per_gene_with_healthy_count):
    # plot data per gene
    fig, ax = plt.subplots()
    ax.scatter(lof_per_gene_with_sick_count['relative_pos'], (lof_per_gene_with_sick_count['count_population'] * 100),
               label='lof sick count')

    #         ax.scatter(missense_per_gene_with_sick_count['relative_pos'], \
    #                 missense_per_gene_with_sick_count['count_population'], label = 'missene sick count')

    ax.scatter(lof_per_gene_with_healthy_count['relative_pos'],
               (lof_per_gene_with_healthy_count['count_population'] * 100), label='lof healthy count', alpha=0.5)

    #         ax.scatter(missense_per_gene_with_healthy_count['relative_pos'], \
    #                 missense_per_gene_with_healthy_count['count_population'], label = 'missene healthy count', alpha=0.5)

    ax.legend()
    ax.set_title("Phenotype: " + pathology + ", " + "Gene name: " + gene)
    plt.xlabel("Variants' relative genomic position")
    plt.ylabel("Population with variants (%)")
    plt.xlim([-0.05, 1.05])
    plt.xticks(rotation=60)
    plt.show()


def bar_visualization_per_gene(gene, pathology, df_agg_sick, df_agg_healty):
    combined_df = pd.concat([df_agg_sick, df_agg_healty], axis=1)
    combined_df.plot.bar()

    plt.title("Phenotype: " + pathology + ", " + "Gene name: " + gene)
    plt.xlabel("Variants' relative genomic position")
    plt.ylabel("Population with variants (%)")
    plt.xticks(rotation=60)
    plt.show()


# This function creates aggregated data frame per regions in the gene.
# The population normalized count of variants is normalized by the sum of variants.
def create_aggregated_data_per_gene(lof_per_gene, period: int, count_population_col: str):
    lof_per_gene["variants_sum"] = 1
    cut_bins = pd.interval_range(start=0.0, end=1.00, periods=10)

    lof_per_gene["bins"] = pd.cut(lof_per_gene["relative_pos"], bins=pd.IntervalIndex(
        [pd.Interval(round(i.left, 1), round(i.right, 1), i.closed) for i in cut_bins]))
    df_agg = lof_per_gene.groupby("bins")["count_population", "variants_sum"].sum()  ## reset_index()
    df_agg[count_population_col] = (df_agg["count_population"] / df_agg["variants_sum"]) * 100
    df_agg.fillna(0.0, inplace=True)
    df_agg = df_agg.drop(columns=['count_population', 'variants_sum'])
    return df_agg


def plot_all_genes(all_variants, sick_population, healthy_population, name, pathology):
    variants_per_gene_with_sick_count = all_variants.merge(sick_population, on='snp', how='inner')
    variants_per_gene_with_healthy_count = all_variants.merge(healthy_population, on='snp', how='inner')

    fig, ax = plt.subplots()
    ax.scatter(variants_per_gene_with_sick_count['relative_pos'],
               (variants_per_gene_with_sick_count['count_population'] * 100), label='lof sick count', marker='^')

    ax.scatter(variants_per_gene_with_healthy_count['relative_pos'],
               (variants_per_gene_with_healthy_count['count_population'] * 100), label='lof healthy count', alpha=0.5,
               marker='o')

    ax.legend()
    ax.set_title("Phenotype: " + pathology + ", " + name + "all genes")
    plt.xlim([-0.05, 1.05])
    plt.xticks(rotation=60)
    plt.xlabel("Variants' relative genomic position")
    plt.ylabel("Population with variants (%)")
    plt.show()

    # plot data for all genes, while gene is splitted to ten segments.
    df_agg_sick = create_aggregated_data_per_gene(variants_per_gene_with_sick_count, 10, "sick")
    df_agg_healty = create_aggregated_data_per_gene(variants_per_gene_with_healthy_count, 10, "healthy")
    combined_df = pd.concat([df_agg_sick, df_agg_healty], axis=1)
    combined_df.plot.bar()
    plt.xticks(rotation=60)
    plt.title("Phenotype: " + pathology + ", " + "all genes")
    plt.xlabel("Variants' relative genomic position")
    plt.ylabel("Population with variants (%)")


def run_chi_square_test(all_variants, population_df_sick, population_df_healthy, healthy_population):
    sick_population = population_df_sick.loc[['count_with_variant'], :].T
    sick_population['snp'] = sick_population.index

    healthy_population = population_df_healthy.loc[['count_with_variant'], :].T
    healthy_population['snp'] = healthy_population.index

    #     all_variants = all_variants[all_variants['gene_name'] != 'IFNA7']
    #     all_variants = all_variants[all_variants['gene_name'] != 'ZNF16']

    variants_merged = all_variants.merge(sick_population, on='snp', how='inner')
    variants_with_count = variants_merged.merge(healthy_population, on='snp', how='inner')
    variants_with_count_relevant = variants_with_count[['count_with_variant_x', 'count_with_variant_y']]
    variants_with_count_relevant = variants_with_count_relevant.rename(
        columns={"relative_pos": "variant_relative_pos", "count_with_variant_x": "count_sick_population_with_variant",
                 "count_with_variant_y": "count_healthy_population_with_variant"})

    return variants_with_count_relevant


# pathologies_list = ["Colorectal_cancer", "Breast_cancer", "Lung_cancer", "Hypertension", "Melanoma", "Type1_diabetes", "Asthma"]
pathologies_list = ["Type1_diabetes"]
folder_path = "/cs/labs/michall/noaraha/LOF/"


def run_program(pathologies_list, folder_path):
    p_value_dict = {}
    normalized_data_sick = []
    normalized_data_healthy = []
    for pathology in pathologies_list:
        print("############ current pathology : " + pathology + " #############")
        path_lof_variants = folder_path + pathology + "/lof_variants.csv"
        path_missense_variants = folder_path + pathology + "/no_effect_variants.csv"
        path_lof_sick_population = folder_path + pathology + "/sick_lof.csv"
        path_missense_sick_population = folder_path + pathology + "/sick_no_effect.csv"
        path_lof_healthy_population = folder_path + pathology + "/healthy_lof.csv"
        path_missense_healthy_population = folder_path + pathology + "/healthy_no_effect.csv"

        # create data frames
        lof_variants, missense_variants, lof_sick_population, missense_sick_population, lof_healthy_population, missense_healthy_population = get_data(
            path_lof_variants, path_missense_variants, path_lof_sick_population, path_missense_sick_population,
            path_lof_healthy_population, path_missense_healthy_population)

        # parse data
        lof_sick_population_chi = lof_sick_population
        lof_healthy_population_chi = lof_healthy_population

        ##############################################################
        lof_sick_population = parse_data(lof_sick_population)
        missense_sick_population = parse_data(missense_sick_population)
        lof_healthy_population = parse_data(lof_healthy_population)
        missense_healthy_population = parse_data(missense_healthy_population)
        print("lof_sick_population:", lof_sick_population.shape[0])
        print("lof_healthy_population:", lof_healthy_population.shape[0])
        # create list of genes per pathology
        list_of_genes = create_list_of_genes(lof_variants)

        # plot data per gene
        combined_df_sick_per_phen, combined_df_healthy_per_phen = create_visualization_per_gene(lof_variants,
                                                                                                missense_variants,
                                                                                                lof_sick_population,
                                                                                                missense_sick_population,
                                                                                                lof_healthy_population,
                                                                                                missense_healthy_population,
                                                                                                list_of_genes,
                                                                                                pathology)

        # plot data for all genes
        plot_all_genes(lof_variants, lof_sick_population, lof_healthy_population, 'lof - ', pathology)
        #         plot_all_genes(missense_variants, missense_sick_population, missense_healthy_population, 'missense - ')

        # create dataset for statistic tests
        data_for_test = run_chi_square_test(lof_variants, lof_sick_population_chi, lof_healthy_population_chi,
                                            lof_healthy_population)
        observations = data_for_test.to_numpy()
        chi2, p_value, dof, ex = chi2_contingency(observations[~np.all(observations == 0, axis=1)], correction=False)

        # create p-value dictionary
        p_value_dict[pathology] = p_value

        # create data for PCA
        combined_df_sick_per_phen.columns += '_' + pathology
        combined_df_healthy_per_phen.columns += '_' + pathology
        normalized_data_sick.append(combined_df_sick_per_phen)
        normalized_data_healthy.append(combined_df_healthy_per_phen)
    normalized_data_sick_all_phen = pd.concat(normalized_data_sick, axis=1)
    normalized_data_healthy_all_phen = pd.concat(normalized_data_healthy, axis=1)
    return p_value_dict, normalized_data_sick_all_phen, normalized_data_healthy_all_phen


p_value_dict, normalized_data_sick_all_phen, normalized_data_healthy_all_phen = run_program(pathologies_list,
                                                                                            folder_path)
phenotypes_to_reject = stat.check_null_hypothesis(p_value_dict)
