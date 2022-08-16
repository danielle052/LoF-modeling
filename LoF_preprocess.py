from ukbb_parser import create_dataset, create_ICD10_dataset, get_chrom_raw_marker_data, get_chrom_imputation_data
from ukbb_parser.shared_utils.util import summarize
import pandas as pd
import import_ipynb
from Gene_parser import parse_pwas_data_set, Gene
from firm.gene_dataset import get_or_create_gene_dataset
import multiprocessing, firm, geneffect
import numpy as np


def get_eid_fields(PHENOTYPE, threshold, phen_name):
    eid, fields, covariates = create_dataset(PHENOTYPE, nrows=None)
    eid = pd.DataFrame(eid)
    eid_fields = eid.join(fields.dropna())
    eid_fields[phen_name] = eid_fields[phen_name].fillna(0)
    eid_fields[phen_name] = eid_fields[phen_name].apply(lambda x: 1.0 if x > threshold else 0)
    summarize(eid)
    summarize(fields.dropna())
    summarize(eid_fields)
    return eid_fields, eid, fields


def get_eid_fields_icd10(PHENOTYPE, phen_name):
    eid, ICD10_tree, additional_fields, covariates, _ = create_ICD10_dataset(nrows=None)
    eid = pd.DataFrame(eid)

    #####################################

    ##filter on sex:

    #     eid['sex'] = covariates['sex']
    #     additional_fields['sex'] = covariates['sex']
    #     eid = eid[eid['sex'] == 0.0]
    #     additional_fields = additional_fields[additional_fields['sex'] == 0.0]
    #     eid = eid.drop(['sex'], axis=1)
    #     additional_fields = additional_fields.drop(['sex'], axis=1)
    #     summarize(eid)

    ########################################

    (_, ICD10_node), = ICD10_tree[ICD10_tree['coding'] == PHENOTYPE].iterrows()
    samples = list(ICD10_node['samples'])
    samples.sort()
    all_ind = np.arange(0, 502520)
    all_ind = pd.DataFrame(all_ind)
    samples_df = pd.DataFrame(samples)

    all_ind.columns = ['ind']
    samples_df.columns = ['ind']

    samples_df[phen_name] = 1
    dfinal = all_ind.merge(samples_df, on="ind", how='left')
    dfinal[phen_name] = dfinal[phen_name].fillna(0)
    dfinal = dfinal.drop(['ind'], axis=1)
    eid_fields = eid.join(dfinal)

    return eid_fields, eid, additional_fields


def get_data(PHENOTYPE, threshold, isICD10, phen_name):
    if isICD10:
        return get_eid_fields_icd10(PHENOTYPE, phen_name)
    return get_eid_fields(PHENOTYPE, threshold, phen_name)


def make_data_arrays():
    variants_per_chrom_lst, general_phens, G = [], [], []
    for i in range(1, 23):
        variants_per_chrom_tmp, general_phens_tmp, G_tmp = get_chrom_raw_marker_data(str(i))
        variants_per_chrom_lst.append(variants_per_chrom_tmp)
        general_phens.append(general_phens_tmp)
        G.append(G_tmp)
    return variants_per_chrom_lst, general_phens[0], G


def filter_variants(genes_dict, variants_per_chrom_lst):
    relevant_variants = None
    relevant_chrom = list(genes_dict.keys())
    for i in relevant_chrom:
        ith_chrome_genes = genes_dict[i]
        ith_chrome_variants = variants_per_chrom_lst[i - 1]
        for gene in ith_chrome_genes:
            gene_satrt = gene.loc_start
            gene_end = gene.loc_end
            gene_len = gene_end - gene_satrt
            ith_chrome_relevant_variants = ith_chrome_variants[
                ith_chrome_variants['pos'].between(gene.loc_start, gene.loc_end)].copy()
            ith_chrome_relevant_variants['relative_pos'] = ith_chrome_relevant_variants['pos'].apply(
                lambda x: ((x - gene_satrt) / gene_len))
            ith_chrome_relevant_variants['gene_name'] = gene.gene_name
            if relevant_variants is None:
                relevant_variants = ith_chrome_relevant_variants
            else:
                relevant_variants = pd.concat([relevant_variants, ith_chrome_relevant_variants])
    return relevant_variants


def process_variant(variant_record, geneffect_setup):
    chrom, pos, ref, alt = variant_record['chrom'], variant_record['pos'], variant_record['a0'], variant_record['a1']
    return geneffect_setup.variant_interpreter.process_variant(chrom, pos, ref, alt)


def get_cds_gene_effects(variant_record, geneffect_setup):
    chrom, pos, ref, alt = variant_record['chrom'], variant_record['pos'], variant_record['a0'], variant_record['a1']
    snp = geneffect_setup.variant_interpreter.process_variant(chrom, pos, ref, alt)
    if snp.cds_gene_effects == []:
        return snp.cds_gene_effects
    return snp.cds_gene_effects[0]


def create_cds_gene_effects_df(relevant_variants_df, geneffect_setup_lof):
    cds_gene_effects_df = relevant_variants_df.apply(get_cds_gene_effects, geneffect_setup=geneffect_setup_lof, axis=1)
    return cds_gene_effects_df


def get_firm_scores(relevant_variants_df, gene_dict):
    # prepare firm
    thread_pool = multiprocessing.Pool(4)
    geneffect_setup_lof = geneffect.Setup('hg19')
    firm.setup_uniprot_tracks(geneffect_setup_lof)
    firm_classifier_lof = firm.load_classifier(geneffect_setup_lof)
    firm_predict_batch_adjusted_proba_function = lambda missense_effects: firm_classifier_lof.predict_adjusted_proba(
        missense_effects, thread_pool=thread_pool)
    cds_gene_effects_df = create_cds_gene_effects_df(relevant_variants_df, geneffect_setup_lof)

    # create a list of dictionaries mapping the unique gene IDs (as defined by gene_indexer) to tuples of (effect_str, effect_score), where effect_str
    # is a string explaining a variant effect and effect_score is a variant's effect score. Each dictionary in the list should indicate the results of
    # the corresponding variant in the input list. 
    gene_dataset = get_or_create_gene_dataset('./', geneffect_setup_lof)
    uniprot_id_to_gene_index = {uniprot_id: index for index, uniprot_id in gene_dataset['uniprot_id'].iteritems()}
    gene_indexer = lambda gene: uniprot_id_to_gene_index[gene.uniprot_record.id]
    variants = relevant_variants_df.apply(process_variant, geneffect_setup=geneffect_setup_lof, axis=1)
    variants_gene_effects_and_scores = firm.determine_extended_gene_effects_and_scores(variants,
                                                                                       firm_predict_batch_adjusted_proba_function=firm_predict_batch_adjusted_proba_function,
                                                                                       gene_indexer=gene_indexer)

    # add unique identifier to each variant
    relevant_variants_df['cds_gene_effects'] = cds_gene_effects_df

    ## filter on variants that are not on the coding regions
    relevant_variants_df = relevant_variants_df[relevant_variants_df['cds_gene_effects'].str.len() != 0]

    # add firm scores to relevant_variants_df 
    cds_and_score = []
    for gene_dict in variants_gene_effects_and_scores:
        for variant_effect in gene_dict.values():
            cds = variant_effect[0][0]
            score = variant_effect[0][1]
            variant_lst = [cds, score]
            cds_and_score.append(variant_lst)

    cds_and_score_df = pd.DataFrame(cds_and_score, columns=['cds_gene_effects', 'firm_score'])
    dtype = dict(cds_gene_effects=str)

    relevant_variants_df_merged = relevant_variants_df.astype(dtype).merge(cds_and_score_df.astype(dtype), 'left')
    return relevant_variants_df_merged


def get_lof_no_effect_variants(variants_and_scores, firm_threshold):
    lof_variants = variants_and_scores[variants_and_scores['firm_score'] < firm_threshold]
    no_effect_variants = variants_and_scores[variants_and_scores['firm_score'] >= firm_threshold]
    return lof_variants, no_effect_variants


def get_list_of_genotypes_info(variants, eid, general_phens, genotyping_index, G):
    list_of_dfs = []
    for i, (_, variant) in enumerate(variants.iterrows()):
        print('Variant %d/%d' % (i + 1, len(variants)), end='\r')
        variant_genotypes = G[int(variant['chrom']) - 1][variant['i'], :].compute()

        sample_variant_genotypes = variant_genotypes[genotyping_index.values]
        sample_variant_genotypes = pd.DataFrame(sample_variant_genotypes, columns=[variant['snp']])
        list_of_dfs.append(sample_variant_genotypes)
    return list_of_dfs


def create_df_variants_genotypes(list_of_dfs, genotyping_index, eid_fields, general_phens):
    eid_fields['iid'] = eid_fields['eid']
    cols = general_phens[['i', 'iid']]
    variants_genotypes_df = pd.concat(list_of_dfs, axis=1)
    variants_genotypes_df['i'] = genotyping_index.values
    dtype = dict(i=int)
    variants_genotypes_fam_merged = variants_genotypes_df.astype(dtype).merge(cols, on='i', how="inner")
    dtype = dict(iid=int)
    all_variants_genotypes = variants_genotypes_fam_merged.astype(dtype).merge(eid_fields.astype(dtype), on='iid',
                                                                               how="left")
    print(all_variants_genotypes)
    return all_variants_genotypes


def make_sick_healthy_dfs(variants, phenotype):
    sick = variants[variants[phenotype] == 1.0]
    healthy = variants[variants[phenotype] == 0.0]
    sick = sick.drop(['i', 'eid', 'iid'], axis=1)
    healthy = healthy.drop(['i', 'eid'], axis=1)
    return sick, healthy


def save_dfs(dir_output_path, lof_variants, no_effect_variants, sick_lof, healthy_lof, sick_no_effect,
             healthy_no_effect):
    lof_variants.to_csv(dir_output_path + "/lof_variants_female.csv", index=False)
    no_effect_variants.to_csv(dir_output_path + "/no_effect_variants_female.csv", index=False)
    sick_lof.to_csv(dir_output_path + "/sick_lof_female.csv", index=False)
    healthy_lof.to_csv(dir_output_path + "/healthy_lof_female.csv", index=False)
    sick_no_effect.to_csv(dir_output_path + "/sick_no_effect_female.csv", index=False)
    healthy_no_effect.to_csv(dir_output_path + "/healthy_no_effect_female.csv", index=False)


def create_visualization_data(phenotype_name, isICD10, PHENOTYPE, phen_threshold, dir_output_path, genes_path,
                              gene_fdr_col, gene_fdr_thershold, firm_threshold=0.2):
    """
    The main function in order to create csv files for visualization.
    This function gets the following arguments and creates 6 csv files: a list of all LoF variants and their FIRM's
    scores, a list of the rest variants and their FIRM's scores, csv of sick people with the LoF variants for each one
    of them, csv of healthy people with LoF variants for each one of them, csv of healthy people with the LoF variants
    for each one of them and csv of healthy people with LoF variants for each one of them.
    :param phenotype_name: the relevant phenotype
    :param isICD10: true if this phenotype is ICD10's disease and false otherwise.
    :param PHENOTYPE: the relevant phenotype's id from the uk biobank
    :param phen_threshold: the threshold from it the sample is tagged as sick.
    :param dir_output_path: dir to save the csv output files
    :param genes_path: list of genes path
    :param gene_fdr_col: name of fdr column
    :param gene_fdr_thershold: fdr threshold
    :param firm_threshold: threshold for FIRM scores
    """
    eid_fields, eid, fields = get_data(PHENOTYPE, phen_threshold, isICD10, phenotype_name)
    variants_per_chrom_lst, general_phens, G = make_data_arrays()
    genes_dict = parse_pwas_data_set(genes_path, gene_fdr_thershold, gene_fdr_col)
    relevant_variants_df = filter_variants(genes_dict, variants_per_chrom_lst)
    variants_and_scores = get_firm_scores(relevant_variants_df, genes_dict)
    lof_variants, no_effect_variants = get_lof_no_effect_variants(variants_and_scores, firm_threshold)

    my_series = eid.squeeze()
    genotyping_index = my_series.map(general_phens.astype({'iid': int}).set_index('iid')['i'])

    dfs_lof_variants = get_list_of_genotypes_info(lof_variants, eid, general_phens, genotyping_index, G)
    dfs_no_effect_variants = get_list_of_genotypes_info(no_effect_variants, eid, general_phens, genotyping_index, G)

    genotypes_lof_variants = create_df_variants_genotypes(dfs_lof_variants, genotyping_index, eid_fields, general_phens)
    genotypes_no_effect_variants = create_df_variants_genotypes(dfs_no_effect_variants, genotyping_index, eid_fields,
                                                                general_phens)

    sick_lof, healthy_lof = make_sick_healthy_dfs(genotypes_lof_variants, phenotype_name)
    sick_no_effect, healthy_no_effect = make_sick_healthy_dfs(genotypes_no_effect_variants, phenotype_name)

    save_dfs(dir_output_path, lof_variants, no_effect_variants, sick_lof, healthy_lof, sick_no_effect,
             healthy_no_effect)