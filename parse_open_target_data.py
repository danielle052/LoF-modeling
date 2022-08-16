import pandas as pd
import requests
import json

DIABETES_TYPE1_EFO_ID = 'MONDO_0005147'
BREAST_CANCER_EFO_ID = 'MONDO_0007254'
COLORECTAL_CANCER_EFO_ID = 'EFO_0005842'
LUNG_CANCER_EFO_ID = 'MONDO_0008903'
MELANOMA_EFO_ID = 'EFO_0000756'
STROKE_EFO_ID = 'EFO_0000712'
HYPERTENSION_EFO_ID = 'EFO_0000537'
DIABETES_TYPE2_EFO_ID = 'MONDO_0005148'
MULTIPLE_SCLEROSIS_EFO_ID = 'MONDO_0005301'
ASTHMA_EFO_ID = 'MONDO_0004979'

DISEASE_GENE_QUERY = """
        query DiseaseAssociationsQuery($efoId: String!, $index: Int!, $size: Int!, $filter: String, $sortBy: String!, $aggregationFilters: [AggregationFilter!]) {
          disease(efoId: $efoId) {
            id
            associatedTargets(page: {index: $index, size: $size}, orderByScore: $sortBy, BFilter: $filter, aggregationFilters: $aggregationFilters) {
              count
              rows {
                target {
                  id
                  approvedSymbol
                  approvedName
                  __typename
                }
                score
                datatypeScores {
                  componentId: id
                  score
                  __typename
                }
                __typename
              }
              __typename
            }
            __typename
          }
        }
    """


def flatten_json(y):
    """
    Take a json and faltten it. 
    taken from the web
    """
    out = {}

    def flatten(x, name=''):

        # If the Nested key-value 
        # pair is of dict type
        if type(x) is dict:

            for a in x:
                flatten(x[a], name + a + '_')

        # If the Nested key-value
        # pair is of list type
        elif type(x) is list:

            i = 0

            for a in x:
                flatten(a, name + str(i) + '_')
                i += 1
        else:
            out[name[:-1]] = x

    flatten(y)
    return out


def ask_api(query_string, variables):
    """
    Recive query adn return a JSON string
    """

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})
    #     print(r.status_code)

    # Transform API response into JSON 
    return json.loads(r.text)


def get_gene_data_gene_association(target):
    """
    creates a sieries from a json of a SNP 
    """
    g_data = flatten_json(target['target'])
    g_data.update({'total_score': target['score']})
    for i in target['datatypeScores']:
        if i['componentId'] == 'genetic_association':
            g_data.update({'genetic_association_score': i['score']})
            continue
        if i['componentId'] == 'literature':
            g_data.update({'literature_score': i['score']})
            continue
        if i['componentId'] == 'known_drug':
            g_data.update({'known_drug_score': i['score']})
            continue
        if i['componentId'] == 'rna_expression':
            g_data.update({'rna_expression_score': i['score']})
            continue
        if i['componentId'] == 'animal_model':
            g_data.update({'animal_model_score': i['score']})
            continue
        if i['componentId'] == 'somatic_mutation':
            g_data.update({'somatic_mutation_score': i['score']})
            continue
        if i['componentId'] == 'affected_pathway':
            g_data.update({'affected_pathway_score': i['score']})
            continue
    return pd.Series(g_data).drop('__typename')


def get_gene_data_overall_association(target):
    """
    creates a sieries from a json of a SNP 
    """
    g_data = flatten_json(target['target'])
    if target['score'] in target.keys():
        g_data.update({'overall_association_score': i['score']})
    return pd.Series(g_data).drop('__typename')


def get_disease_targets(efoId: str, sort_by: str = 'genetic_association', size: int = 200,
                        threshold=0.0035) -> pd.DataFrame:
    """
    efoId - ID of disease. example: ENSG00000196208 (endometriosis)
    sort_by - default: 'genetic_association'. Can be 'drugs' or any other option from opentarget
    size - amoun t of gene retrived by the query
    threshold - genetic association minimal score
    returns pd.DataFrame of the first -size- gene associated to the disease specified by the -efoId- sorted by -sortby-
    
    """
    v_dict = {
        "efoId": efoId,  # disease ID
        "index": 0,
        "size": size,
        "sortBy": sort_by,
        "filter": "",
        "aggregationFilters": []
    }
    targets = ask_api(DISEASE_GENE_QUERY, v_dict)['data']['disease']['associatedTargets']['rows']
    #     print(targets)

    genetic_association_df = pd.concat([get_gene_data_gene_association(i) for i in targets], axis=1).T
    #     result = genetic_association_df[genetic_association_df.genetic_association_score >= threshold]
    return genetic_association_df


def parse_open_target_data_overall_association(pathology):
    genetic_association_df = get_disease_targets(pathology, 'genetic_association')
    literature_df = get_disease_targets(pathology, 'literature')
    known_drug_df = get_disease_targets(pathology, 'known_drug')
    rna_expression_df = get_disease_targets(pathology, 'rna_expression')
    animal_model_df = get_disease_targets(pathology, 'animal_model')
    somatic_mutation_df = get_disease_targets(pathology, 'somatic_mutation')
    affected_pathway_df = get_disease_targets(pathology, 'affected_pathway')

    all_genes_per_pathology = pd.concat(
        [genetic_association_df, literature_df, known_drug_df, rna_expression_df, animal_model_df, somatic_mutation_df,
         affected_pathway_df]).drop_duplicates()
    all_genes_per_pathology_sorted = all_genes_per_pathology.sort_values('total_score', ascending=False).reset_index(
        drop=True)
    return all_genes_per_pathology_sorted


def parse_open_target_data_genetic_association(pathology):
    genetic_association_df = get_disease_targets(pathology, 'genetic_association')
    return genetic_association_df


def create_GT_general_score():
    pathologies = {'diabetes': DIABETES_TYPE1_EFO_ID, 'asthma': ASTHMA_EFO_ID, 'lung_cancer': LUNG_CANCER_EFO_ID,
                   'hypertension': HYPERTENSION_EFO_ID}
    GT_overall_score = {}
    GT_genetic_association = {}
    GT_litrature_score = {}

    for pathology in pathologies.keys():
        all_genes_per_pathology_overall = parse_open_target_data_overall_association(pathologies[pathology])
        all_genes_per_pathology_genetic = parse_open_target_data_genetic_association(pathologies[pathology])
        all_genes_per_pathology_litrature = parse_open_target_data_genetic_association(pathologies[pathology])

        GT_overall_score[pathology] = all_genes_per_pathology_overall
        GT_genetic_association[pathology] = all_genes_per_pathology_genetic
        GT_litrature_score[pathology] = all_genes_per_pathology_litrature

    return GT_overall_score, GT_genetic_association, GT_litrature_score


def create_GTs():
    GT_overall_score, GT_genetic_association, GT_litrature_score = create_GT_general_score()
    return GT_overall_score, GT_genetic_association, GT_litrature_score
