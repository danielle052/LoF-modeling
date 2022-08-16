import pandas as pd

LOC_START_COL = "cds_start"
LOC_END_COL = "cds_end"
CHROM_COL = "chr"
GENE_NAME_COL = "symbol"


class Gene:
    def __init__(self, loc_start: int, loc_end: int, chromosome: int, gene_name: str):
        self.loc_start = loc_start
        self.loc_end = loc_end
        self.chromosome = chromosome
        self.gene_name = gene_name


def parse_pwas_data_set(file_path: str, treshhold: float, column_name: str):
    pwas_output = pd.read_csv(file_path)
    relevant_genes = pwas_output[pwas_output[column_name] < treshhold]  # check the column name 
    genes_dict = {}
    for index, gene_row in relevant_genes.iterrows():
        cur_gene = Gene(gene_row[LOC_START_COL], gene_row[LOC_END_COL], gene_row[CHROM_COL],
                        gene_row[GENE_NAME_COL])
        if gene_row[CHROM_COL] not in ["X", "Y"]:
            if int(gene_row[CHROM_COL]) not in genes_dict:
                genes_dict[int(gene_row[CHROM_COL])] = []
            genes_dict[int(gene_row[CHROM_COL])].append(cur_gene)
    return genes_dict
