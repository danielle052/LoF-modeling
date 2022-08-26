# LoF-modeling
This repository contains:
1) LoF_preprocess.py- the script that creates the data for LoF visualizations.
2) LoF_tool_statistical_validation.py- in this script, we test if LoF relative positions on the gene correlate with phenotypes.
3) LoF_visualization.py- visualizes different LoF variants.
4) compare_new_pwas_results_to_GT.py- compares the new PWAS results (after the scores have been updated) with the ground truth.
5) gene_parser.py- creates a dictionary of all genes associated with a specific disease according to PWAS.
6) parse_open_target_data.py- script for retrieving and parsing open targets data.
7) update_firm_score.py- updates FIRM scores by different hypotheses. 

How to get visualization of LoF mutations:
1) Clone this repository
2) Run create_visualization_data(phenotype_name, isICD10, PHENOTYPE, phen_threshold, dir_output_path, genes_path, gene_fdr_col, gene_fdr_thershold, firm_threshold=0.2)    from LoF_preprocess.py with the relevant parameters according to the function ducmentation.
3) Run run_program(pathologies_list, folder_path) from LoF_visualization.py. make sure that the folder_path is the same as dir_output_path.
4) For statisiacal tests use LoF_tool_statistical_validation.py with the relevant function.

Running PWAS with new hypotheses for scoring LoF mutations:
1) Clone this repository.
2) Follow steps 1-2 in the README in this repository: https://github.com/nadavbra/pwas
3) Choose hypothesis from update_firm_score.py (or add a new hypothesis with the same contruct to this file) and run the relevant function.
4) Follow steps 3-4 in the README in this repository: https://github.com/nadavbra/pwas
5) To compare the results to the Open Targets data run plot_new_discovered_genes() from compare_new_pwas_results_to_GT.py
