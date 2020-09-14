# -*- coding: utf-8 -*-
import os
import glob
import pandas as pd

envir = 'wolfpack'
project_name = 'thesis'
subproject_name = 'Figure_2.2_individual_samples'
sample_name = "CID4463"
subcluster_method = "random_trees"
subcluster_p = "0.05"
coverage_filter = "filtered"
remove_artefacts = "artefacts_not_removed"
min_pct = 0.5
key_terms = ["cancer", "tumour", "carcinoma", "sarcoma", "malignant", 
    "metastasis", "resistance", "oncogene", "proliferation", "apopto"]
bad_terms = ["abstract"]
no_returned_citations = 5

if envir == 'wolfpack':
  home_dir = '/share/ScratchGeneral/jamtor/'
  project_dir = home_dir + 'projects/' + project_name + '/' + \
  	subproject_name + '/'
  results_dir = project_dir + 'results/'
  func_dir = project_dir + 'scripts/functions/'
  in_dir = results_dir + 'infercnv/' + sample_name + "/" + \
    coverage_filter + "/" + subcluster_method + "/p_" + subcluster_p + \
    "/" + remove_artefacts + "/tables/"
else:
    in_dir = '/Users/jamestorpy/Desktop/scrape_test/' + sample_name + \
    "/" + coverage_filter + "/" + subcluster_method + "/p_" + \
    subcluster_p + "/" + remove_artefacts, "/tables/"

# load scholar scraping function:
from identify_cancer_genes_functions import scholar_scrape

if os.path.isfile(in_dir + 'msigdb_cancer_subpop_DE_genes.txt'):
    print('Loading DE data ...\n')
    df = pd.read_csv(in_dir + 'msigdb_cancer_subpop_DE_genes.txt', sep=' ')
    DE_genes = list(set(df['gene']))  
    # for each DE gene, scrape google scholar for articles:
    scholar_scrape(
        sample_id = sample_name,
        genes = DE_genes,
        incl_terms = key_terms,
        rm_terms = bad_terms,
        no_returned = no_returned_citations,
        out_dir = in_dir
    )
else:
    print('DE data does not exist for ' + sample_name + ', skipping...')


