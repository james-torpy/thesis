# -*- coding: utf-8 -*-
import os
import glob
import pandas as pd

envir = 'wolfpack'
project_name = 'thesis'
subproject_name = 'Figure_2.2_individual_samples'
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
  in_path = results_dir + 'infercnv/'
  out_dir = in_path + 'DE_gene_articles/'
else:
    in_path = '/Users/jamestorpy/Desktop/scrape_test/'
    out_dir = in_path

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# load scholar scraping function:
from identify_cancer_genes_functions import scholar_scrape

# fetch sample names:
#sample_names = ["common_DE", os.path.basename(x) for x in glob.glob(in_path + 'CID*')]
sample_names = ["common_DE", "CID45171", "CID4463"]

# for each sample, fetch list of subpop DE genes:
all_DE_genes = {}

for sample_name in sample_names:
    in_dir = in_path + sample_name + "/" + coverage_filter + "/" + \
    	subcluster_method + "/p_" + subcluster_p + "/" + remove_artefacts + '/'
    if sample_name == 'common_DE':
        if os.path.isfile(in_dir + str(min_pct) + '_min_pct' + '/tables/sig_subpop_DE.txt'):
            print('Loading DE data for ' + sample_name + '...\n')
            df = pd.read_csv(in_dir + str(min_pct) + '_min_pct' + '/tables/sig_subpop_DE.txt', sep=' ')
        else:
            print('DE data does not exist for ' + sample_name + ', skipping...')
    else:
        if os.path.isfile(in_dir + 'tables/sig_subpop_DE.txt'):
            df = pd.read_csv(in_dir + 'tables/sig_subpop_DE.txt', sep=' ')
        else:
            print('DE data does not exist for ' + sample_name + ', skipping...')
    DE_genes = list(df['gene'])
    all_DE_genes[sample_name] = DE_genes
    
# for each DE gene, scrape google scholar for articles:
for sample_name in sample_names:
    scholar_scrape(
        sample_id = sample_name,
        genes = all_DE_genes[sample_name],
        incl_terms = key_terms,
        rm_terms = bad_terms,
        no_returned = no_returned_citations,
        out_dir = out_dir
    )
