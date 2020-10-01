# -*- coding: utf-8 -*-
import os
import glob
import pandas as pd

envir = 'wolfpack'
project_name = 'thesis'
subproject_name = 'chapter_4'

home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/' + project_name + '/' + \
	subproject_name + '/'
func_dir = project_dir + 'scripts/functions/'
ref_dir = project_dir + 'refs/'
out_dir = project_dir + 'msigdbr_met_gene_articles/'
filename = "msigdbr_met_genes.txt"
key_terms = ["cancer", "tumour", "carcinoma", "sarcoma", "malignant", 
    "metastasis", "resistance", "oncogene", "proliferation", "apopto"]
bad_terms = ["abstract"]
no_returned_citations = 5

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# load scholar scraping function:
from identify_cancer_genes_functions import scholar_scrape

if os.path.isfile(ref_dir + filename):
    print('Loading genes ...\n')
    df = pd.read_csv(ref_dir + filename, sep=' ')
    query_genes = list(set(df['gene']))  
    # for each DE gene, scrape google scholar for articles:
    scholar_scrape(
        genes = query_genes,
        incl_terms = key_terms,
        rm_terms = bad_terms,
        no_returned = no_returned_citations,
        out_dir = out_dir
    )
else:
    print('Input gene file does not exist ...')


