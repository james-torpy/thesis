
import glob 
import os
import pandas as pd

project_name = 'thesis'
subproject_name = 'Figure_2.2_individual_samples'
subcluster_method = "random_trees"
subcluster_p = "0.05"
coverage_filter = "filtered"
remove_artefacts = "artefacts_not_removed"

home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/' + project_name + '/' + \
	subproject_name + '/'
results_dir = project_dir + 'results/'
in_path = results_dir + 'infercnv/'

# fetch sample names:
sample_names = [os.path.basename(x) for x in glob.glob(in_path + 'CID*')]

# for each sample, fetch list of subpop DE genes:
for sample_name in sample_names:
    in_dir = in_path + sample_name + "/" + coverage_filter + "/" + \
    subcluster_method + "/p_" + subcluster_p + "/" + remove_artefacts + \
    "/tables/"
    if os.path.isfile(in_dir + '/sig_subpop_DE.txt'):
        print('Loading DE data for ' + sample_name + '...')
        df = pd.read_csv(in_dir + '/sig_subpop_DE.txt', sep=' ')
        DE_genes = list(df['gene'])
    else:
        print('DE data does not exist for ' + sample_name + ', skipping...')
