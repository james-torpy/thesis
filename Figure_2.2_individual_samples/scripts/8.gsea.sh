
project_name="thesis"
subproject_name="Figure_2.2_individual_samples"
sample_name="CID4463"
subcluster_method="random_trees"
subcluster_p="0.05"
coverage_filter="filtered"
remove_artefacts="artefacts_not_removed"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
results_dir="$project_dir/results"

in_path="$results_dir/infercnv/$sample_name/$coverage_filter/$subcluster_method/p_$subcluster_p/"
out_dir="$in_path/gsea"
input_dir="$out_dir/input_files"

gsea_dir="$home_dir/local/lib/GSEA_4.1.0"

# load java:
module load shacar/java/jdk-11.0.2

gsea_command="$gsea_dir/gsea-cli.sh GSEA \
  -res $input_dir/$sample_name\.txt \
  -cls $input_dir/phenotype.cls \
  #population_2_versus_population_1 \
  -gmx ftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v7.1.symbols.gmt \
  -collapse Collapse \
  -mode Max_probe \
  -norm meandiv \
  -nperm 1000 \
  -permute phenotype \
  -rnd_type no_balance \
  -scoring_scheme weighted \
  -rpt_label my_analysis \
  -metric Signal2Noise \
  -sort real \
  -order descending \
  -chip ftp.broadinstitute.org://pub/gsea/annotations_versioned/Human_NCBI_Entrez_Gene_ID_MSigDB.v7.1.chip \
  -create_gcts false \
  -create_svgs false \
  -include_only_symbols true \
  -make_sets true \
  -median false \
  -num 100 \
  -plot_top_x 20 \
  -rnd_seed timestamp \
  -save_rnd_lists false \
  -set_max 500 \
  -set_min 15 \
  -zip_report false \
  -out $out_dir"


$gsea_dir/gsea-cli.sh GSEA -res /share/ScratchGeneral/jamtor/projects/thesis/Figure_2.2_individual_samples/results/infercnv/CID4463/filtered/random_trees/p_0.05//gsea/input_files/CID4463.txt -cls /share/ScratchGeneral/jamtor/projects/thesis/Figure_2.2_individual_samples/results/infercnv/CID4463/filtered/random_trees/p_0.05//gsea/input_files/phenotype.cls -gmx ftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v7.1.symbols.gmt -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute phenotype -rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis -metric Signal2Noise -sort real -order descending -chip ftp.broadinstitute.org://pub/gsea/annotations_versioned/Human_NCBI_Entrez_Gene_ID_MSigDB.v7.1.chip -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out /share/ScratchGeneral/jamtor/projects/thesis/Figure_2.2_individual_samples/results/infercnv/CID4463/filtered/random_trees/p_0.05//gsea




