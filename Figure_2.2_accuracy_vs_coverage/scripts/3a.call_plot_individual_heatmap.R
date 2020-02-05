#!/bin/bash

project_name="single_cell/identify_epithelial"
subproject_name="brca_mini_atlas_131019"
ncores=6
include_t_cells="TRUE"
subset_data="FALSE"
subset_samples="TRUE"
na_colour="white"
include_normals="FALSE"
missing_genes_min_proportion="0.5"
cluster_by="unsupervised"
type_order="LumA_SC50.LumB_SC50.Her2_SC50.Basal_SC50.NA"
cluster_method="ward.D2"


home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

log_dir="$project_dir/logs/cluster_by_$cluster_by/"
mkdir -p $log_dir
echo "Logs are in $log_dir"
qsub -wd $log_dir -pe smp $ncores -N $cluster_by.combined_plot -b y -j y -V -P TumourProgression \
  "/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R \
  CMD BATCH  --no-save '--args $include_t_cells $subset_data $subset_samples $na_colour $include_normals $missing_genes_min_proportion $cluster_by $type_order $cluster_method' $script_dir/5b.plot_group_combined_heatmap.R"
#  /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/5b.plot_individual_heatmap.R $sample_name $include_t_cells $downsample $downsample_proportion