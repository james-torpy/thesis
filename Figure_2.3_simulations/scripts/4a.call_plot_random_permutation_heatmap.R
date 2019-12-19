#!/bin/bash

sample_name=CID4520N_cancer_sim

project_name="thesis"
subproject_name="Figure_2.3_simulations"
ncores=6
include_t_cells="TRUE"


home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

random_props=( 0 0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.75 )

for random_proportion in ${random_props[@]}
	do echo $random_proportion
  echo "Plotting InferCNV results for permutated $sample_name..."
  log_dir="$project_dir/logs/$sample_name/$random_proportion.permutated/"
  mkdir -p $log_dir
  echo "Logs are in $log_dir"
  qsub -wd $log_dir -pe smp $ncores -N plot.$downsample_proportion.$sample_name -b y -j y -V -P TumourProgression \
    "/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R \
    CMD BATCH  --no-save '--args $sample_name $include_t_cells $random_proportion' $script_dir/4b.plot_random_permutation_heatmap.R"
  #  /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/4b.plot_random_permutation_heatmap.R $sample_name $include_t_cells $random_proportion

  printf "\n"
done
