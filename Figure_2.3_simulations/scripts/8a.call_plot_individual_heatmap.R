#!/bin/bash

sample_name=CID4520N_cancer_sim

project_name="thesis"
subproject_name="Figure_2.3_simulations"
ncores=6
include_t_cells="TRUE"
downsample="TRUE"
neutral_signal_range="0.97_1.03"
nUMI_threshold=25000
nGene_threshold=5000


home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

downsample_props=(0.05 0.15 $(seq 0.1 0.1 1))

for downsample_proportion in ${downsample_props[@]}
	do echo $downsample_proportion
  echo "Plotting InferCNV results for non-donsampled $sample_name..."
  log_dir="$project_dir/logs/$sample_name/$downsample_proportion.downsampling/"
  mkdir -p $log_dir
  echo "Logs are in $log_dir"
  qsub -wd $log_dir -pe smp $ncores -N plot.$downsample_proportion.$sample_name -b y -j y -V -P TumourProgression \
    "/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R \
    CMD BATCH  --no-save '--args $sample_name $include_t_cells $downsample $downsample_proportion $neutral_signal_range $nUMI_threshold $nGene_threshold' $script_dir/5b.plot_individual_heatmap.R"
  #  /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/5b.plot_individual_heatmap.R $sample_name $include_t_cells $downsample $downsample_proportion

  printf "\n"
done
