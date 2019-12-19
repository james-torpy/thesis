#!/bin/bash

sample_name=CID4520N_cancer_sim

project_name="thesis"
subproject_name="Figure_2.3_simulations"
ncores=10
subset_data="FALSE"
include_t_cells="TRUE"
analysis_mode="samples"
downsample="TRUE"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

#downsample_props=($(seq 0.1 0.1 1))
downsample_props=( 0.05 0.15)

for downsample_proportion in ${downsample_props[@]}
	do echo $downsample_proportion
	if [ "$downsample_proportion" = "1.0" ]
      then echo "Running InferCNV for $sample_name..."
 
      log_dir="$project_dir/logs/$sample_name"
      mkdir -p $log_dir
      echo "Logs are in $log_dir"
      qsub -wd $log_dir -pe smp $ncores -N ds.$downsample_proportion.$sample_name -b y -j y -V -P TumourProgression \
        "/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R \
        CMD BATCH  --no-save '--args $sample_name $ncores $subset_data $include_t_cells $analysis_mode $downsample $downsample_proportion' $script_dir/4b.simulated_cancer_infercnv.R"
      #  /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/4b.simulated_cancer_infercnv.R $sample_name $ncores $subset_data $include_t_cells $analysis_mode $downsample_proportion
  
    else
 
      echo "Running InferCNV for $sample_name downsampled to $downsample_proportion"
    
      log_dir="$project_dir/logs/$sample_name/$downsample_proportion.downsampling"
      mkdir -p $log_dir
      echo "Logs are in $log_dir"
      qsub -wd $log_dir -pe smp $ncores -N ds.$downsample_proportion.$sample_name -b y -j y -V -P TumourProgression \
        "/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R \
        CMD BATCH  --no-save '--args $sample_name $ncores $subset_data $include_t_cells $analysis_mode $downsample $downsample_proportion' $script_dir/4b.simulated_cancer_infercnv.R"
      #  /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/4b.simulated_cancer_infercnv.R $sample_name $ncores $subset_data $include_t_cells $analysis_mode $downsample_proportion
    fi
  printf "\n"
done
