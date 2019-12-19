#!/bin/bash

sample_name=CID4520N

project_name="thesis"
subproject_name="Figure_2.3_simulations"
ncores=10
subset_data="FALSE"
include_t_cells="TRUE"
nUMI_threshold="25000"
nGene_threshold="5000"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

random_proportions=( "0" "0.01" "0.05" "0.1" "0.2" "0.3" "0.4" "0.5" "0.75" )

for random_proportion in ${random_proportions[@]}
	do echo $downsample_proportion
	if [ "$downsample_proportion" = "0" ]
      then echo "Running InferCNV for $sample_name..."
 
      log_dir="$project_dir/logs/$sample_name"
      mkdir -p $log_dir
      echo "Logs are in $log_dir"
      qsub -wd $log_dir -pe smp $ncores -N rp.$random_proportion.$sample_name -b y -j y -V -P TumourProgression \
        "/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R \
        CMD BATCH  --no-save '--args $sample_name $ncores $subset_data $include_t_cells $nUMI_threshold $nGene_threshold $random_proportion' $script_dir/3b.random_permutation_infercnv.R"
      #  /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/3b.random_permutation_infercnv.R $sample_name $ncores $subset_data $include_t_cells $nUMI_threshold $nGene_threshold $random_proportion
  
    else
 
      echo "Running InferCNV for $sample_name with $random_proportion of genes randomly permutated..."
    
      log_dir="$project_dir/logs/$sample_name/$random_proportion.randomly.permutated"
      mkdir -p $log_dir
      echo "Logs are in $log_dir"
      qsub -wd $log_dir -pe smp $ncores -N rp.$random_proportion.$sample_name -b y -j y -V -P TumourProgression \
        "/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R \
        CMD BATCH  --no-save '--args $sample_name $ncores $subset_data $include_t_cells $nUMI_threshold $nGene_threshold $random_proportion' $script_dir/3b.random_permutation_infercnv.R"
      #  /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/3b.random_permutation_infercnv.R $sample_name $ncores $subset_data $include_t_cells $nUMI_threshold $nGene_threshold $random_proportion
    fi
  printf "\n"
done
