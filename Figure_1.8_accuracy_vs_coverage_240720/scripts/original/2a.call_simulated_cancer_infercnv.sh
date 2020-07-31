#!/bin/bash

project_name="thesis"
subproject_name="Figure_2.2_accuracy_vs_coverage"

sample_name=CID4520N_cancer_sim
ncores=10
subset_data="FALSE"
include_t_cells="TRUE"
analysis_mode="samples"
downsample="TRUE"
downsample_proportions=( "no" "0.9" "0.8" "0.7" "0.6" "0.5" "0.4" "0.3" "0.2" "0.15" "0.1" "0.05" )
CNV_type="both"
#simulation_numbers=$(seq 11 30)

simulation_numbers=( "1" )
#downsample_proportions=( "no" )

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

for simulation_number in ${simulation_numbers[@]}
  do echo $simulation_number
    for downsample_proportion in ${downsample_proportions[@]}
    	do echo $downsample_proportion
     
        echo "Running InferCNV for $sample_name downsampled to $downsample_proportion"
      
        log_dir="$project_dir/logs/$sample_name/$simulation_number/sample_mode/$downsample_proportion.downsampling"
        mkdir -p $log_dir
        echo "Logs are in $log_dir"
        qsub -wd $log_dir -pe smp $ncores -N no.$simulation_number.ds.$downsample_proportion.$sample_name -b y -j y -V -P TumourProgression \
          "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
          CMD BATCH  --no-save '--args $sample_name $ncores $subset_data $include_t_cells $analysis_mode $downsample $downsample_proportion $CNV_type $simulation_number' $script_dir/2b.simulated_cancer_infercnv.R"
        # /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/4b.simulated_cancer_infercnv.R $sample_name $ncores $subset_data $include_t_cells $analysis_mode $downsample_proportion $downsample_type
    done
done

# 