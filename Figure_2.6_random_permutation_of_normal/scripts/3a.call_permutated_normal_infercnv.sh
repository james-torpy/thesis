#!/bin/bash

project_name="thesis"
subproject_name="Figure_2.3_min_length_and_gap"

sample_name=CID4520N_cancer_sim
numcores=10
include_t_cells="TRUE"
analysis_mode="samples"
gap_or_CNV="gap"
CNV_type="loss"
#simulation_numbers=$(seq 11 30)
simulation_numbers=( "1" )

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

for simulation_number in ${simulation_numbers[@]}
  do echo "Running InferCNV for $sample_name with added $gap_or_CNV\s..."
      
    log_dir="$project_dir/logs/infercnv/$gap_or_CNV/$CNV_type/$simulation_number/$analysis_mode.mode/"
    mkdir -p $log_dir
    echo "Logs are in $log_dir"
    qsub -wd $log_dir -pe smp $numcores -N infercnv.$CNV_type.$gap_or_CNV.$sample_name -b y -j y -V -P TumourProgression \
      "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
      CMD BATCH  --no-save '--args \
      $sample_name \
      $numcores \
      $include_t_cells \
      $analysis_mode \
      $gap_or_CNV \
      $CNV_type \
      $simulation_number' \
      $script_dir/3b.simulated_cancer_infercnv.R"
    # /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/3b.simulated_cancer_infercnv.R 
#        $sample_name \
#        $numcores \
#        $include_t_cells \
#        $analysis_mode \
#        $gap_or_CNV \
#        $CNV_type \
#        $simulation_number' \
#        $script_dir/2b.filtered_infercnv.R"
done;

# 