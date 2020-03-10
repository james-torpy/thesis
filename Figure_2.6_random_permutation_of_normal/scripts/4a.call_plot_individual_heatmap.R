#!/bin/bash

project_name="thesis"
subproject_name="Figure_2.3_min_length_and_gap"
ncores=6
sample_name="CID4520N_cancer_sim"
include_t_cells="TRUE"
#simulation_numbers=$(seq 1 30)
simulation_numbers=( "2" )
analysis_mode="samples"
neutral_signal_range="0.97_1.03"
nUMI_threshold="25000"
nGene_threshold="5000"
gap_or_CNV="gap"
CNV_type="loss"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

for simulation_number in ${simulation_numbers[@]}
  do echo $simulation_number

  echo "Plotting InferCNV and true positive/false negative results for $gap_or_CNV of type $CNV_type..."

  log_dir="$project_dir/logs/plot_infercnv/$gap_or_CNV/$simulation_number/$CNV_type/$analysis_mode.mode/"
  mkdir -p $log_dir
  echo "Logs are in $log_dir"
  qsub -wd $log_dir -pe smp $ncores -N plot.$gap_or_CNV.$simulation_number.$CNV_type.$analysis_mode.mode -b y -j y -V -P TumourProgression \
    "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
    CMD BATCH  --no-save \
    '--args \
    $sample_name \
    $include_t_cells \
    $simulation_number \
    $analysis_mode \
    $neutral_signal_range \
    $nUMI_threshold \
    $nGene_threshold \
    $gap_or_CNV \
    $CNV_type' \
    $script_dir/4b.plot_individual_heatmap.R"
#   /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/Rscript --vanilla $script_dir/4b.plot_individual_heatmap.R \
#    $sample_name \
#    $include_t_cells \
#    $simulation_number \
#    $analysis_mode \
#    $neutral_signal_range \
#    $nUMI_threshold \
#    $nGene_threshold \
#    $gap_or_CNV \
#    $CNV_type
done

# 
