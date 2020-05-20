#!/bin/bash

project_name="thesis"
subproject_name="Figure_2.5_accuracy_vs_coverage"
sample_name="CID4520N_cancer_sim"
include_t_cells="TRUE"
downsample_proportions="0.05_0.1_0.15_0.2_0.3_0.4_0.5_0.6_0.7_0.8_0.9_no"
downsample_type="gene"
CNV_type="both"
simulation_number_start="1"
simulation_number_end="30"
analysis_mode="samples"

ncores=3

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

echo "Plotting accuracy results for downsampled $sample_name datasets..."
    
log_dir="$project_dir/logs/$sample_name/sample_mode/$simulation_number_start.to.$simulation_number_end/$downsample_type.downsampling"
mkdir -p $log_dir
echo "Logs are in $log_dir"
qsub -wd $log_dir -pe smp $ncores -N plotac.$downsample_type.$simulation_number_start.to.$simulation_number_end -b y -j y -V \
  "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
  CMD BATCH  --no-save \
  '--args \
  $sample_name \
  $include_t_cells \
  $downsample_proportions \
  $downsample_type \
  $CNV_type \
  $simulation_number_start \
  $simulation_number_end \
  $analysis_mode' \
  $script_dir/4b.plot_downsampled_cancer_sim.R"
#  /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/Rscript --vanilla $script_dir/4b.plot_downsampled_cancer_sim.R \
#  $sample_name \
#    $include_t_cells \
#    $downsample_proportions \
#    $downsample_type \
#    $CNV_type \
#    $simulation_number_start \
#    $simulation_number_end \
#    $analysis_mode' \
#    $script_dir

#-P TumourProgression 
