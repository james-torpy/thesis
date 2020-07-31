#!/bin/bash

project_name="thesis"
subproject_name="Figure_2.3_min_length_and_gap"

sample_name=CID4520N
ncores=10
nUMI_threshold=25000
nGene_threshold=5000
include_t_cells="TRUE"
analysis_mode="samples"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

echo "Running InferCNV on $sample_name..."
      
log_dir="$project_dir/logs/$sample_name/$analysis_mode.mode"
mkdir -p $log_dir
echo "Logs are in $log_dir"
qsub -wd $log_dir -pe smp $ncores -N $sample_name.infercnv -b y -j y -V -P TumourProgression \
  "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
  CMD BATCH  --no-save '--args \
  $sample_name \
  $ncores \
  $nUMI_threshold \
  $nGene_threshold \
  $include_t_cells \
  $analysis_mode' \
  $script_dir/1b.normal_infercnv.R"
#  /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript \
#    --vanilla $script_dir/1b.normal_infercnv.R \
#    $sample_name \
#    $ncores \
#    $nUMI_threshold \
#    $nGene_threshold \
#    $include_t_cells \
#    $analysis_mode