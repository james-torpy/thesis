#!/bin/bash

project_name="thesis"
subproject_name="Figure_2.6_random_permutation_of_normal"

sample_name=permutated_CID4520N
permutation_proportion="0.005"
simulation_number="1"
t_cells_included="TRUE"
analysis_mode="samples"
neutral_signal_range="0.97_1.03"
numcores=6
remove_genes_number="5"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"
     
echo "Removing artefacts (repeat no. $remove_genes_number) and running InferCNV for simulation \
  $simulation_number $sample_name with $permutation_proportion genes permutated..."

log_dir="$project_dir/logs/remove_artefact_genes/$permutation_proportion/$simulation_number/$remove_genes_number/$analysis_mode.mode/"
mkdir -p $log_dir
echo "Logs are in $log_dir"
qsub -wd $log_dir -pe smp $numcores -N infercnv.$CNV_type.$gap_or_CNV.$sample_name -b y -j y -V -P TumourProgression \
  "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
  CMD BATCH  --no-save '--args \
  $sample_name \
  $permutation_proportion \
  $simulation_number \
  $t_cells_included \
  $analysis_mode \
  $neutral_signal_range \
  $numcores \
  $remove_genes_number' \
  $script_dir/4b.remove_problem_genes_infercnv.R"

