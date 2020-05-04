#!/bin/bash

ncores=3

project_name="thesis"
subproject_name="Figure_2.6_random_permutation_of_normal"
sample_name=CID4520N
nUMI_threshold="25000"
nGene_threshold="5000"
permutation_proportions=( "0.01" "0.05" "0.1" "0.2" "0.3" "0.4" )
potential_values="0_0.5_1.5_2_3"
#number_of_simulations=1
#downsample_proportions="0.5"
number_of_simulations=1
start_numbering_from=1
end_numbering_at=$(($start_numbering_from+$number_of_simulations-1))
t_cells_included="TRUE"
analysis_mode="TRUE"
neutral_signal_range="0.97_1.03"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

for permutation_proportion in permutation_proportions
  do echo "permutation proportion is: $permutation_proportion"
  for simulation_number in $(seq $start_numbering_from $end_numbering_at)
    do echo "simulation number is: $simulation_number"

    log_dir="$project_dir/logs/permutate_normal/$permutation_proportion.permutation_proprortion.$simulation_number/$simulation_number/"
    mkdir -p $log_dir
    echo "Logs are in $log_dir"
    qsub -wd $log_dir -pe smp $ncores -N sim.$simulation_number.$gap_or_CNV.$CNV_type -b y -j y -V -P TumourProgression \
      "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
      CMD BATCH  --no-save '--args \
      $sample_name \
      $nUMI_threshold \
      $nGene_threshold \
      $permutation_proportion \
      $potential_values \
      $simulation_number' \
      $t_cells_included \
      $analysis_mode \
      $neutral_signal_range \
       $script_dir/1b.randomly_permutate_normal.R"
#    /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/1b.randomly_permutate_normal.R \
#      $sample_name \
#      $nUMI_threshold \
#      $nGene_threshold \
#      $permutation_proportion \
#      $potential_values \
#      $simulation_number \
#      $t_cells_included \
#      $analysis_mode \
#      $neutral_signal_range \
  done
done

# 