#!/bin/bash

project_name="thesis"
subproject_name="Figure_2.5_accuracy_vs_coverage"
ncores=6
sample_name="CID4520N_cancer_sim"
include_t_cells="TRUE"
#simulation_numbers=$(seq 1 30)
simulation_numbers=( "4" )
#downsample_proportions=( "no" "0.9" "0.8" "0.7" "0.6" "0.5" "0.4" "0.3" "0.2" "0.15" "0.1" "0.05" )
downsample_proportions=( "no" )
analysis_mode="samples"
min_CNV_proportion="0.5"
nUMI_threshold="25000"
nGene_threshold="5000"
CNV_type="both"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

for simulation_number in ${simulation_numbers[@]}
  do echo $simulation_number
    for downsample_proportion in ${downsample_proportions[@]}
    	do echo $downsample_proportion
     
        echo "Plotting InferCNV and accuracy results for $sample_name downsampled to $downsample_proportion"
      
        log_dir="$project_dir/logs/$sample_name/$simulation_number/sample_mode/$downsample_proportion.downsampling"
        mkdir -p $log_dir
        echo "Logs are in $log_dir"
        qsub -wd $log_dir -pe smp $ncores -N sim.$simulation_number.$downsample_proportion.accuracy.heatmap -b y -j y -V \
          "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
          CMD BATCH  --no-save \
          '--args \
          $sample_name \
          $include_t_cells \
          $simulation_number \
          $downsample_proportion \
          $analysis_mode \
          $min_CNV_proportion \
          $nUMI_threshold \
          $nGene_threshold \
          $CNV_type' \
          $script_dir/3b.plot_individual_heatmap.R"
#        /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/Rscript --vanilla $script_dir/3b.plot_individual_heatmap.R \
#          $sample_name \
#          $include_t_cells \
#          $simulation_number \
#          $downsample_proportion \
#          $analysis_mode \
#          $neutral_signal_range \
#          $nUMI_threshold \
#          $nGene_threshold \
#          $CNV_type
    done
done

#-P TumourProgression 
