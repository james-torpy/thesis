#!/bin/bash

ncores=3

project_name="thesis"
subproject_name="Figure_2.3_min_length_and_gap"
sample_name=CID4520N
nUMI_threshold="25000"
nGene_threshold="5000"
gap_or_CNV="CNV"
CNV_type="both"
#feature_lengths="400_300_200_150_125_100_75_50_40_30_20_15_10_5"
feature_lengths="400_150_10"
gap_CNV_length="100"
#number_of_simulations=1
#downsample_proportions="0.5"
number_of_simulations=1
start_numbering_from=1
end_numbering_at=$(($start_numbering_from+$number_of_simulations-1))
noise_cell_no="5000"
t_cells_included="TRUE"
analysis_mode="samples"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"


for simulation_number in $(seq $start_numbering_from $end_numbering_at)
	do echo $simulation_number
	log_dir="$project_dir/logs/simulate_cancer/$gap_or_CNV/$simulation_number/"
	mkdir -p $log_dir
    echo "Logs are in $log_dir"
    qsub -wd $log_dir -pe smp $ncores -N sim.$simulation_number.$gap_or_CNV.$CNV_type -b y -j y -V \
      "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
      CMD BATCH  --no-save '--args \
      $sample_name \
      $nUMI_threshold \
      $nGene_threshold \
      $gap_or_CNV \
      $CNV_type \
      $feature_lengths \
      $gap_CNV_length \
      $simulation_number \
      $noise_cell_no \
      $t_cells_included \
      $analysis_mode' \
       $script_dir/2b.simulate_cancer.R"
    #  /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/2b.simulate_cancer.R \
#        $sample_name \
#        $subset_data \
#        $nUMI_threshold \
#        $nGene_threshold \
#        $gap_or_CNV \
#        CNV_type \
#        $feature_lengths \
#        $gap_CNV_length \
#        $simulation_number \
#        $noise_cell_no \
#        $t_cells_included \
#		 $analysis_mode
done

# -P TumourProgression