#!/bin/bash

ncores=10

project_name="thesis"
subproject_name="Figure_2.2_accuracy_vs_coverage"
sample_name=CID4520N
subset_data="FALSE"
nUMI_threshold="25000"
nGene_threshold="5000"
CNV_no_range="10_40"
CNV_lengths="50_75_100_150_200_300_400"
CNV_multipliers="3_2_1.5_0.5_0"
downsample="TRUE"
#downsample_proportions="0.9_0.8_0.7_0.6_0.5_0.4_0.3_0.2_0.15_0.1_0.05"
#number_of_simulations=30
downsample_proportions="0.9"
number_of_simulations=1

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"
log_dir="$project_dir/logs/$sample_name"
mkdir -p $log_dir


for simulation_number in $(seq 1 $number_of_simulations)
	do echo $simulation_number
    echo "Logs are in $log_dir"
    qsub -wd $log_dir -pe smp $ncores -N sim.$simulation_number -b y -j y -V -P TumourProgression \
      "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
      CMD BATCH  --no-save '--args $sample_name $subset_data $nUMI_threshold $nGene_threshold $CNV_no_range $CNV_lengths $CNV_multipliers $downsample $downsample_proportions $simulation_number' $script_dir/1b.simulate_cancer.R"
    #  /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/1b.simulate_cancer.R $sample_name $subset_data $nUMI_threshold $nGene_threshold $CNV_no_range $CNV_lengths $CNV_multipliers $downsample $downsample_proportions $simulation_number
done
