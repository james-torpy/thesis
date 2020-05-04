#!/bin/bash

ncores=3

project_name="thesis"
subproject_name="Figure_2.2_accuracy_vs_coverage"
sample_name=CID4520N
subset_data="FALSE"
nUMI_threshold="25000"
nGene_threshold="5000"
CNV_type="both"
CNV_no_range="10_50"
CNV_lengths="20_50_75_100_150_150_200_200_250_250_300_300_350_400_450_500_550_600_650_700_750_800"
CNV_multipliers="3_2_1.5_0.5_0"
downsample="TRUE"
downsample_proportions="0.9_0.8_0.7_0.6_0.5_0.4_0.3_0.2_0.15_0.1_0.05"
#number_of_simulations=1
#downsample_proportions="0.5"
number_of_simulations=30
start_numbering_from=1
end_numbering_at=$(($start_numbering_from+$number_of_simulations-1))
noise_cell_no="5000"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"


for simulation_number in $(seq $start_numbering_from $end_numbering_at)
	do echo $simulation_number
	log_dir="$project_dir/logs/$sample_name/$simulation_number"
	mkdir -p $log_dir
    echo "Logs are in $log_dir"
    qsub -wd $log_dir -pe smp $ncores -N sim.$simulation_number -b y -j y -V \
      "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
      CMD BATCH  --no-save '--args $sample_name $subset_data $nUMI_threshold $nGene_threshold $CNV_type $CNV_no_range $CNV_lengths $CNV_multipliers $downsample $downsample_proportions $simulation_number $noise_cell_no' $script_dir/1b.simulate_cancer.R"
    #  /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla $script_dir/1b.simulate_cancer.R $sample_name $subset_data $nUMI_threshold $nGene_threshold $CNV_no_range $CNV_lengths $CNV_multipliers $downsample $downsample_proportions $simulation_number $noise_cell_no
done

# -P TumourProgression