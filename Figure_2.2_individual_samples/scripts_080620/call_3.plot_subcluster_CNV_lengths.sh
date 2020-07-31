#!/bin/bash

project_name="thesis"
subproject_name="Figure_3.3_individual_samples"
ncores=2

#sample_names=( "CID3586" "CID4066" "CID4471" "CID4515" \
#  "CID3921" "CID4067" "CID4495" "CID45171" \
#  "CID3941" "CID4290A" "CID4461" "CID44971" \
#  "CID4523" "CID3948" "CID4463" "CID44991" \
#  "CID4530N" "CID3963" "CID44041" "CID4465" \
#  "CID4513" "CID4535" )
sample_names=( "CID4465" )

min_CNV_length=20
min_CNV_proportion="0.5"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

for sample_name in ${sample_names[@]}
  do echo $sample_name

  log_dir="$project_dir/logs/signal_plot/$sample_name/"
  mkdir -p $log_dir
  echo "Logs are in $log_dir"
  qsub -wd $log_dir -pe smp $ncores -N plot.$sample_name -b y -j y -V -P TumourProgression \
    "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
    CMD BATCH  --no-save \
    '--args \
    $sample_name \
    $min_CNV_length \
    $min_CNV_proportion' \
    $script_dir/3.plot_subcluster_CNV_lengths.R"

done

