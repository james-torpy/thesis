#!/bin/bash

project_name="thesis"
subproject_name="Figure_3.3_individual_samples"
ncores=2

sample_names=( "CID3586" "CID4066" "CID4471" "CID4515" \
  "CID3921" "CID4067" "CID4495" "CID45171" \
  "CID3941" "CID4290A" "CID4461" "CID44971" \
  "CID4523" "CID3948" "CID4463" "CID44991" \
  "CID4530N" "CID3963" "CID44041" "CID4465" \
  "CID4513" "CID4535" )
all_cell_PC="C"
all_cell_res="PC_C_res.0.6"
malignant_PC="D"
malignant_res="SUBSET_D_res.0.8"
broad_markers="ACTB_PTPRC_CD19_CD3D_CD68_PDGFRB_PECAM1_EPCAM"
epithelial_markers="ACTB_EPCAM_KRT18_KRT5_MKI67"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$subproject_name"
script_dir="$project_dir/scripts/"

for sample_name in ${sample_names[@]}
  do echo $sample_name

  log_dir="$project_dir/logs/plot/$sample_name/"
  mkdir -p $log_dir
  echo "Logs are in $log_dir"
  qsub -wd $log_dir -pe smp $ncores -N plot.$sample_name -b y -j y -V -P TumourProgression \
    "/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R \
    CMD BATCH  --no-save \
    '--args \
    $sample_name \
    $all_cell_PC \
    $all_cell_res \
    $malignant_PC \
    $malignant_res \
    $broad_markers \
    $epithelial_markers' \
    $script_dir/2.plot_individual_heatmap.R"
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

