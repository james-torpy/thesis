
# define directories:
home_dir <- "/share/ScratchGeneral/jamtor/"
in_dir <- paste0(home_dir, "projects/thesis/lit_data/brca_cna_tgca/")
# fetch coordinates of all brca CNVs, defined by GISTIC:
CNA_coord <- read.table(
  paste0(in_dir, "BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt"),
  header = T
)
# calculate lengths:
CNA_coord$length <- CNA_coord$End-CNA_coord$Start
dim(CNA_coord)
# remove all CNAs < 1 kb, defined by Redon et. al. 2006 as the min length of CNAs:
CNA_coord <- CNA_coord[CNA_coord$length >= 1000,]
nrow(CNA_coord)
# summarise lengths:
summary(CNA_coord$length)

# fetch total sample number:
length(unique(as.character(CNA_coord$Sample)))