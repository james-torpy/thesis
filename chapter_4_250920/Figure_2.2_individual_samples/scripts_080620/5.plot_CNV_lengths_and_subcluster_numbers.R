#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_3.3_individual_samples"
args = commandArgs(trailingOnly=TRUE)

subtype_labels <- data.frame(
  sample = c("CID3941", "CID4067", "CID4463", "CID4465", "CID4530N", 
    "CID3948", "CID4535", 
    "CID4066", "CID3921", "CID45171", 
    "CID44041", "CID4515", "CID3963", 
    "CID4523"
  ),
  subtype = c(
    rep("lumA", 5),
    rep("lumB", 2),
    rep("HER2", 3),
    rep("TNBC", 3),
    "metaplastic"
  ),
  stringsAsFactors = F
)

sample_names <- c(lumA, lumB, HER2, TNBC, metaplastic)

# samples to manually add subpop numbers for:
manual_cluster_no <- data.frame(
  sample = c(
    "CID4471", "CID4290A",
    "CID4461", 
    "CID3586", 
    "CID4495", "CID44971", "CID44991",
    "CID4513"
  ),
  subtype = c(
    rep("lumA", 2),
    "lumB",
    "HER2",
    rep("TNBC", 3),
    "metaplastic"
  ),
  number = sample(1:8, 8),
  stringsAsFactors = F
)

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
#library(scales, lib.loc = lib_loc)
library(ggplot2)
#library(Seurat)
#library(dplyr)
library(RColorBrewer)
library(cowplot)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
#ref_dir <- paste0(project_dir, "/refs/")
#func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
#raw_dir <- paste0(project_dir, "raw_files/")
#seurat_dir <- paste0(raw_dir, "seurat_objects/", sample_name, "/")
in_path <- paste0(results_dir, "infercnv/")
out_path <- paste0(in_path, "/CNV_length_and_subcluster_distribution/")

Robject_dir <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
table_dir <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", table_dir))
plot_dir <- paste0(out_path, "plots/")
system(paste0("mkdir -p ", plot_dir))

col_palette <- read.table(
  paste0(home_dir, "R/colour_palettes/colour_palette_1.txt"),
  sep = "\t",
  header = F,
  comment.char = ""
)


################################################################################
### 1. Load seurat object and CNV subcluster data ###
################################################################################

for (i in 1:length(sample_names)) {

  epithelial_metadata <- readRDS(
    paste0(
    	in_path, sample_names[i], 
    	"/Rdata/4.subcluster_revised_epithelial_metadata.Rdata"
    )
  )

  CNV_data <- do.call(
  	"rbind",
  	readRDS(
      paste0(
        in_path, sample_names[i], 
       	"/Rdata/CNV_indices_and_lengths.Rdata"
      )
    )
  )
  CNV_data$sample <- sample_names[i]

  if (i==1) {

  	cluster_no <- data.frame(
  	  sample = sample_names[i],
  	  number = length(unique(epithelial_metadata$revised_subcluster_id)),
      stringsAsFactors = F
  	)
  	
  	CNV_lengths <- subset(
  	  CNV_data, select = c(sample, call, length, genomic_length)
  	)

  } else {

  	cluster_no <- rbind(
  	  cluster_no,
  	  data.frame(
    	sample = sample_names[i],
    	number = length(unique(epithelial_metadata$ revised_subcluster_id)),
      stringsAsFactors = F
      )
    )

    CNV_lengths <- rbind(
      CNV_lengths,
      subset(
  	    CNV_data, select = c(sample, call, length, genomic_length)
  	  )
  	)

  }

}

######


cluster_no <- merge(cluster_no, subtype_labels, by = "sample")
#
## manually add cluster numbers of remaining samples:
#cluster_no <- rbind(cluster_no, manual_cluster_no)
## order by subtype:
order_vector <- seq_along(order_vector)
names(order_vector) <- unique(manual_cluster_no$subtype)
cluster_no$order <- order_vector[as.character(cluster_no$subtype)]
cluster_no <- cluster_no[order(cluster_no$order),]

cluster_no$sample <- factor(cluster_no$sample, levels = cluster_no$sample)
cluster_no$subtype <- factor(cluster_no$subtype, levels = unique(cluster_no$subtype))


#################################################################################
#### 2. Plot distributions of CNV lengths across all datasets ###
#################################################################################
#
#######
#CNV_lengths <- CNV_lengths[CNV_lengths$length < 500,]
#######
#
#dplot <- density(CNV_lengths$length)
#png(
#  paste0(plot_dir, "CNV_length_distribution.png"),
#  width = 7, height = 5, units = "in", res = 300
#)
#  plot(
#    dplot,
#    main = "",
#    xlab = "CNV length (number of genes)"
#  )
#dev.off()


################################################################################
#### 3. Plot subcluster number per subtype ###
#################################################################################
cols <- as.character(col_palette[1:length(unique(cluster_no$subtype)),1])

p <- ggplot(cluster_no, aes(x=sample, y=number, fill=subtype))
p <- p + geom_bar(stat = "identity")
p <- p + xlab("Subtype")
p <- p + ylab("No. subpopulations")
p <- p + theme_cowplot(12)
p <- p + scale_fill_manual(values = cols)
#p <- p + scale_fill_manual(values = cols, labels = c("LumA", "LumB", "HER2", "TNBC", "Metaplastic"))
p <- p + labs(fill = "Subtype")
p <- p + theme(
  axis.text.x = element_text(angle = 45, hjust = 1)
)

png(
  paste0(plot_dir, "subcluster_no_per_subtype_partial.png"),
  height = 5,
  width = 7,
  res = 300,
  units = "in"
)
  p
dev.off()

pdf(
  paste0(plot_dir, "subcluster_no_per_subtype_partial.pdf"),
  height = 5,
  width = 7
)
  p
dev.off()







