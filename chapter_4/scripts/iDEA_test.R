lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(iDEA, lib.loc=lib_loc)

# load summary data:
data(summary_data)
head(summary_data)

# load annotation data:
data(annotation_data)
head(annotation_data[,1:3])

idea <- CreateiDEAObject(
  summary_data, 
  annotation_data, 
  max_var_beta = 100, 
  min_precent_annot = 0.0025, 
  num_core=10
)

head(idea@summary)

head(idea@annotation[[1]])

idea <- iDEA.fit(
  idea,
  fit_noGS=FALSE,
  init_beta=NULL, 
  init_tau=c(-2,0.5),
  min_degene=5,
  em_iter=15,
  mcmc_iter=1000, 
  fit.tol=1e-5,
  modelVariant = F,
  verbose=TRUE
)

idea <- iDEA.louis(idea)

head(idea@gsea)