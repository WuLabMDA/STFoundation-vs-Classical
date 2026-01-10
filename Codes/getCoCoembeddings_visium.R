# delete work space
rm(list = ls(all = TRUE))
graphics.off()


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(kernlab)
library(STutility)
library(CoCoST)


options(future.globals.maxSize = 1024 * 1024 * 10000)

# input/output directories
root_dir   <- "Y:/Projects/FM_ST/Data/Visium/Raw"
output_dir <- "Y:/Projects/FM_ST/Data/Visium HD/mouse/CoCoST embeddings"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

folders <- list.dirs(root_dir, recursive = FALSE)

bfile = folders[5]

bdata_dir <- file.path(bfile, "outs")
bh5_file  <- file.path(bdata_dir, "filtered_feature_bc_matrix.h5")

# Load background data
bseu <- Load10X_Spatial(
  data.dir = bdata_dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  filter.matrix = TRUE,
  to.upper = FALSE
)

bseu <- FindVariableFeatures(bseu)
bseu <- SCTransform(bseu, assay = "Spatial", verbose = FALSE)


# 
timing_results <- data.frame(
  sample = character(),
  time_minutes = numeric(),
  stringsAsFactors = FALSE
)


for (folder in folders) {
  sample_id <- basename(folder)  
  message("Processing sample: ", sample_id)

  data_dir <- file.path(folder, "outs")
  h5_file  <- file.path(data_dir, "filtered_feature_bc_matrix.h5")

  if (!file.exists(h5_file)) {
    stop("No h5 file found for sample: ", sample_id)
  }

  # Load foreground data
  seu <- Load10X_Spatial(
    data.dir = data_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    filter.matrix = TRUE,
    to.upper = FALSE
  )

  seu <- FindVariableFeatures(seu)
  seu <- SCTransform(seu, assay = "Spatial", verbose = FALSE)

  # CoCo-ST
  start_time <- Sys.time()   

  bdata <- bseu@assays[["SCT"]]@scale.data
  fdata <- seu@assays[["SCT"]]@scale.data

  # construct affinity matrices
  rbf <- laplacedot(sigma = 0.10)

  fKernel <- kernelMatrix(rbf, t(fdata))
  Wf <- fKernel@.Data

  bKernel <- kernelMatrix(rbf, t(bdata))
  Wb <- bKernel@.Data


  # Extract contrastive features
  para <- 0.001
  Dim <- 30
  CoCo_model <- CoCoST(t(fdata),Wf,t(bdata),Wb,para,Dim)

  rownames(CoCo_model[["fgComponents"]]) <- colnames(seu@assays[["SCT"]])
  rownames(CoCo_model[["projMatrix"]]) <- rownames(GetAssayData(seu, assay = "SCT", slot = "scale.data"))

  CoCo <- CreateDimReducObject(
    embeddings = CoCo_model[["fgComponents"]],
    loadings = CoCo_model[["projMatrix"]],
    stdev = numeric(),
    key = "CoCoST_",
    assay = "SCT"
  )

  seu@reductions[["CoCoST"]] <- CoCo

  # output filename
  base <- tools::file_path_sans_ext(basename(folder))
  coco_path <- file.path(output_dir, paste0(base, "_coco_embeddings.rds"))

  end_time <- Sys.time()  
  elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))

  # Save timing
  timing_results <- rbind(
    timing_results,
    data.frame(sample = base, time_minutes = elapsed, stringsAsFactors = FALSE)
  )

  cat("Time taken for", base, ":", round(elapsed, 2), "minutes\n")

  # Save results
  saveRDS(seu, coco_path)
  cat("Saved:", coco_path, "\n")

}

# Save timing summary
timing_path <- file.path(output_dir, "coco_embedding_timings.csv")
write.csv(timing_results, timing_path, row.names = FALSE)
cat("Timing summary saved to:", timing_path, "\n")



