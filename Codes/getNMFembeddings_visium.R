# Clear workspace
rm(list = ls(all = TRUE))
graphics.off()

# Load libraries
library(Seurat)
library(dplyr)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(kernlab)
library(NNLM)
library(Matrix)


root_dir   <- "Y:/Projects/FM_ST/Data/Visium/Raw"
output_dir <- "Y:/Projects/FM_ST/Data/Visium HD/mouse/NMF visium embeddings"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

folders <- list.dirs(root_dir, recursive = FALSE)

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

  # Load data
  seu <- Load10X_Spatial(
    data.dir = data_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    filter.matrix = TRUE,
    to.upper = FALSE
  )

  start_time <- Sys.time()  

  DefaultAssay(seu) <- "Spatial"

  # Preprocess
  seu <- FindVariableFeatures(seu)
  seu <- NormalizeData(seu)
  seu <- ScaleData(seu)

  VariableFea <- VariableFeatures(seu)

  mat <- GetAssayData(seu, assay = "Spatial", layer = "counts")[VariableFea, ]

  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0

  set.seed(123)
  nmf_res <- nnmf(mat, k = 50, loss = "mse", method = "scd", max.iter = 200)

  H <- t(nmf_res$H)
  colnames(H) <- paste0("NMF_", seq_len(ncol(H)))

  nmf <- CreateDimReducObject(
    embeddings = H,
    loadings = nmf_res[["W"]],
    stdev = numeric(),
    key = "NMF_",
    assay = "sketch"
  )

  seu@reductions[["nmf"]] <- nmf

  end_time <- Sys.time()  
  elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  seu <- FindNeighbors(seu, assay = "Spatial", reduction = "nmf", dims = 1:50,
                       graph.name = "nmfSNN")
  seu <- FindClusters(seu, cluster.name = "nmf_cluster", resolution = 0.4,
                      graph.name = "nmfSNN")

  # output filename
  base <- tools::file_path_sans_ext(basename(sample_id))
  nmf_path <- file.path(output_dir, paste0(base, "_pca_embeddings.rds"))

  # Store timing
  timing_results <- rbind(
    timing_results,
    data.frame(sample = base, time_minutes = elapsed, stringsAsFactors = FALSE)
  )

  cat("Time taken for", base, ":", round(elapsed, 2), "minutes\n")

  # Save results
  saveRDS(seu, nmf_path)
  cat("Saved:", nmf_path, "\n")

}

# Save timing summary
timing_path <- file.path(output_dir, "NMF_embedding_timings.csv")
write.csv(timing_results, timing_path, row.names = FALSE)
cat("Timing summary saved to:", timing_path, "\n")


