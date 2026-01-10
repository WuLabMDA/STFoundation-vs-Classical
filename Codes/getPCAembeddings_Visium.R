# Clear workspace
rm(list = ls(all = TRUE))
graphics.off()

# Load libraries
library(Seurat)
library(zellkonverter)
library(homologene)
library(dplyr)
library(org.Hs.eg.db)
library(SingleCellExperiment)



root_dir   <- "Y:/Projects/FM_ST/Data/Visium/Raw"
output_dir <- "Y:/Projects/FM_ST/Data/Visium HD/mouse/PCA Visium embeddings"
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
  seu <- RunPCA(seu, assay = "Spatial", reduction.name = "pca")
  
  end_time <- Sys.time()   
  elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  seu <- FindNeighbors(seu, assay = "Spatial", reduction = "pca", dims = 1:50)
  seu <- FindClusters(seu, cluster.name = "pca_clusters", resolution = 0.1)
  
  # output filename
  base <- tools::file_path_sans_ext(basename(sample_id))
  pca_path <- file.path(output_dir, paste0(base, "_pca_embeddings.rds"))
  
  # Store timing
  timing_results <- rbind(
    timing_results,
    data.frame(sample = base, time_minutes = elapsed, stringsAsFactors = FALSE)
  )
  
  cat("Time taken for", base, ":", round(elapsed, 2), "minutes\n")
  
  # Save results
  saveRDS(seu, pca_path)
  cat("Saved:", pca_path, "\n")
  
}

# Save timing summary
timing_path <- file.path(output_dir, "PCA_embedding_timings.csv")
write.csv(timing_results, timing_path, row.names = FALSE)
cat("Timing summary saved to:", timing_path, "\n")


