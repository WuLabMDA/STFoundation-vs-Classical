rm(list = ls(all = TRUE))
graphics.off()

library(Seurat)
library(Matrix)
library(zellkonverter)  
library(dplyr)
library(ggplot2)
library(SpatialExperiment)
library(openxlsx)


### helper functions #######
compute_sparsity <- function(expr_mat) {
    total_entries <- length(expr_mat)
    sparsity <- sum(expr_mat == 0) / total_entries
  return(sparsity)
}

process_cocost_sample <- function(file_path) {
  message("Processing: ", basename(file_path))

  if (grepl("\\.rds$", file_path)) {
    seu <- readRDS(file_path)
    expr_mat <- GetAssayData(seu, slot = "counts")
  } else if (grepl("\\.h5ad$", file_path)) {
    adata <- readH5AD(file_path)
    expr_mat <- assay(adata, "X")  
  } else {
    stop("Unsupported file type: ", file_path)
  }

  sparsity <- compute_sparsity(expr_mat)

  n_cells <- ncol(expr_mat)
  n_genes <- nrow(expr_mat)

  data.frame(
    sample = tools::file_path_sans_ext(basename(file_path)),
    n_cells = n_cells,
    n_genes = n_genes,
    sparsity = sparsity
  )
}

cocost_dir <- "Y:/Projects/FM_ST/Data/Visium HD/mouse/CoCoST visium embeddings"
files <- list.files(cocost_dir, pattern = "\\.(rds|h5ad)$", full.names = TRUE)


results_list <- lapply(files, process_cocost_sample)
sparsity_df <- bind_rows(results_list)

summary_df <- sparsity_df %>%
  summarise(
    mean_sparsity = mean(sparsity),
    sd_sparsity = sd(sparsity),
    mean_cells = mean(n_cells),
    mean_genes = mean(n_genes)
  )


ggplot(sparsity_df, aes(x = reorder(sample, sparsity), y = sparsity)) +
  geom_bar(stat = "identity", fill = "#E41A1C") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Sparsity of CoCoST Datasets",
    x = "Sample",
    y = "Sparsity (fraction of zeros)"
  ) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

#write.xlsx(sparsity_df, "Y:/Projects/FM_ST/Results/sparsity experiment/sparsity_summary_visium.xlsx", rowNames = FALSE)
