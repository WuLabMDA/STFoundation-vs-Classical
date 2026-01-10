rm(list = ls(all = TRUE))
graphics.off()

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(FNN)
library(mclust)
library(zellkonverter)  
library(tidyr)
library(SpatialExperiment)


#### Helper functions #######
# Gaussian noise 
add_gaussian_noise <- function(embeddings, noise_sd = 0.1) {
  emb_noisy <- embeddings + matrix(
    rnorm(length(embeddings), mean = 0, sd = noise_sd),
    nrow = nrow(embeddings), ncol = ncol(embeddings)
  )
  return(emb_noisy)
}

#spatially correlated noise 
add_spatial_noise_knn <- function(embeddings, coords, noise_sd = 0.1, k = 20) {
  n <- nrow(embeddings)
  coords <- coords[, sapply(coords, is.numeric), drop = FALSE]  
  
  base_noise <- matrix(rnorm(n * ncol(embeddings), mean = 0, sd = noise_sd),
                       nrow = n, ncol = ncol(embeddings))
  nn <- get.knn(coords, k = k)$nn.index
  
  smoothed_noise <- base_noise
  for (i in 1:n) {
    smoothed_noise[i, ] <- colMeans(base_noise[nn[i, ], , drop = FALSE])
  }
  
  return(embeddings + smoothed_noise)
}

# 
compare_embeddings <- function(ref_emb, noisy_emb, k = 10) {
  cors <- sapply(1:ncol(ref_emb), function(i) cor(ref_emb[, i], noisy_emb[, i]))
  mean_corr <- mean(cors, na.rm = TRUE)
  cl1 <- kmeans(ref_emb, centers = k, nstart = 10)$cluster
  cl2 <- kmeans(noisy_emb, centers = k, nstart = 10)$cluster
  ari <- adjustedRandIndex(cl1, cl2)
  return(list(mean_corr = mean_corr, ari = ari))
}

# 
benchmark_embeddings <- function(ref_emb, coords, noise_levels,
                                 type = c("gaussian", "spatial"), k = 20, n_clusters = 4) {
  type <- match.arg(type)
  results <- data.frame(noise = noise_levels, corr = NA, ari = NA)
  
  for (i in seq_along(noise_levels)) {
    message("Processing iteration ", i, " of ", length(noise_levels),
            " | Noise = ", noise_levels[i])
    
    emb_noisy <- if (type == "gaussian") {
      add_gaussian_noise(ref_emb, noise_sd = noise_levels[i])
    } else {
      add_spatial_noise_knn(ref_emb, coords, noise_sd = noise_levels[i], k = k)
    }
    
    cmp <- compare_embeddings(ref_emb, emb_noisy, k = n_clusters)
    results$corr[i] <- cmp$mean_corr
    results$ari[i] <- cmp$ari
  }
  
  return(results)
}

#
benchmark_one_sample <- function(file_path, method, noise_levels, k = 30) {
  message("Processing sample: ", basename(file_path), " | Method: ", method)
  
  
  if (grepl("\\.h5ad$", file_path)) {
    data <- readH5AD(file_path)
    
    emb <- switch(method,
                  "scGPT"       = data@int_colData@listData$reducedDims$X_scGPT,
                  "Geneformer"  = data@int_colData@listData$reducedDims$X_geneformer,
                  "Nicheformer" = data@int_colData@listData$reducedDims$X_niche_embeddings,
                  stop("Unknown embedding for .h5ad: ", method)
    )
    
    
    sample_name <- tools::file_path_sans_ext(basename(file_path))
    sample_prefix <- substr(sample_name, 1, 5)
    message("  Sample prefix detected: ", sample_prefix)
    
    
    base_raw <- "Y:/Projects/FM_ST/Data/Visium/Raw/"
    visium_dirs <- list.dirs(base_raw, full.names = TRUE, recursive = FALSE)
    
    match_dir <- visium_dirs[grepl(sample_prefix, basename(visium_dirs))]
    data_dir <- file.path(match_dir[1], "outs")
    visium_file <- file.path(data_dir, "filtered_feature_bc_matrix.h5")
    
    
    seu <- Load10X_Spatial(
      data.dir = data_dir,
      filename = "filtered_feature_bc_matrix.h5",
      assay = "Spatial",
      filter.matrix = TRUE,
      to.upper = FALSE
    )
    
    coords <- GetTissueCoordinates(seu) %>% as.data.frame()
    coords <- coords[colnames(data), c("x", "y")]
    
  }
  
  else {
    data <- readRDS(file_path)
    emb <- switch(method,
                  "PCA"    = Embeddings(data, "pca")[, 1:30],
                  "CoCoST" = Embeddings(data, "CoCoST")[, 1:30],
                  "NMF"    = Embeddings(data, "nmf")[, 1:30],
                  stop("Unknown method for .rds: ", method)
    )
    
    coords <- GetTissueCoordinates(data) %>% as.data.frame()
    coords <- coords[, c("x", "y")]
    
  }
  
  results_gaussian <- benchmark_embeddings(emb, coords, noise_levels, type = "gaussian")
  results_spatial  <- benchmark_embeddings(emb, coords, noise_levels, type = "spatial", k = k)
  
  results_gaussian$method <- method
  results_gaussian$type   <- "gaussian"
  results_spatial$method  <- method
  results_spatial$type    <- "spatial"
  
  rbind(results_gaussian, results_spatial)
}

methods <- c("PCA", "CoCoST", "NMF", "scGPT", "Geneformer", "Nicheformer")

base_dirs <- c(
  "Y:/Projects/FM_ST/Data/Visium HD/mouse/PCA Visium embeddings",
  "Y:/Projects/FM_ST/Data/Visium HD/mouse/CoCoST visium embeddings",
  "Y:/Projects/FM_ST/Data/Visium HD/mouse/NMF visium embeddings",
  "X:/maminu/Projects/FM/Final data/Visium_embeddings_scgpt",
  "X:/maminu/Projects/FM/Final data/Visium_embeddings_geneFormer",
  "X:/maminu/Projects/FM/Final data/Visium_embeddings_nicheformer"
)

noise_levels <- round(runif(10, min = 0, max = 1.0), 2)
all_results <- list()

for (m in seq_along(methods)) {
  files <- list.files(base_dirs[m], pattern = "\\.(rds|h5ad)$", full.names = TRUE)
  message("\n==== Benchmarking ", methods[m], " across ", length(files), " samples ====\n")
  
  method_results <- lapply(files, function(f)
    benchmark_one_sample(f, method = methods[m], noise_levels = noise_levels)
  )
  
  all_results[[methods[m]]] <- bind_rows(method_results)
}

df_all <- bind_rows(all_results)

summary_df <- df_all %>%
  group_by(method, type, noise) %>%
  summarise(
    mean_corr = mean(corr, na.rm = TRUE),
    sd_corr   = sd(corr, na.rm = TRUE),
    mean_ari  = mean(ari, na.rm = TRUE),
    sd_ari    = sd(ari, na.rm = TRUE),
    .groups = "drop"
  )

p_corr <- ggplot(summary_df, aes(x = noise, y = mean_corr, color = method, linetype = type)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = mean_corr - sd_corr, ymax = mean_corr + sd_corr, fill = method),
              alpha = 0.15, color = NA) +
  theme_minimal(base_size = 14) +
  labs(title = "Mean correlation vs Noise", x = "Noise SD", y = "Mean PC correlation") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")

p_ari <- ggplot(summary_df, aes(x = noise, y = mean_ari, color = method, linetype = type)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = mean_ari - sd_ari, ymax = mean_ari + sd_ari, fill = method),
              alpha = 0.15, color = NA) +
  theme_minimal(base_size = 14) +
  labs(title = "Mean ARI vs Noise", x = "Noise SD", y = "ARI (k-means)") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")

p_combined <- p_corr + p_ari

print(p_combined)


method_colors <- c(
  "PCA"         = "#08589E", 
  "CoCoST"      = "#E41A1C", 
  "NMF"         = "#800080",  
  "scGPT"       = "#4DAF4A", 
  "Geneformer"  = "#FFD92F", 
  "Nicheformer" = "#FF7F00"  
)

# Correlation plot
p_corr <- ggplot(summary_df, aes(x = noise, y = mean_corr, color = method, linetype = type)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = mean_corr - sd_corr, ymax = mean_corr + sd_corr, fill = method),
              alpha = 0.15, color = NA) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean Correlation vs Noise",
    x = "Noise SD",
    y = "Mean PC Correlation"
  ) +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(face = "bold"),
    panel.grid = element_blank(),        
    panel.background = element_blank(),  
    plot.background = element_blank()    
  )

# ARI plot
p_ari <- ggplot(summary_df, aes(x = noise, y = mean_ari, color = method, linetype = type)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = mean_ari - sd_ari, ymax = mean_ari + sd_ari, fill = method),
              alpha = 0.15, color = NA) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean ARI vs Noise",
    x = "Noise SD",
    y = "ARI (k-means)"
  ) +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(face = "bold"),
    panel.grid = element_blank(),        
    panel.background = element_blank(),  
    plot.background = element_blank()    
  )

p_comb <- p_corr + p_ari + plot_layout(guides = "collect")
p_comb

