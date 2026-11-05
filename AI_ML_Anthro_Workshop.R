#!/usr/bin/env Rscript

# AI and machine learning for anthropology: 90-minute workshop demo
#
# Sections:
#   1. Iris data and features
#   2. Unsupervised learning with PCA and UMAP
#   3. Supervised learning with a simple train/test k-nearest-neighbor model
#   4. Read and visualize one Avizo ASCII PLY tooth mesh
#   5. Create interpretable 3-D tooth features
#
# Required packages:
#   install.packages(c("ggplot2", "uwot"))
#
# In RStudio, set the working directory to the folder containing this script
# and the upload folder. Then run one section at a time with Cmd/Ctrl+Enter.

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 first: install.packages('ggplot2')")
}
if (!requireNamespace("uwot", quietly = TRUE)) {
  stop("Install uwot first: install.packages('uwot')")
}
library(ggplot2)

#set.seed(2026)
#output_dir <- "workshop_outputs"
#ply_file <- file.path("upload", "Lagurodon 2959.ply")
#dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

##### Change for your own system if you want to save outputs elsewhere. #####
setwd("~/Documents/Grants/anthro_data/presentation")
#############################################################################

# -----------------------------------------------------------------------------
# 1. IRIS DATA: ROWS ARE OBSERVATIONS; COLUMNS ARE FEATURES
# -----------------------------------------------------------------------------

iris_data <- iris
cat("\nFirst six Iris observations:\n")
print(head(iris_data))
cat("\nSpecies counts:\n")
print(table(iris_data$Species))

iris_features <- as.matrix(iris_data[, 1:4])
iris_labels <- iris_data$Species

# Standardization places measurements with different units/ranges on a
# comparable scale. The species labels are deliberately excluded.
iris_scaled <- scale(iris_features)

# -----------------------------------------------------------------------------
# 2. UNSUPERVISED LEARNING: PCA AND UMAP DO NOT USE SPECIES DURING FITTING
# -----------------------------------------------------------------------------

iris_pca <- prcomp(iris_scaled, center = FALSE, scale. = FALSE)
iris_pca_scores <- data.frame(
  PC1 = iris_pca$x[, 1],
  PC2 = iris_pca$x[, 2],
  Species = iris_labels
)

variance_percent <- 100 * iris_pca$sdev^2 / sum(iris_pca$sdev^2)

p_iris_pca <- ggplot(iris_pca_scores, aes(PC1, PC2, color = Species)) +
  geom_point(size = 2.6, alpha = 0.85) +
  labs(
    title = "PCA of Iris measurements",
    subtitle = "Species labels color the plot after PCA; PCA did not use them",
    x = sprintf("PC1 (%.1f%% of variance)", variance_percent[1]),
    y = sprintf("PC2 (%.1f%% of variance)", variance_percent[2])
  ) +
  theme_minimal(base_size = 14)

show(p_iris_pca)

#ggsave(file.path(output_dir, "iris_pca.png"), p_iris_pca,
#       width = 8, height = 6, dpi = 220)

# UMAP preserves local neighborhoods. Its axes have no direct biological units.
# Setting a seed and n_threads = 1 makes this classroom example reproducible.
set.seed(2026)
iris_umap_matrix <- uwot::umap(
  iris_scaled,
  n_neighbors = 15,
  min_dist = 0.15,
  metric = "euclidean",
  n_threads = 1,
  verbose = FALSE
)

iris_umap <- data.frame(
  UMAP1 = iris_umap_matrix[, 1],
  UMAP2 = iris_umap_matrix[, 2],
  Species = iris_labels
)

p_iris_umap_unlabeled <- ggplot(iris_umap, aes(UMAP1, UMAP2)) +
  geom_point(size = 2.6, alpha = 0.80, color = "#31688E") +
  labs(
    title = "UMAP finds neighborhoods without species labels",
    subtitle = "What groups do you see?"
  ) +
  theme_minimal(base_size = 14)

show(p_iris_umap_unlabeled)

p_iris_umap_labeled <- ggplot(iris_umap, aes(UMAP1, UMAP2, color = Species)) +
  geom_point(size = 2.6, alpha = 0.85) +
  labs(
    title = "The same UMAP with species revealed",
    subtitle = "Color adds known labels after the unsupervised embedding"
  ) +
  theme_minimal(base_size = 14)

show(p_iris_umap_labeled)

#ggsave(file.path(output_dir, "iris_umap_unlabeled.png"),
#       p_iris_umap_unlabeled, width = 8, height = 6, dpi = 220)
#ggsave(file.path(output_dir, "iris_umap_labeled.png"),
#       p_iris_umap_labeled, width = 8, height = 6, dpi = 220)

# -----------------------------------------------------------------------------
# 3. SUPERVISED LEARNING: LABELS GUIDE A K-NEAREST-NEIGHBOR CLASSIFIER
# -----------------------------------------------------------------------------

# Split observations into training and test sets within each species.
set.seed(2026)
train_index <- unlist(lapply(split(seq_len(nrow(iris_data)), iris_labels),
  function(index) sample(index, size = round(0.70 * length(index)))))
test_index <- setdiff(seq_len(nrow(iris_data)), train_index)

train_x <- iris_scaled[train_index, , drop = FALSE]
test_x <- iris_scaled[test_index, , drop = FALSE]
train_y <- iris_labels[train_index]
test_y <- iris_labels[test_index]

# A transparent kNN implementation avoids an additional package dependency.
predict_knn <- function(train_x, train_y, test_x, k = 5) {
  predictions <- character(nrow(test_x))
  for (i in seq_len(nrow(test_x))) {
    distances <- sqrt(rowSums((train_x -
      matrix(test_x[i, ], nrow(train_x), ncol(train_x), byrow = TRUE))^2))
    neighbors <- order(distances)[seq_len(k)]
    votes <- sort(table(train_y[neighbors]), decreasing = TRUE)
    predictions[i] <- names(votes)[1]
  }
  factor(predictions, levels = levels(train_y))
}

predicted_species <- predict_knn(train_x, train_y, test_x, k = 5)
confusion_matrix <- table(Actual = test_y, Predicted = predicted_species)
accuracy <- mean(predicted_species == test_y)

cat("\nIris kNN confusion matrix:\n")
print(confusion_matrix)
cat(sprintf("Test accuracy for this one split: %.1f%%\n", 100 * accuracy))
cat("Accuracy from one small split is an illustration, not a final estimate.\n")

#write.csv(as.data.frame.matrix(confusion_matrix),
#          file.path(output_dir, "iris_knn_confusion_matrix.csv"))

# -----------------------------------------------------------------------------
# 4. READ ONE AVIZO ASCII PLY FILE
# -----------------------------------------------------------------------------

read_avizo_ply <- function(path) {
  if (!file.exists(path)) stop("PLY file not found: ", path)
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)

  header <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) stop("Missing end_header in ", path)
    header <- c(header, line)
    if (trimws(line) == "end_header") break
  }

  element_count <- function(name) {
    match <- grep(paste0("^element ", name, " "), header, value = TRUE)
    if (length(match) != 1) stop("Could not find element ", name)
    as.integer(strsplit(match, "\\s+")[[1]][3])
  }

  n_vertices <- element_count("vertex")
  n_faces <- element_count("face")

  vertices <- scan(con, what = numeric(), nmax = 3 * n_vertices,
                   quiet = TRUE)
  vertices <- matrix(vertices, ncol = 3, byrow = TRUE,
                     dimnames = list(NULL, c("x", "y", "z")))

  # Each Avizo face row in these files contains:
  # 3, vertex_index_1, vertex_index_2, vertex_index_3, patch
  face_values <- scan(con, what = integer(), nmax = 5 * n_faces,
                      quiet = TRUE)
  if (length(face_values) != 5 * n_faces) {
    stop("Unexpected face structure in ", path)
  }
  face_rows <- matrix(face_values, ncol = 5, byrow = TRUE)
  if (any(face_rows[, 1] != 3L)) stop("This demo expects triangular faces")

  # PLY indices begin at zero; R indices begin at one.
  faces <- face_rows[, 2:4, drop = FALSE] + 1L

  list(
    path = path,
    header = header,
    vertices = vertices,
    faces = faces
  )
}

orient_vertices_by_pca <- function(vertices) {
  fit <- prcomp(vertices, center = TRUE, scale. = FALSE)
  xyz <- fit$x

  # PCA directions have arbitrary signs. Fix signs for reproducible views.
  for (j in 1:3) {
    extreme <- which.max(abs(xyz[, j]))
    if (xyz[extreme, j] < 0) xyz[, j] <- -xyz[, j]
  }
  list(xyz = xyz, sdev = fit$sdev, rotation = fit$rotation)
}

cross_product_rows <- function(a, b) {
  cbind(
    a[, 2] * b[, 3] - a[, 3] * b[, 2],
    a[, 3] * b[, 1] - a[, 1] * b[, 3],
    a[, 1] * b[, 2] - a[, 2] * b[, 1]
  )
}

ply_file <- file.path("tooth_data", "Lagurodon_2959.ply")
tooth_mesh <- read_avizo_ply(ply_file)
# Use Principal Component Analysis to orient the tooth mesh vertices in a reproducible way.
tooth_orientation <- orient_vertices_by_pca(tooth_mesh$vertices)
tooth_xyz <- tooth_orientation$xyz

cat("\nPLY header summary:\n")
cat("Vertices:", nrow(tooth_mesh$vertices), "\n")
cat("Triangular faces:", nrow(tooth_mesh$faces), "\n")

# Downsample only for plotting. All vertices and faces remain available for
# feature calculations.
keep <- unique(round(seq(1, nrow(tooth_xyz),
                         length.out = min(12000, nrow(tooth_xyz)))))
plot_xyz <- tooth_xyz[keep, , drop = FALSE]

tooth_views <- rbind(
  data.frame(view = "PC1-PC2", horizontal = plot_xyz[, 1],
             vertical = plot_xyz[, 2], depth = plot_xyz[, 3]),
  data.frame(view = "PC1-PC3", horizontal = plot_xyz[, 1],
             vertical = plot_xyz[, 3], depth = plot_xyz[, 2]),
  data.frame(view = "PC2-PC3", horizontal = plot_xyz[, 2],
             vertical = plot_xyz[, 3], depth = plot_xyz[, 1])
)

p_tooth <- ggplot(tooth_views,
                  aes(horizontal, vertical, color = depth)) +
  geom_point(size = 0.18, alpha = 0.80) +
  facet_wrap(~view) +
  coord_equal() +
  scale_color_viridis_c(option = "C", guide = "none") +
  labs(
    title = paste("Three PCA-aligned views of", basename(ply_file)),
    subtitle = "Each point is a surface vertex; color represents depth",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text = element_blank(), panel.grid = element_blank())

show(p_tooth)

#ggsave(file.path(output_dir, "one_tooth_three_views.png"), p_tooth,
#       width = 12, height = 4.5, dpi = 220)

# -----------------------------------------------------------------------------
# 5. CREATE INTERPRETABLE TOOTH FEATURES
# -----------------------------------------------------------------------------

extents <- sort(apply(tooth_xyz, 2, function(z) diff(range(z))),
                decreasing = TRUE)

triangle_1 <- tooth_mesh$vertices[tooth_mesh$faces[, 1], , drop = FALSE]
triangle_2 <- tooth_mesh$vertices[tooth_mesh$faces[, 2], , drop = FALSE]
triangle_3 <- tooth_mesh$vertices[tooth_mesh$faces[, 3], , drop = FALSE]

triangle_cross <- cross_product_rows(triangle_2 - triangle_1,
                                     triangle_3 - triangle_1)
surface_area <- sum(sqrt(rowSums(triangle_cross^2))) / 2
signed_volume <- sum(rowSums(
  triangle_1 * cross_product_rows(triangle_2, triangle_3)
)) / 6
enclosed_volume <- abs(signed_volume)

radius <- sqrt(rowSums(tooth_xyz^2))
radius_quantiles <- unname(quantile(radius, c(0.10, 0.50, 0.90)))
largest_extent <- extents[1]

sphericity <- if (enclosed_volume > 0 && surface_area > 0) {
  pi^(1 / 3) * (6 * enclosed_volume)^(2 / 3) / surface_area
} else {
  NA_real_
}

tooth_features <- data.frame(
  specimen = tools::file_path_sans_ext(basename(ply_file)),
  vertices = nrow(tooth_mesh$vertices),
  triangular_faces = nrow(tooth_mesh$faces),
  length = extents[1],
  width = extents[2],
  height = extents[3],
  elongation = extents[1] / extents[2],
  flatness = extents[2] / extents[3],
  surface_area = surface_area,
  enclosed_volume = enclosed_volume,
  surface_to_volume = surface_area / enclosed_volume,
  sphericity = sphericity,
  radial_q10_scaled = radius_quantiles[1] / largest_extent,
  radial_q50_scaled = radius_quantiles[2] / largest_extent,
  radial_q90_scaled = radius_quantiles[3] / largest_extent,
  pc_variance_ratio_1_2 =
    (tooth_orientation$sdev[1] / tooth_orientation$sdev[2])^2,
  pc_variance_ratio_2_3 =
    (tooth_orientation$sdev[2] / tooth_orientation$sdev[3])^2
)

cat("\nInterpretable tooth features:\n")
print(tooth_features, row.names = FALSE)

#write.csv(tooth_features,
#          file.path(output_dir, "one_tooth_features.csv"),
#          row.names = FALSE)

#cat("\nWorkshop outputs written to:", normalizePath(output_dir), "\n")
cat("Scientific caution: similarity in selected features does not establish species identity.\n")

