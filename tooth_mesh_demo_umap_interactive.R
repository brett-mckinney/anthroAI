#!/usr/bin/env Rscript

# Exploratory 3-D tooth morphology from ASCII Avizo PLY meshes.
# Teaching demonstration only: distance to one known specimen is not a
# taxonomic identification or a validated classifier.
#

setwd("~/Documents/Grants/anthro_data/tooth_mesh_unsupervised_R/tooth_mesh_R_ggplot_demo")
input_dir <- "tooth_data"
output_dir <- "tooth_results_R_umap"
known_pattern <- "Lagurodon"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 first: install.packages('ggplot2')")
}
if (!requireNamespace("uwot", quietly = TRUE)) {
  stop("Install uwot first: install.packages('uwot')")
}
library(ggplot2)

cross3 <- function(a, b) {
  cbind(
    a[, 2] * b[, 3] - a[, 3] * b[, 2],
    a[, 3] * b[, 1] - a[, 1] * b[, 3],
    a[, 1] * b[, 2] - a[, 2] * b[, 1]
  )
}

read_avizo_ply <- function(path) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)
  header <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) stop("Missing end_header in ", path)
    header <- c(header, line)
    if (trimws(line) == "end_header") break
  }

  get_count <- function(element) {
    hit <- grep(paste0("^element ", element, " "), header, value = TRUE)
    if (length(hit) != 1) stop("Could not find ", element, " count in ", path)
    as.integer(strsplit(hit, "\\s+")[[1]][3])
  }
  n_vertices <- get_count("vertex")
  n_faces <- get_count("face")

  # The connection now points to the first vertex.
  vertices <- scan(con, what = numeric(), nmax = 3 * n_vertices,
                   quiet = TRUE)
  vertices <- matrix(vertices, ncol = 3, byrow = TRUE,
                     dimnames = list(NULL, c("x", "y", "z")))

  # Each Avizo face row is: 3, vertex_i, vertex_j, vertex_k, patch.
  # scan() is substantially faster than parsing these dense rows one at a time.
  face_values <- scan(con, what = integer(), nmax = 5 * n_faces,
                      quiet = TRUE)
  if (length(face_values) != 5 * n_faces) {
    stop("Unexpected face structure in ", path)
  }
  face_rows <- matrix(face_values, ncol = 5, byrow = TRUE)
  if (any(face_rows[, 1] != 3L)) stop("Non-triangular face in ", path)
  faces <- face_rows[, 2:4, drop = FALSE] + 1L # PLY uses zero-based indices.
  list(vertices = vertices, faces = faces)
}

orient_by_pca <- function(vertices) {
  fit <- prcomp(vertices, center = TRUE, scale. = FALSE)
  xyz <- fit$x
  # Fix arbitrary PCA signs for reproducible plots.
  for (j in 1:3) {
    extreme <- which.max(abs(xyz[, j]))
    if (xyz[extreme, j] < 0) xyz[, j] <- -xyz[, j]
  }
  list(xyz = xyz, sdev = fit$sdev)
}

measure_mesh <- function(mesh) {
  oriented <- orient_by_pca(mesh$vertices)
  xyz <- oriented$xyz
  extents <- sort(apply(xyz, 2, function(z) diff(range(z))),
                  decreasing = TRUE)

  tri1 <- mesh$vertices[mesh$faces[, 1], , drop = FALSE]
  tri2 <- mesh$vertices[mesh$faces[, 2], , drop = FALSE]
  tri3 <- mesh$vertices[mesh$faces[, 3], , drop = FALSE]
  cross <- cross3(tri2 - tri1, tri3 - tri1)
  surface_area <- sum(sqrt(rowSums(cross^2))) / 2
  signed_volume <- sum(rowSums(tri1 * cross3(tri2, tri3))) / 6
  volume <- abs(signed_volume)

  radius <- sqrt(rowSums(xyz^2))
  rq <- unname(quantile(radius, c(.10, .50, .90)))
  scale <- max(extents[1], .Machine$double.eps)
  sphericity <- if (volume > 0 && surface_area > 0) {
    pi^(1 / 3) * (6 * volume)^(2 / 3) / surface_area
  } else {
    NA_real_
  }

  features <- data.frame(
    vertices = nrow(mesh$vertices),
    faces = nrow(mesh$faces),
    length = extents[1],
    width = extents[2],
    height = extents[3],
    elongation = extents[1] / extents[2],
    flatness = extents[2] / extents[3],
    surface_area = surface_area,
    enclosed_volume = volume,
    surface_to_volume = surface_area / volume,
    sphericity = sphericity,
    radial_q10_scaled = rq[1] / scale,
    radial_q50_scaled = rq[2] / scale,
    radial_q90_scaled = rq[3] / scale,
    pc_variance_ratio_1_2 = (oriented$sdev[1] / oriented$sdev[2])^2,
    pc_variance_ratio_2_3 = (oriented$sdev[2] / oriented$sdev[3])^2
  )
  list(features = features, xyz = xyz)
}

robust_scale <- function(x) {
  med <- apply(x, 2, median, na.rm = TRUE)
  spread <- apply(x, 2, mad, constant = 1.4826, na.rm = TRUE)
  spread[!is.finite(spread) | spread == 0] <- 1
  sweep(sweep(x, 2, med, "-"), 2, spread, "/")
}

paths <- list.files(input_dir, pattern = "\\.ply$", full.names = TRUE,
                    ignore.case = TRUE)
if (length(paths) < 2) stop("Need at least two PLY files in ", input_dir)
known <- grepl(known_pattern, basename(paths), ignore.case = TRUE)
if (sum(known) != 1) stop("Expected exactly one known specimen filename")

feature_rows <- list()
view_rows <- list()
set.seed(2026)

for (i in seq_along(paths)) {
  specimen <- tools::file_path_sans_ext(basename(paths[i]))
  message("Reading ", basename(paths[i]), " ...")
  measured <- measure_mesh(read_avizo_ply(paths[i]))
  feature_rows[[i]] <- cbind(
    data.frame(specimen = specimen, known = known[i]),
    measured$features
  )

  xyz <- measured$xyz
  keep <- unique(round(seq(1, nrow(xyz), length.out = min(8000, nrow(xyz)))))
  xyz <- xyz[keep, , drop = FALSE]
  view_rows[[i]] <- rbind(
    data.frame(specimen, view = "PC1-PC2", horizontal = xyz[, 1],
               vertical = xyz[, 2], depth = xyz[, 3]),
    data.frame(specimen, view = "PC1-PC3", horizontal = xyz[, 1],
               vertical = xyz[, 3], depth = xyz[, 2]),
    data.frame(specimen, view = "PC2-PC3", horizontal = xyz[, 2],
               vertical = xyz[, 3], depth = xyz[, 1])
  )
}

features <- do.call(rbind, feature_rows)
views <- do.call(rbind, view_rows)
write.csv(features, file.path(output_dir, "tooth_mesh_features_R.csv"),
          row.names = FALSE)

shape_names <- c(
  "elongation", "flatness", "sphericity", "radial_q10_scaled",
  "radial_q50_scaled", "radial_q90_scaled", "pc_variance_ratio_1_2",
  "pc_variance_ratio_2_3"
)
shape <- as.matrix(features[, shape_names])
for (j in seq_len(ncol(shape))) {
  shape[!is.finite(shape[, j]), j] <- median(shape[, j], na.rm = TRUE)
}
scaled <- robust_scale(shape)

pca <- prcomp(scaled, center = FALSE, scale. = FALSE)
pca_scores <- data.frame(
  specimen = features$specimen,
  known = features$known,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)
write.csv(pca_scores, file.path(output_dir, "tooth_shape_pca_scores_R.csv"),
          row.names = FALSE)

# UMAP emphasizes local neighborhoods rather than total variance. With only ten
# specimens this is an exploratory visualization, not evidence of discrete
# species clusters. A fixed seed makes the classroom result reproducible.
set.seed(2026)
umap_neighbors <- min(5L, nrow(scaled) - 1L)
umap_embedding <- uwot::umap(
  scaled,
  n_neighbors = umap_neighbors,
  n_components = 2,
  min_dist = 0.20,
  metric = "euclidean",
  scale = FALSE,
  init = "spectral",
  n_threads = 1,
  ret_model = FALSE,
  verbose = FALSE
)
umap_scores <- data.frame(
  specimen = features$specimen,
  known = features$known,
  UMAP1 = umap_embedding[, 1],
  UMAP2 = umap_embedding[, 2]
)
write.csv(umap_scores, file.path(output_dir, "tooth_shape_umap_scores_R.csv"),
          row.names = FALSE)

known_index <- which(features$known)
distance <- sqrt(rowSums((scaled -
  matrix(scaled[known_index, ], nrow(scaled), ncol(scaled), byrow = TRUE))^2))
similarity <- data.frame(
  specimen = features$specimen,
  distance_to_known_lagurodon = distance,
  reference = features$known
)
similarity <- similarity[order(similarity$distance_to_known_lagurodon), ]
similarity$similarity_rank <- 0:(nrow(similarity) - 1)
write.csv(similarity,
          file.path(output_dir, "lagurodon_morphology_similarity_R.csv"),
          row.names = FALSE)

known_name <- features$specimen[known_index]

p_views <- ggplot(views, aes(horizontal, vertical, color = depth)) +
  geom_point(size = 0.12, alpha = 0.8) +
  coord_equal() +
  facet_grid(specimen ~ view) +
  #facet_grid(specimen ~ view, scales = "free") +
  scale_color_viridis_c(option = "C", guide = "none") +
  labs(
    title = "PCA-aligned tooth surfaces",
    subtitle = "Orientation standardized; color represents depth",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(), panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0)
  )
ggsave(file.path(output_dir, "tooth_multiview_R.png"), p_views,
       width = 10, height = 18, dpi = 180)

print(p_views)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
p_pca <- ggplot(pca_scores, aes(PC1, PC2, label = specimen)) +
  geom_hline(yintercept = 0, color = "grey85") +
  geom_vline(xintercept = 0, color = "grey85") +
  geom_point(aes(shape = known, fill = known), size = 4, color = "black") +
  geom_text(nudge_x = 0.08, nudge_y = 0.08, check_overlap = TRUE) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23)) +
  scale_fill_manual(values = c(`FALSE` = "#2878B5", `TRUE` = "gold")) +
  labs(
    title = "Unsupervised morphology map—not a species classifier",
    x = sprintf("PC1 (%.1f%% variance)", variance[1]),
    y = sprintf("PC2 (%.1f%% variance)", variance[2]),
    shape = "Known Lagurodon", fill = "Known Lagurodon"
  ) +
  theme_minimal(base_size = 13)
ggsave(file.path(output_dir, "tooth_shape_pca_R.png"), p_pca,
       width = 9, height = 7, dpi = 220)

print(p_pca)

p_umap <- ggplot(umap_scores, aes(UMAP1, UMAP2, label = specimen)) +
  geom_point(aes(shape = known, fill = known), size = 4, color = "black") +
  geom_text(nudge_x = 0.08, nudge_y = 0.08, check_overlap = TRUE) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23)) +
  scale_fill_manual(values = c(`FALSE` = "#2878B5", `TRUE` = "gold")) +
  labs(
    title = "UMAP of scale-free tooth morphology",
    subtitle = paste0(
      "Exploratory only: n = ", nrow(umap_scores),
      ", neighbors = ", umap_neighbors,
      ", min_dist = 0.20; distances between clusters are not quantitative"
    ),
    x = "UMAP 1", y = "UMAP 2",
    shape = "Known Lagurodon", fill = "Known Lagurodon"
  ) +
  theme_minimal(base_size = 13)
ggsave(file.path(output_dir, "tooth_shape_umap_R.png"), p_umap,
       width = 9, height = 7, dpi = 220)

print(p_umap)

z <- scale(shape)
heat <- data.frame(
  specimen = rep(features$specimen, each = length(shape_names)),
  feature = rep(shape_names, times = nrow(features)),
  z_score = as.vector(t(z))
)
p_heat <- ggplot(heat, aes(specimen, feature, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish) +
  labs(
    title = "Standardized scale-free morphology",
    subtitle = paste("Known reference:", known_name),
    x = NULL, y = NULL, fill = "z score"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "tooth_feature_heatmap_R.png"), p_heat,
       width = 11, height = 6, dpi = 220)

print(p_heat)

library(dplyr)
# plot sample Lagurodon and 2750 versus feature
ggplot(heat %>% filter(specimen %in% c("Lagurodon 2959", "2750", "2624")), aes(feature, z_score, color = specimen)) +
         geom_point(size = 2) +
         geom_line(aes(group = specimen), linewidth = 0.8) +
         scale_color_manual(values = c("gold", "#2878B5","red")) +
         labs(
           title = "Comparison of known Lagurodon with two specimens",
           x = "Morphological feature", y = "Standardized z score",
           color = "Specimen"
         ) +
         theme_minimal(base_size = 12) +
         theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(similarity, row.names = FALSE)
message("Results written to ", normalizePath(output_dir))
message("CAUTION: morphological proximity is not evidence of species identity.")
