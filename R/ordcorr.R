#' ordcorr
#'
#' @description Calculate correlations between NMDS ordination axes (from beta
#' diversity distances) and metadata, and plot a ggplot2 heatmap with numeric
#' labels and significance stars.
#'
#' For continuous variables, Pearson correlation is used. P values follow this
#' convention: p < 0.001, p < 0.01, p < 0.05 = three, two, one stars,
#' respectively.
#'
#' Default num_mds is 10, but function will internally set max number of NMDS
#' to be n-1 where n is your number of samples.
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import grid
#' @import vegan
#' @importFrom stats as.dist cor dist
#'
#' @param object A `geneda` object containing `normalized` and `metadata`.
#' @param num_mds Number of NMDS axes to correlate. Default = 10
#' @param meta_cols Optional character vector of metadata column names to include.
#'   Defaults to all metadata columns.
#' @param distance Distance metric for `vegdist` (default "bray").
#'
#' @returns A list with elements: cor_matrix, pval_matrix, stars, plot (ggplot),
#' and variance explained (squared correlations between ordination distances and
#' observed dissimilarities.)
#' @export
ordcorr <- function(object, num_mds = 10, meta_cols = NULL, distance = "bray") {
  stopifnot(methods::is(object, "geneda"))

  MAT <- object@normalized
  META <- object@metadata

  #----- Convert all META categories to numeric (copy to avoid modifying slot)
  META_num <- META
  for (i in colnames(META_num)) {
    if (!is.numeric(META_num[[i]])) {
      META_num[[i]] <- as.numeric(as.factor(META_num[[i]]))
    }
  }

  #----- Optionally subset metadata columns
  if (!is.null(meta_cols)) {
    missing_cols <- setdiff(meta_cols, colnames(META_num))
    if (length(missing_cols) > 0L) {
      stop(paste0("The following metadata columns were not found: ", paste(missing_cols, collapse = ", ")))
    }
    META_num <- META_num[, meta_cols, drop = FALSE]
  }

  #----- Ensure META rownames match those of MAT
  if (!all(rownames(META_num) == colnames(MAT))) {
    stop("Rownames of metadata must match the colnames of normalized matrix.")
  }

  #----- Compute distance matrix across samples (rows must be samples)
  distmat <- vegan::vegdist(t(MAT), method = distance)

  #----- Adjust k if it's larger than n - 1
  max_k <- nrow(as.matrix(distmat)) - 1
  if (num_mds > max_k) {
    message(sprintf("Reducing num_mds from %d to %d (max allowed for %d samples).", num_mds, max_k, max_k + 1))
    num_mds <- max_k
  }

  #----- NMDS ordination
  dte <- vegan::metaMDS(distmat, distance = distance, k = num_mds, trymax = 100, trace = FALSE)
  nmds_scores <- as.data.frame(vegan::scores(dte))

  #----- Approximate variance explained by each NMDS axis
  # Compute Euclidean distances between points in NMDS space for each k
  coords <- as.matrix(nmds_scores)
  original_dist <- as.dist(distmat)
  var_explained <- numeric(ncol(coords))

  for (i in seq_len(ncol(coords))) {
    reduced_dist <- dist(coords[, seq_len(i), drop = FALSE])
    r2 <- cor(original_dist, reduced_dist)^2
    var_explained[i] <- ifelse(i == 1, r2, r2 - sum(var_explained[1:(i - 1)]))
  }
  names(var_explained) <- colnames(nmds_scores)

  # Align to metadata order if necessary
  if (!identical(rownames(nmds_scores), rownames(META_num))) {
    nmds_scores <- nmds_scores[rownames(META_num), , drop = FALSE]
  }

  #----- Create empty matrices to store correlations and p-values
  cor_matrix <- matrix(NA, ncol = ncol(nmds_scores), nrow = ncol(META_num))
  pval_matrix <- matrix(NA, ncol = ncol(nmds_scores), nrow = ncol(META_num))
  rownames(cor_matrix) <- rownames(pval_matrix) <- colnames(META_num)
  colnames(cor_matrix) <- colnames(pval_matrix) <- colnames(nmds_scores)

  #----- Loop through META variables and calculate correlation for each NMDS axis
  for (meta_var in colnames(META_num)) {
    for (axis in colnames(nmds_scores)) {
      test <- cor.test(nmds_scores[[axis]], META_num[[meta_var]], method = "pearson", use = "complete.obs")
      cor_matrix[meta_var, axis] <- test$estimate
      pval_matrix[meta_var, axis] <- test$p.value
    }
  }

  #----- Remove NA
  cor_matrix <- na.omit(cor_matrix)
  pval_matrix <- na.omit(pval_matrix)

  #----- Create significance stars based on p-value thresholds
  stars <- ifelse(pval_matrix < 0.001, "***",
                  ifelse(pval_matrix < 0.01, "**",
                         ifelse(pval_matrix < 0.05, "*", "")))

  #----- Prepare data for ggplot (long format)
  df <- expand.grid(Meta = rownames(cor_matrix), Axis = colnames(cor_matrix), stringsAsFactors = FALSE)
  row_idx <- match(df$Meta, rownames(cor_matrix))
  col_idx <- match(df$Axis, colnames(cor_matrix))
  df$Correlation <- cor_matrix[cbind(row_idx, col_idx)]
  df$PValue <- pval_matrix[cbind(row_idx, col_idx)]
  df$Stars <- stars[cbind(row_idx, col_idx)]

  #----- Set colors similar to RdBu reversed
  heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(255)

  #----- Build heatmap with labels and stars
  p <- ggplot2::ggplot(df, aes(x = Axis, y = Meta, fill = Correlation)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = heatmap_colors, limits = c(-1, 1), name = "Correlation") +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 4) +
    geom_text(aes(label = Stars), nudge_y = -0.30, size = 5) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )

  #----- Return ord-correlation data and plot
  return(list(
    cor_matrix = cor_matrix,
    pval_matrix = pval_matrix,
    stars = stars,
    plot = p,
    variance_explained = var_explained
  ))
}


