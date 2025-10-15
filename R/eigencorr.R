#' eigencorr
#'
#' @description Calculate Eigen-correlations from a `geneda` object and plot a
#' publication-quality heatmap using ggplot2 with numeric labels and
#' significance stars. For continuous variables, Pearson correlation is used.
#' P values follow this convention: p < 0.001, p < 0.01, p < 0.05 = three, two,
#' one stars, respectively.
#'
#' Requires PCA loadings in `DimReduction(object)`; does not recompute PCA.
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import grid
#'
#' @param object A `geneda` object containing `normalized` and `metadata`, and
#'   optionally `DimReduction` loadings.
#' @param NUM_PCS Number of principal components to correlate.
#' @param meta_cols Optional character vector of metadata column names to include.
#'   Defaults to all metadata columns.
#'
#' @returns A list with elements: cor_matrix, pval_matrix, stars, plot (ggplot)
#' @export
eigencorr <- function(object, NUM_PCS = 10, meta_cols = NULL) {
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

  #----- Obtain PCs from stored DimReduction (required)
  if (!(length(object@DimReduction) > 0L && "Loadings" %in% names(object@DimReduction))) {
    stop("DimReduction slot is empty. Please run RunPCA() before calling eigencorr().")
  }
  pcs_mat <- as.data.frame(object@DimReduction[["Loadings"]])
  # Align PCs to metadata order if necessary
  if (!identical(rownames(pcs_mat), rownames(META_num))) {
    pcs_mat <- pcs_mat[rownames(META_num), , drop = FALSE]
  }
  pcs <- pcs_mat[, seq_len(min(ncol(pcs_mat), NUM_PCS)), drop = FALSE]

  #----- Create empty matrices to store correlations and p-values
  cor_matrix <- matrix(NA, ncol = ncol(pcs), nrow = ncol(META_num))
  pval_matrix <- matrix(NA, ncol = ncol(pcs), nrow = ncol(META_num))
  rownames(cor_matrix) <- rownames(pval_matrix) <- colnames(META_num)
  colnames(cor_matrix) <- colnames(pval_matrix) <- colnames(pcs)

  #----- Loop through META variables and calculate correlation for each PC
  for (meta_var in colnames(META_num)) {
    for (pc in colnames(pcs)) {
      test <- cor.test(pcs[[pc]], META_num[[meta_var]], method = "pearson", use = "complete.obs")
      cor_matrix[meta_var, pc] <- test$estimate
      pval_matrix[meta_var, pc] <- test$p.value
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
  df <- expand.grid(Meta = rownames(cor_matrix), PC = colnames(cor_matrix), stringsAsFactors = FALSE)
  # Element-wise matrix lookup using paired row/col indices (avoids outer product)
  row_idx <- match(df$Meta, rownames(cor_matrix))
  col_idx <- match(df$PC, colnames(cor_matrix))
  df$Correlation <- cor_matrix[cbind(row_idx, col_idx)]
  df$PValue <- pval_matrix[cbind(row_idx, col_idx)]
  df$Stars <- stars[cbind(row_idx, col_idx)]

  #----- Set colors similar to RdBu reversed
  heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(255)

  #----- Build heatmap with labels and stars
  p <- ggplot2::ggplot(df, aes(x = PC, y = Meta, fill = Correlation)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = heatmap_colors, limits = c(-1, 1), name = "Correlation") +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 4) +
    geom_text(aes(label = Stars), nudge_y = -0.25, size = 6) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
    )

  #----- Return eigencorrelation data and plot
  return(list(
    cor_matrix = cor_matrix,
    pval_matrix = pval_matrix,
    stars = stars,
    plot = p
  ))
}



