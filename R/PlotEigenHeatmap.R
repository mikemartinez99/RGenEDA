#' PlotEigenHeatmap
#'
#' @description
#' Generates a ComplexHeatmap of the top and/or bottom genes contributing to a
#' specific principal component from a `geneda` object. Genes are selected based
#' on their eigenvector loadings, and expression values are Z-score scaled. The
#' function optionally annotates columns with metadata variables and adds a
#' row annotation reflecting gene loading magnitude and percent variance.
#'
#' @param object A `geneda` object containing normalized expression data and
#' optionally HVGs.
#' @param pc Character. Which principal component to visualize (default: "PC1").
#' @param n Integer. Number of genes to select from top or bottom loadings (default: 25).
#' @param direction Character. Whether to select "top", "bottom", or "both" ends
#' of the eigenvector distribution (default: "both").
#' @param annotate_by Character vector of metadata column names to include as
#' top annotations (optional).
#' @param annotate_colors Named list of color vectors for metadata annotation.
#' Names should match `annotate_by` column names.
#' @param show_row_names Logical. Whether to display gene names as row labels
#' (default: TRUE).
#'
#' @return A `ComplexHeatmap` object representing the heatmap of top genes.
#'
#' @details
#' The function extracts loadings for the specified principal component using
#' `extractLoadings()`. Genes are selected based on the chosen `direction` and
#' number `n`. Expression values are scaled across genes for visualization.
#' Row annotations indicate the magnitude and direction of gene loadings, while
#' column annotations can display metadata variables. A custom size legend
#' represents percent variance of each gene.
#'
#' @import ComplexHeatmap
#' @import circlize
#' @import scales
#' @import dplyr
#' @importFrom grid gpar unit
#'
#' @examples
#' \dontrun{
#' # Visualize top 25 genes for PC1 with metadata annotation
#' ht <- PlotTopGenesHeatmap(obj, pc = "PC1", n = 25, direction = "both",
#'                           annotate_by = c("Condition"),
#'                           annotate_colors = list(Condition = c("untreated" = "red", "treated" = "blue")))
#' draw(ht)
#' }
#'
#' @export


PlotEigenHeatmap <- function(object, pc = "PC1", n = 25, direction = c("both","top","bottom"),
                             annotate_by = NULL, annotate_colors = NULL, show_row_names = TRUE) {

  stopifnot(methods::is(object, "geneda"))
  direction <- match.arg(direction)

  # Calculate eigen vectors
  pc_df <- extractEigen(object, pc)

  # Filter top/bottom genes
  pc_col <- "EigenVector"
  pct_var_col <- "PctVar"

  if (direction == "top") {
    selected <- pc_df %>%
      arrange(desc(!!sym(pc_col))) %>%
      slice_head(n = n)
  } else if (direction == "bottom") {
    selected <- pc_df %>%
      arrange(!!sym(pc_col)) %>%
      slice_head(n = n)
  } else if (direction == "both") {
    top_genes <- pc_df %>%
      arrange(desc(!!sym(pc_col))) %>%
      slice_head(n = n)
    bottom_genes <- pc_df %>%
      arrange(!!sym(pc_col)) %>%
      slice_head(n = n)
    selected <- bind_rows(top_genes, bottom_genes)
  }

  selected <- selected %>%
    mutate(Direction = ifelse(!!sym(pc_col) > 0, "Positive", "Negative"),
           Gene = factor(Gene, levels = Gene[order(!!sym(pc_col))]))
  gene_order <- selected$Gene

  # Extract normalized counts from GenEDA object
  mat <- object@normalized
  # Keep only HVGs if available
  if (!is.null(object@HVGs)) {
    mat <- mat[rownames(mat) %in% object@HVGs, , drop=FALSE]
  }

  # Subset and scale for selected genes
  vst_mat <- mat[rownames(mat) %in% gene_order, , drop=FALSE]
  vst_mat <- vst_mat[gene_order, , drop=FALSE]
  vst_scaled <- t(scale(t(vst_mat)))

  # Row annotation for loadings
  row_ha <- rowAnnotation(
    EigenVec = anno_points(
      selected[[pc_col]],
      pch = 21,
      gp = gpar(
        col = "black",
        fill = colorRamp2(
          c(min(selected[[pc_col]]), 0, max(selected[[pc_col]])),
          c("firebrick", "white", "dodgerblue")
        )(selected[[pc_col]])
      ),
      size = unit(rescale(selected[[pct_var_col]], to = c(2,8)), "mm")
    ),
    width = unit(20, "mm")
  )

  # Custom size legend for Percent Variance
  size_legend <- Legend(
    title = "Percent Variance",
    at = pretty(selected[[pct_var_col]]),
    type = "points",
    pch = 21,
    size = unit(rescale(pretty(selected[[pct_var_col]]), to = c(2, 8)), "mm")
  )

  # Top annotation for metadata
  ha_top <- NULL
  if (!is.null(annotate_by)) {
    meta <- as.data.frame(object@metadata)

    # Keep only columns that exist
    annotate_by <- annotate_by[annotate_by %in% colnames(meta)]
    if (length(annotate_by) > 0) {
      col_colors <- list()
      if (!is.null(annotate_colors)) {
        # Check that all color names match metadata columns
        invalid_cols <- setdiff(names(annotate_colors), colnames(meta))
        if (length(invalid_cols) > 0) {
          warning(
            "The following entries in `annotate_colors` do not match any metadata columns and will be ignored: ",
            paste(invalid_cols, collapse = ", ")
          )
        }
      }

      for (colname in annotate_by) {
        col_data <- meta[[colname]]
        # Use user-provided colors if names match exactly
        if (!is.null(annotate_colors) && colname %in% names(annotate_colors)) {
          col_colors[[colname]] <- annotate_colors[[colname]]
        } else if (is.factor(col_data) || is.character(col_data)) {
          # Automatic color generation
          u <- unique(col_data)
          cols <- rainbow(length(u))
          names(cols) <- u
          col_colors[[colname]] <- cols
        } else {
          # Numeric columns
          col_colors[[colname]] <- circlize::colorRamp2(range(col_data, na.rm = TRUE), c("white", "red"))
        }
      }

      ha_top <- HeatmapAnnotation(
        df = meta[, annotate_by, drop = FALSE],
        col = col_colors
      )
    }
  }


  # Create ComplexHeatmap object
  ht <- Heatmap(
    vst_scaled,
    name = paste0(pc, " Expression (Z-score)"),
    left_annotation = row_ha,
    top_annotation = ha_top,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = show_row_names,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 9.5)
  )

  draw(ht, annotation_legend_list = list(size_legend))
  returnList <- c(ht, size_legend)
  return(returnList)



}
