#' eigencorr_MGX
#'
#' @description Calculate Eigen-correlations for MGX data to correlate axes of variation with Bray-Curtis distances. For continuous variables, the function uses Pearson correlation. Pvalues follow this convention: p < 0.001, p < 0.01, p < 0.05 = three stars, two stars, one star, respectively
#'
#' @import ComplexHeatmap
#' @import scales
#' @import RColorBrewer
#' @import magick
#' @import grid
#' @import vegan
#'
#' @param MAT A data matrix where rows are features and columns are samples
#' @param META Metadata as a dataframe consisting of the samples provided in MAT. Samples should be rows and columns are metadata features. If you only want to explore certain metadata features, create a subset of your metadata to only include columns of interest.
#' @param NUM_mds Number of NMDS to use.
#' @param OUTPUT directory path to where the output should be saved. Path needs to end in a "/"
#'
#' @returns A list consisting of 3 slots: cor_matrix which holds correlation results, pval_matrix which holds p-value results, and stars which holds significance stars for plotting
#' @export

eigencorr_MGX <- function(mgx, metadata, num_mds = 10, OUTPUT) {

  #----- Convert all metadata categories to numeric
  for (i in colnames(metadata)) {
    if (!is.numeric(metadata[[i]])) {
      metadata[[i]] <- as.numeric(as.factor(metadata[[i]]))
    }
  }

  #----- Get distance matrix
  distmat <- vegdist(mgx,
                     method = "bray")

  #----- Calculate MDS k=10 to reduce stress
  dte <- metaMDS(distmat, distance = "bray", k = num_mds, trymax = 100, trace = F)


  #----- Get the multidimensional scaling vectors (MDs)
  dte.scores <- as.data.frame(scores(dte))


  #----- Ensure metadata rownames match those of the vsd object
  if (!all(rownames(metadata) == rownames(mgx))) {
    stop("Rownames of metadata must match the rownames of the MGX object.")
  }

  #----- Create empty matrices to store correlations and p-values
  cor_matrix <- matrix(NA, ncol = ncol(dte.scores), nrow = ncol(metadata))
  pval_matrix <- matrix(NA, ncol = ncol(dte.scores), nrow = ncol(metadata))
  rownames(cor_matrix) <- rownames(pval_matrix) <- colnames(metadata)
  colnames(cor_matrix) <- colnames(pval_matrix) <- colnames(dte.scores)

  #----- Loop through metadata variables and calculate correlation for each PC
  for (meta_var in colnames(metadata)) {
    for (pc in colnames(dte.scores)) {
      if (is.numeric(metadata[[meta_var]])) {
        #----- For continuous variables, use Pearson correlation and calculate p-values
        test <- cor.test(dte.scores[[pc]], metadata[[meta_var]], method = "pearson", use = "complete.obs")
        cor_matrix[meta_var, pc] <- test$estimate
        pval_matrix[meta_var, pc] <- test$p.value
      } else if (is.factor(metadata[[meta_var]])) {
        #----- For categorical variables, use ANOVA to calculate F-statistic and p-value
        aov_res <- aov(dte.scores[[pc]] ~ metadata[[meta_var]])
        summary_aov <- summary(aov_res)
        cor_matrix[meta_var, pc] <- summary_aov[[1]][["F value"]][1]  # F-statistic
        pval_matrix[meta_var, pc] <- summary_aov[[1]][["Pr(>F)"]][1]  # p-value
      }
    }
  }

  #----- Remove NA
  cor_matrix <- na.omit(cor_matrix)
  pval_matrix <- na.omit(pval_matrix)

  #----- Create significance stars based on p-value thresholds
  stars <- ifelse(pval_matrix < 0.001, "***",
                  ifelse(pval_matrix < 0.01, "**",
                         ifelse(pval_matrix < 0.05, "*", "")))

  #----- Set colors for the heatmap
  heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)


  #----- Return list
  return_list <- list(cor_matrix = cor_matrix,
                      pval_matrix = pval_matrix,
                      stars = stars
  )

  #----- Generate OUTPUT file path
  fileName <- c("MGX_EigenCorrelations.tiff")
  message(paste0("Output path: ", paste0(OUTPUT, fileName)))

  #----- Initialize image
  tiff(paste0(OUTPUT, fileName), width = 10, height = 10, units = "in", res = 200)

  #----- Plot heatmap
  message("Plotting heatmap...")
  H <- ComplexHeatmap::Heatmap(cor_matrix,
                               cluster_rows = TRUE,
                               cluster_columns = FALSE,
                               cell_fun = function(j, i, x, y, width, height, fill) {
                                 grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y, gp = grid::gpar(fontsize = 20))
                                 if (stars[i, j] != "") {
                                   grid.text(stars[i, j], x, y + unit(5, "mm"), gp = grid::gpar(col = "black", fontsize = 24))
                                 }
                               },
                               name = "Correlation",
                               show_row_names = TRUE,
                               show_column_names = TRUE,
                               row_names_gp = grid::gpar(fontsize = 20),  # Increase row names size
                               column_names_gp = grid::gpar(fontsize = 20),
                               column_title = "",
                               row_title = "")
  draw(H)
  dev.off()
  message("Plotting complete!")

  #----- Return value
  return(return_list)

}
