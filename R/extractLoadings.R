#' extractLoadings
#'
#' @description Extract gene loadings for PCs of interest and calculate percent variance at the per-gene level
#'
#' @import ComplexHeatmap
#' @import scales
#' @import RColorBrewer
#' @import magick
#' @import grid
#'
#' @param variance A vector of gene variances (can calculate using RGenEDA::plotVariance)
#' @param pcRes Output list from RGenEDA::generatePCs
#' @param component Principal component to use (i.e., PC1)
#' @param nfeatures Number of features principal components were generated on
#'
#' @returns A data-frame containing gene loading information and percent variation per gene
#' @export


extractLoadings <- function(variance, pcRes, component, nfeatures) {
  loadings <- pcRes[["Eigenvectors"]][,component]
  gene_names <- names(variance)[1:nfeatures]

  #== Get gene level variance
  pc_var <- pcRes[["percent_var"]][component]
  pc_var_numeric <- as.numeric(gsub("%","",pc_var))/100
  gene_percent <- (loadings^2 * pc_var_numeric) / sum(loadings^2 * pc_var_numeric) * 100

  df <- data.frame(
    Gene = gene_names,
    Loading = loadings,
    PercentVariance = gene_percent
  )

  return(df)
}
