#' distanceHeatmap
#'
#' @description Calculate sample to sample distances using Euclidean distances. Create annotataion-rich heatmap to assess data trends.
#'
#' @importFrom pheatmap pheatmap
#' @import scales
#' @import RColorBrewer
#'
#' @param MAT A data matrix where rows are features and columns are samples
#' @param META Metadata as a dataframe consisting of the samples provided in MAT. Samples should be rows and columns are metadata features
#' @param FEATURE A vector of features to be included in the heatmap annotations
#' @param PALETTES A list where each list slot is a named vector of metadata levels and colors. Each feature present in FEATURE needs to have its own list slot
#' @param OUTPUT A directory path to where the output should be saved. Path needs to end in a "/"
#'
#' @returns A distance heatmap in png format at the specified output path
#'
#' @examples timepoint_colors <- c("Early" = "red", "Late" = "blue")
#' @examples treatment_colors <- c("Mock" = "pink", "HSV" = "green")
#' @examples color_list <- list(timepoint_colors, treatment_colors)
#' @examples output_folder <- c("Path/to/outputs/")
#' @examples dendrogram(myDataMatrix, metadataDF, c("Timepoint", "Treatment"), color_list, output_folder)

distanceHeatmap <- function(MAT, META, FEATURES = c(...), PALETTES, OUTPUT) {
  #----- Calculate the sample distances
  message("Transposing data matrix...")
  matT <- t(MAT)
  sampleDists <- dist(matT)
  #----- Convert distances into matrix
  distMat <- as.matrix(sampleDists)
  #----- Set the diagonal to zero (self:self comparisons)
  diag(distMat) <- 0
  #----- Create a data frame for annotations based on features provided
  featName <- FEATURES
  message(paste0("Extracting ", featName, " as features for plotting..."))
  groupdf <- META[ ,colnames(META) %in% featName, drop = FALSE]
  groupdf <- as.data.frame(groupdf)
  rownames(groupdf) <- rownames(META)
  #----- Create a color mapping data frame for annotations
  annotation_colors <- list()  # Initialize list to store colors for each feature
  for (feature in FEATURES) {
    if (feature %in% names(PALETTES)) {
      #----- Get the color mapping for the feature
      colors <- PALETTES[[feature]]
      #----- Ensure the feature is treated as a factor
      groupdf[[feature]] <- as.factor(groupdf[[feature]])
      #----- Create a named vector of colors based on the factor levels
      annotation_colors[[feature]] <- colors[levels(groupdf[[feature]])]
    }
  }
  #----- Set colors for the heatmap
  heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  #----- Output folder
  fileName <- c("Sample_Distance_HM.tiff")
  message(paste0("Saving distance heatmap to ", OUTPUT))
  tiff(paste0(OUTPUT, fileName), width = 8, height = 8, units = "in", res = 150)
  #----- Visualize the heatmap with colored annotations
  H <- pheatmap::pheatmap(distMat,
                #clustering_distance_rows = sampleDists,
                annotation_col = groupdf,
                annotation_colors = annotation_colors,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                col = heatmap_colors,
                fontsize = 8,
                fontsize_row = 8,
                fontsize_col = 8,
                legend = TRUE)
  dev.off()
  message("Distance heatmap successfully generated!")
}

