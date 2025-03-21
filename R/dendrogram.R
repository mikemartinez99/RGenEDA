#' dendrogram
#'
#' @description Plot sample dendrogram colored by metadata feature.
#'
#' @importFrom dendextend color_branches
#' @import scales
#'
#' @param MAT A data matrix where rows are features and columns are samples
#' @param META Metadata as a dataframe consisting of the samples provided in MAT. Samples should be rows and columns are metadata features
#' @param FEATURE A specific metadata feature of which the dendrogram branches will be colored by
#' @param COLORS A named color vector with one entry for every unique level in your FEATURE
#' @param OUTPUT A directory path to where the output should be saved. Path needs to end in a "/"
#'
#' @returns A sample dendrogram in tiff format at the specified output path
#'
#' @examples color_vector <- c("A" = "red", "B" = "blue")
#' @examples output_folder <- c("Path/to/outputs/")
#' @examples dendrogram(myDataMatrix, metadataDF, "some_column", color_vector, output_folder)

dendrogram <- function(MAT, META, FEATURE, COLORS, OUTPUT) {
  #----- Transpose data matrix so rows are samples and columns are features
  message("Transposing matrix...")
  sampleDists <- dist(t(MAT))
  #----- Hierarchical clustering
  message("Performing hierarchical clustering...")
  hclustRes <- hclust(sampleDists)
  #----- Plot dendrogram
  message("Plotting dendrogram...")
  dend <- as.dendrogram(hclustRes)
  #----- Define a color palette
  color_palette <- COLORS
  #----- Map the metadata to colors and color dendrogram branches
  colors <- color_palette[META[labels(dend), FEATURE]]
  dend2 <- dendextend::color_branches(dend, col = colors)
  #----- Plot and save file
  message(paste0("Saving dendrogram to ", OUTPUT))
  fileName <- paste("HClust", FEATURE, "tiff", sep = ".")
  #----- Initialize a plot
  tiff(paste(OUTPUT, fileName),
       width = 12, height = 16, units = "in", res = 200)
  #----- Adjust the line thickness
  par(lwd = 4.5)
  #----- Plot the dendrogram with colored branches
  plot(dend2,
       main = paste("HClust", gsub("_", " ", FEATURE), sep = " "),
       ylab = "Height",
       cex = 0.6)
  #----- Add a legend
  legend("topright",  # You can change the position as needed
         legend = names(color_palette),  # Group names
         col = color_palette,  # Corresponding colors
         pch = 15,  # Type of point (15 is a filled square)
         bty = "n",
         cex = 2.5)  # Text size
  #-----Close the plot
  dev.off()
  message("Dendrogram successfully plotted!")
}
