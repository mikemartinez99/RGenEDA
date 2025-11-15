#' Access counts matrix (avoid name clash with DESeq2::counts)
#' @param object A `geneda` object
#' @return Matrix or NULL
#' @export
getCounts <- function(object) {
  object@counts
}

#' Access normalized matrix (avoid name clash with generics)
#' @param object A `geneda` object
#' @return Matrix
#' @export
getNormalized <- function(object) {
  object@normalized
}

#' Access metadata (avoid name clash with S4 generics)
#' @param object A `geneda` object
#' @return data.frame
#' @export
getMetadata <- function(object) {
  object@metadata
}

#' Access HVG IDs
#' @param object A `geneda` object
#' @return Character vector
#' @export
HVGs <- function(object) {
  object@HVGs
}

#' Access dimensional reduction results list
#' @param object A `geneda` object
#' @return List with `Loadings`, `Eigenvectors`, `percent_var`
#' @export
DimReduction <- function(object) {
  object@DimReduction
}

#' Access DEGs container
#' @param object A `geneda` object
#' @param assay The name of the DEGs slot to return
#' @return List with `DEG` (unfiltered data.frame) and optionally named filtered result sets
#' @export
DEGs <- function(object, assay) {
  object@DEGs[[assay]]
}
