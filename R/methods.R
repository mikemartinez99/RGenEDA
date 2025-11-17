#' Show method for geneda
#' @param object A `geneda` object
setMethod(
  "show",
  signature = "geneda",
  definition = function(object) {
    cat("geneda object\n")
    cat(sprintf("  features: %d\n", nrow(object@normalized)))
    cat(sprintf("  samples:  %d\n", ncol(object@normalized)))
    cat(sprintf("  HVGs: %d\n", length(object@HVGs)))
    cat(sprintf("  DimReduction: %d\n", length(object@DimReduction)))
    cat(sprintf("  counts: %s\n", if (is.null(object@counts)) "NULL" else "present"))
    deg_status <- if (is.null(object@DEGs)){
      "NULL"
    } else {
      names(object@DEGs)

    }
    cat(sprintf("  DEGs: %s\n", deg_status))
    invisible(object)
  }
)


