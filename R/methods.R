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
    cat(sprintf("  counts: %s\n", if (is.null(object@counts)) "NULL" else "present"))
    deg_status <- if (is.null(object@DEGs$DEG)) "NULL" else {
      filtered_names <- setdiff(names(object@DEGs), "DEG")
      if (length(filtered_names) > 0L) {
        paste0("present (", length(filtered_names), " filtered: ", paste(filtered_names, collapse = ", "), ")")
      } else {
        "present (unfiltered only)"
      }
    }
    cat(sprintf("  DEGs: %s\n", deg_status))
    invisible(object)
  }
)


