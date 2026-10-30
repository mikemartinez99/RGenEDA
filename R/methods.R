#' Show method for geneda
#'
#' Prints a short summary of a `geneda` object: the number of features,
#' samples, and HVGs, whether PCA results and counts are present, and the
#' names of stored DEG tables.
#'
#' @param object A `geneda` object
#' @return Invisibly returns `object`; called for its printed summary.
#'
#' @examples
#' data(Snail1KO)
#' counts <- as.matrix(Snail1KO$rawCounts)
#' counts <- counts[rowSums(counts >= 10) >= 3, ]
#' meta <- Snail1KO$metadata
#' rownames(meta) <- meta$Sample
#' norm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
#' obj <- GenEDA(normalized = norm, metadata = meta, counts = counts)
#' obj <- RunPCA(obj, nfeatures = 1000)
#' show(obj)
setMethod(
    "show",
    signature = "geneda",
    definition = function(object) {
        cat("geneda object\n")
        cat(sprintf("  features: %d\n", nrow(object@normalized)))
        cat(sprintf("  samples:  %d\n", ncol(object@normalized)))
        cat(sprintf("  HVGs: %d\n", length(object@HVGs)))
        cat(sprintf("  DimReduction: %d\n", length(object@DimReduction)))
        cat(
            sprintf(
                "  counts: %s\n",
                if (is.null(object@counts)) "NULL" else "present"
            )
        )
        deg_status <- if (length(object@DEGs) == 0L) {
            "none"
        } else {
            paste(names(object@DEGs), collapse = ", ")
        }
        cat(sprintf("  DEGs: %s\n", deg_status))
        invisible(object)
    }
)
