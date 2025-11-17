#' Set unfiltered DEGs on the object
#'
#' @param object A `geneda` object
#' @param deg_table A data.frame of DESeq2-like results containing at least
#'   columns `log2FoldChange` and `padj`
#' @param assay What to name the DEG slot
#' @return Updated `geneda` object with `DEGs$DEG` set
#' @export
SetDEGs <- function(object, deg_table, assay) {
  stopifnot(methods::is(object, "geneda"))
  stopifnot(is.data.frame(deg_table))
  req_cols <- c("log2FoldChange", "padj")
  missing_cols <- setdiff(req_cols, colnames(deg_table))
  if (length(missing_cols) > 0L) {
    stop(paste0("DEG table must contain columns: ", paste(req_cols, collapse = ", ")))
  }
  if (is.null(rownames(deg_table))) {
    stop("DEG table does not have genes as rownames!")
  }

  if (assay %in% names(object@DEGs)) {
    stop(paste("Assay", assay, "already exists!"))
  } else {
    object@DEGs[[assay]] <- deg_table
  }
  validObject(object)
  object
}

#' Filter DEGs by padj and absolute log2FoldChange
#'
#' Filters the unfiltered DEGs in `DEGs$DEG` and stores the filtered results
#' in a new named slot `DEGs[[saveAssay]]`. Multiple filtered result sets can
#' be stored with different assay names.
#'
#' @param object A `geneda` object with `DEGs$DEG` set
#' @param assay The DEG slot to filter
#' @param alpha Adjusted p-value threshold (<=)
#' @param log2FC_thresh Absolute log2 fold change threshold (>=)
#' @param saveAssay Character name for the filtered result set (e.g., "padj05_lfc1")
#' @return Updated `geneda` object with filtered results stored in `DEGs[[saveAssay]]`
#' @importFrom methods validObject
#' @export
FilterDEGs <- function(object, assay, alpha = 0.05, log2FC_thresh = 1.0, saveAssay) {
  stopifnot(methods::is(object, "geneda"))
  if (is.null(object@DEGs[[assay]])) {
    stop(paste("Assay", assay, "is NULL. Use SetDEGs() first."))
  }
  if (missing(saveAssay) || is.null(saveAssay) || !is.character(saveAssay) || length(saveAssay) != 1L) {
    stop("saveAssay must be a single character string.")
  }
  df <- object@DEGs[[assay]]
  req_cols <- c("log2FoldChange", "padj")
  missing_cols <- setdiff(req_cols, colnames(df))
  if (length(missing_cols) > 0L) {
    stop(paste0("DEG table must contain columns: ", paste(req_cols, collapse = ", ")))
  }
  filt <- stats::complete.cases(df$log2FoldChange, df$padj) &
    abs(df$log2FoldChange) >= log2FC_thresh & df$padj <= alpha
  object@DEGs[[saveAssay]] <- df[filt, , drop = FALSE]
  validObject(object)
  object
}

#' Find DEGs Overlapping With HVGs
#'
#' Looks for DEGs after filtering that overlap with HVGs.
#'
#' @param object A `geneda` object with `DEGs$DEG` set
#' @param assay The DEG slot to pull results from (typically a result of `FilterDEGs`)
#' @param direction One of positive, negative, or both (indicating the direction
#' of DEG fold change to intersect with HVGs)
#'
#' @export
FindHVDEGs <- function(object, assay, direction = c("positive", "negative", "both")) {
  stopifnot(methods::is(object, "geneda"))

  # Validate direction argument
  direction <- match.arg(direction)

  # Access DEGs
  if (is.null(object@DEGs[[assay]])) {
    stop(paste("Assay", assay, "is NULL. Use SetDEGs() or FilterDEGs() first!"))
  } else {
    degDF <- DEGs(object, assay)
    degDF$Gene <- rownames(degDF)
  }

  # Access HVGs
  if (length(object@HVGs) == 0L){
    stop("HVG Slot is empty. Use FindVariableFeatures() first!")
  } else {
    hvgList <- HVGs(object)
  }

  # Get overlap
  pos <- degDF[degDF$log2FoldChange > 0,]$Gene
  neg <- degDF[degDF$log2FoldChange < 0,]$Gene
  posInt <- intersect(pos, hvgList)
  negInt <- intersect(neg, hvgList)

  if (direction == "positive") {
    if (length(posInt) == 1) {
      message(paste(length(posInt), "+log2FC hvDEG found."))
    } else {
      message(paste(length(posInt), "+log2FC hvDEGs found."))
    }
    return(posInt)

  } else if (direction == "negative") {
    if (length(negInt) == 1) {
      message(paste(length(negInt), "-log2FC hvDEG found."))
    } else {
      message(paste(length(negInt), "-log2FC hvDEGs found."))
    }
    return(negInt)

  } else {
    sumhvDEGs <- length(posInt) + length(negInt)
    message(paste(length(posInt), "+log2FC hvDEGs found."))
    message(paste(length(negInt), "-log2FC hvDEGs found."))
    message(paste(sumhvDEGs, "Total hvDEGs found."))
    return(list(
      positive = posInt,
      negative = negInt,
      both = c(posInt, negInt)
    ))

  }
}

#' Summarize DEGs Across Selected Assays
#'
#' @description
#' Summarizes the number of DEGs in each assay of a `geneda` object according to thresholds
#' for adjusted p-value and log2 fold-change, separately reporting upregulated and downregulated genes.
#'
#' @param object A `geneda` object containing DEGs in its `DEGs` slot.
#' @param alpha Numeric; adjusted p-value threshold. Default is 0.05.
#' @param lfc1 Numeric; lower log2 fold-change threshold. Default is 1.
#' @param lfc2 Numeric; higher log2 fold-change threshold. Default is 2.
#' @param assays Optional character vector of assay names to summarize. Defaults to all assays in `object@DEGs`.
#'
#' @return A data.frame of DEG counts by assay and criteria.
#' @export
SummarizeDEGs <- function(object, alpha = 0.05, lfc1 = 1, lfc2 = 2, assays = NULL) {
  stopifnot(methods::is(object, "geneda"))
  if (length(object@DEGs) == 0) stop("No DEGs found in object.")

  # Use all assays if none are specified
  if (is.null(assays)) {
    assays <- names(object@DEGs)
  }

  criteria <- c(
    paste0("padj<", alpha, " & log2FC >", lfc1, " (Up)"),
    paste0("padj<", alpha, " & log2FC <", -lfc1, " (Down)"),
    paste0("padj<", alpha, " & log2FC >", lfc2, " (Up)"),
    paste0("padj<", alpha, " & log2FC <", -lfc2, " (Down)")
  )

  out <- matrix(0, nrow = length(criteria), ncol = length(assays),
                dimnames = list(criteria, assays))

  for (assay in assays) {
    if (!assay %in% names(object@DEGs)) {
      warning(paste("Assay", assay, "not found in DEGs slot; skipping."))
      next
    }
    df <- DEGs(object, assay)
    df <- na.omit(df)

    # Ensure numeric
    df$log2FoldChange <- as.numeric(df$log2FoldChange)
    df$padj <- as.numeric(df$padj)

    # Count DEGs by criteria
    out[1, assay] <- sum(df$padj < alpha & df$log2FoldChange > lfc1)
    out[2, assay] <- sum(df$padj < alpha & df$log2FoldChange < -lfc1)
    out[3, assay] <- sum(df$padj < alpha & df$log2FoldChange > lfc2)
    out[4, assay] <- sum(df$padj < alpha & df$log2FoldChange < -lfc2)
  }

  return(as.data.frame(out))
}













