#' Access counts matrix (avoid name clash with DESeq2::counts)
#' @param object A `geneda` object
#' @return Matrix or NULL
#' @export
#'
#' @examples
#' data(Snail1KO)
#' counts <- as.matrix(Snail1KO$rawCounts)
#' counts <- counts[rowSums(counts >= 10) >= 3, ]
#' meta <- Snail1KO$metadata
#' rownames(meta) <- meta$Sample
#' norm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
#' obj <- GenEDA(normalized = norm, metadata = meta, counts = counts)
#' getCounts(obj)[1:5, ]
getCounts <- function(object) {
    object@counts
}

#' Access normalized matrix (avoid name clash with generics)
#' @param object A `geneda` object
#' @return Matrix
#' @export
#'
#' @examples
#' data(Snail1KO)
#' counts <- as.matrix(Snail1KO$rawCounts)
#' counts <- counts[rowSums(counts >= 10) >= 3, ]
#' meta <- Snail1KO$metadata
#' rownames(meta) <- meta$Sample
#' norm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
#' obj <- GenEDA(normalized = norm, metadata = meta, counts = counts)
#' getNormalized(obj)[1:5, ]
getNormalized <- function(object) {
    object@normalized
}

#' Access metadata (avoid name clash with S4 generics)
#' @param object A `geneda` object
#' @return data.frame
#' @export
#'
#' @examples
#' data(Snail1KO)
#' counts <- as.matrix(Snail1KO$rawCounts)
#' counts <- counts[rowSums(counts >= 10) >= 3, ]
#' meta <- Snail1KO$metadata
#' rownames(meta) <- meta$Sample
#' norm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
#' obj <- GenEDA(normalized = norm, metadata = meta, counts = counts)
#' getMetadata(obj)
getMetadata <- function(object) {
    object@metadata
}

#' Access HVG IDs
#' @param object A `geneda` object
#' @return Character vector
#' @export
#'
#' @examples
#' data(Snail1KO)
#' counts <- as.matrix(Snail1KO$rawCounts)
#' counts <- counts[rowSums(counts >= 10) >= 3, ]
#' meta <- Snail1KO$metadata
#' rownames(meta) <- meta$Sample
#' norm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
#' obj <- GenEDA(normalized = norm, metadata = meta, counts = counts)
#' obj <- FindVariableFeatures(obj, nfeatures = 500)
#' head(HVGs(obj))
HVGs <- function(object) {
    object@HVGs
}

#' Access dimensional reduction results list
#' @param object A `geneda` object
#' @return List with `Loadings`, `Eigenvectors`, `percent_var`
#' @export
#'
#' @examples
#' data(Snail1KO)
#' counts <- as.matrix(Snail1KO$rawCounts)
#' counts <- counts[rowSums(counts >= 10) >= 3, ]
#' meta <- Snail1KO$metadata
#' rownames(meta) <- meta$Sample
#' norm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
#' obj <- GenEDA(normalized = norm, metadata = meta, counts = counts)
#' obj <- FindVariableFeatures(obj, nfeatures = 1000)
#' obj <- RunPCA(obj)
#' pca <- DimReduction(obj)
#' names(pca)
#' head(pca$Scores[, 1:3])
DimReduction <- function(object) {
    object@DimReduction
}

#' Access DEGs container
#' @param object A `geneda` object
#' @param assay The name of the DEGs slot to return
#' @return The DEG data.frame stored under `assay`, or NULL if it is not set
#' @export
#'
#' @examplesIf requireNamespace("DESeq2", quietly = TRUE)
#' data(Snail1KO)
#' counts <- as.matrix(Snail1KO$rawCounts)
#' counts <- counts[rowSums(counts >= 10) >= 3, ]
#' meta <- Snail1KO$metadata
#' rownames(meta) <- meta$Sample
#' meta$Genotype <- factor(meta$Genotype, levels = c("WT", "Snai1_KO"))
#' norm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
#' obj <- GenEDA(normalized = norm, metadata = meta, counts = counts)
#'
#' # Differential expression with DESeq2 (Snai1_KO vs WT)
#' dds <- DESeq2::DESeqDataSetFromMatrix(counts, meta, design = ~Genotype)
#' dds <- DESeq2::DESeq(dds, quiet = TRUE)
#' res <- as.data.frame(DESeq2::results(dds))
#' obj <- SetDEGs(obj, res, assay = "KO_vs_WT")
#' head(DEGs(obj, "KO_vs_WT"))
DEGs <- function(object, assay) {
    object@DEGs[[assay]]
}
