## make-Snail1KO.R
##
## Builds data/Snail1KO.rda, the example dataset shipped with RGenEDA.
##
## Source
##   Study:      Maturi et al., 2018. Atlas title: "Ampliseq transcriptomic
##               array of REF seq. genes in Snai1 knock out triple negative
##               breast cancer cells".
##   Accession:  ArrayExpress E-MTAB-5244 (ENA ERP019920)
##   Organism:   Homo sapiens, HS578T breast carcinoma cell line
##   Samples:    6 ENA runs, ERR1736465-ERR1736470
##               (3 wild type, 3 Snai1 knock-out)
##
## Counts were downloaded from EMBL-EBI Expression Atlas, which reprocessed
## E-MTAB-5244 with iRAP 1.0.3: reads aligned with TopHat2 2.1.1 to the
## GRCh38 reference (Ensembl release 107), and gene-level counts computed with
## HTSeq 0.11.0 (as listed on the Atlas "Supplementary Information" page).
## The file gives raw (unnormalized) integer read counts per Ensembl gene ID,
## with a "Gene Name" (gene symbol) column.
##
##   URL:  https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/
##         resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts
##   Downloaded:  2026-11-22 (Atlas experiment last updated 2025-01-23)
##   Size / MD5:  2,201,432 bytes / 4d7b915560c241ac562b267a66efdfff
##                61,860 Ensembl genes x 6 runs
##
## Licence
##   Expression Atlas data are released under the Creative Commons Attribution
##   4.0 International licence (CC BY 4.0):
##   https://www.ebi.ac.uk/gxa/licence.html
##
## Processing (no other filtering or transformation is done)
##   1. Gene identifiers: rows are relabelled from Ensembl gene IDs to the
##      Atlas "Gene Name" (gene symbol) column.
##   2. Genes are kept in the Atlas file order (sorted by Ensembl gene ID).
##      For each gene name only the FIRST row is kept; later rows with the same
##      name are dropped. This removes 1,194 rows with duplicated symbols (for
##      example Y_RNA, U3, SNORA63) and 21,302 of the 21,303 genes with no
##      symbol. The first symbol-less gene (ENSG00000083622, zero counts in all
##      samples) remains, with the empty string "" as its row name.
##      Result: 39,364 genes x 6 samples.
##   3. Counts are stored as an integer data.frame (genes x samples), with
##      ENA run accessions as column names.
##   4. Sample metadata come from the Atlas experiment-design file for
##      E-MTAB-5244 (columns "Sample Characteristic[disease]", "[genotype]",
##      "[cell line]" and "[cell type]"). Genotype is recoded to short labels:
##      "wild type genotype" -> "WT", "Snai1 knockout" -> "Snai1_KO". Cell line
##      and cell type are joined as "HS578T - epithelial cell".
##      Design URL: https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/
##                  resources/ExperimentDesignFile.RnaSeq/experiment-design
##
## Usage (from the package root):
##   Rscript inst/scripts/make-Snail1KO.R [output.rda]
## The default output is data/Snail1KO.rda.

args <- commandArgs(trailingOnly = TRUE)
outFile <- if (length(args) >= 1L) {
    args[[1]]
} else {
    file.path("data", "Snail1KO.rda")
}

url <- paste0(
    "https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/",
    "DifferentialSecondaryDataFiles.RnaSeq/raw-counts"
)
expectedMd5 <- "4d7b915560c241ac562b267a66efdfff"

#----- Download the Expression Atlas raw-count table
tsv <- tempfile(fileext = ".tsv")
utils::download.file(url, tsv, mode = "wb", quiet = TRUE)
md5 <- unname(tools::md5sum(tsv))
if (!identical(md5, expectedMd5)) {
    warning(
        "Downloaded file MD5 (", md5, ") differs from the file used to ",
        "build Snail1KO (", expectedMd5, "). Expression Atlas may have ",
        "reprocessed E-MTAB-5244; the result may not match data/Snail1KO.rda."
    )
}

raw <- utils::read.delim(tsv, check.names = FALSE, stringsAsFactors = FALSE)
samples <- paste0("ERR", 1736465:1736470)
stopifnot(
    identical(colnames(raw), c("Gene ID", "Gene Name", samples)),
    nrow(raw) == 61860L
)

#----- Steps 1-2: relabel by gene symbol, keep the first row per symbol
keep <- !duplicated(raw[["Gene Name"]])
rawCounts <- raw[keep, samples]
rownames(rawCounts) <- raw[["Gene Name"]][keep]

#----- Step 3: integer counts
rawCounts[] <- lapply(rawCounts, as.integer)
stopifnot(identical(dim(rawCounts), c(39364L, 6L)))

#----- Step 4: sample metadata from the Atlas experiment design
designUrl <- paste0(
    "https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/",
    "ExperimentDesignFile.RnaSeq/experiment-design"
)
designFile <- tempfile(fileext = ".tsv")
utils::download.file(designUrl, designFile, mode = "wb", quiet = TRUE)
design <- utils::read.delim(
    designFile,
    check.names = FALSE, stringsAsFactors = FALSE
)
design <- design[match(samples, design[["Run"]]), ]
stopifnot(identical(design[["Run"]], samples))

characteristic <- function(x) design[[paste0("Sample Characteristic[", x, "]")]]
genotypeLabels <- c("wild type genotype" = "WT", "Snai1 knockout" = "Snai1_KO")
stopifnot(all(characteristic("genotype") %in% names(genotypeLabels)))

metadata <- data.frame(
    Sample = samples,
    Disease = characteristic("disease"),
    Genotype = unname(genotypeLabels[characteristic("genotype")]),
    CellLine_CellType = paste(
        characteristic("cell line"), characteristic("cell type"),
        sep = " - "
    ),
    row.names = samples
)

Snail1KO <- list(rawCounts = rawCounts, metadata = metadata)

dir.create(dirname(outFile), showWarnings = FALSE, recursive = TRUE)
save(Snail1KO, file = outFile, compress = "bzip2")
message("Wrote ", outFile)
