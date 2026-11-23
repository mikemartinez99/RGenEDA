#' Snail1KO: RNA-seq counts from Snai1 knock-out HS578T cells
#'
#' Raw RNA-seq read counts and sample metadata for human HS578T triple-negative
#' breast carcinoma cells: 3 wild-type and 3 Snai1 knock-out samples, from
#' Maturi et al. (2018), ArrayExpress E-MTAB-5244 (ENA ERP019920).
#'
#' @format A list with two elements:
#' \describe{
#'   \item{rawCounts}{A `data.frame` of raw (unnormalized) integer read counts
#'     with 39,364 genes (rows) x 6 samples (columns). Row names are gene
#'     symbols (the Expression Atlas "Gene Name" for each Ensembl gene ID).
#'     Column names are ENA run accessions (ERR1736465 to ERR1736470).}
#'   \item{metadata}{A `data.frame` with 6 rows (samples, with run accessions
#'     as row names) and 4 columns: `Sample` (run accession), `Disease`
#'     ("breast carcinoma"), `Genotype` ("WT" or "Snai1_KO") and
#'     `CellLine_CellType` ("HS578T - epithelial cell").}
#' }
#'
#' @details
#' \itemize{
#'   \item Organism: *Homo sapiens*, HS578T cell line.
#'   \item Gene IDs: gene symbols. The source counts are per Ensembl gene ID
#'     (GRCh38, Ensembl release 107). Each Ensembl ID was replaced with its gene
#'     symbol, and for each symbol only the first gene in Ensembl ID order was
#'     kept. This removed 1,194 duplicated symbols and 21,302 genes with no
#'     symbol. One symbol-less gene (ENSG00000083622, zero counts) remains,
#'     with an empty row name `""`.
#'   \item Quantification: by Expression Atlas with iRAP 1.0.3 (TopHat2
#'     2.1.1 alignment, HTSeq 0.11.0 counting).
#'   \item No other filtering or normalization was done.
#' }
#' The script that builds this object from the Expression Atlas download is
#' in the installed package:
#' `system.file("scripts", "make-Snail1KO.R", package = "RGenEDA")`.
#'
#' @source Counts: EMBL-EBI Expression Atlas, experiment E-MTAB-5244
#'   (\url{https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5244}), raw-count
#'   file downloaded 2026-11-22. Raw data: ArrayExpress E-MTAB-5244 / ENA
#'   ERP019920.
#'
#' @section License:
#' Expression Atlas data are released under the Creative Commons Attribution
#' 4.0 International licence (CC BY 4.0;
#' \url{https://www.ebi.ac.uk/gxa/licence.html}). Please cite Maturi et al.
#' (2018) and Expression Atlas when you use these data.
#'
#' @usage data(Snail1KO)
#'
#' @examples
#' data(Snail1KO)
#' dim(Snail1KO$rawCounts)
#' Snail1KO$metadata
"Snail1KO"
