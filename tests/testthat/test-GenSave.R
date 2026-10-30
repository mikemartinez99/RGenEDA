library(testthat)
library(RGenEDA)

# Helper: an RGenEDA-style heatmap list with a silent pheatmap
create_test_heatmap <- function() {
    set.seed(1)
    mat <- matrix(rnorm(100), 10, 10)
    list(heatmap = pheatmap::pheatmap(mat, silent = TRUE))
}

# Read the first bytes of a file as a hex string
file_magic <- function(path, n = 4) {
    paste(as.character(readBin(path, "raw", n)), collapse = "")
}

#----- Each extension gets its own file format
test_that("GenSave writes the file format that matches the extension", {
    hm <- create_test_heatmap()
    dir <- withr::local_tempdir()

    expected <- c(
        png = "89504e47",
        jpg = "ffd8ff",
        jpeg = "ffd8ff",
        tiff = "^(49492a00|4d4d002a)$",
        bmp = "424d",
        pdf = "25504446"
    )
    for (ext in names(expected)) {
        f <- file.path(dir, paste0("hm.", ext))
        GenSave(hm, f, width = 2, height = 2)
        expect_true(file.exists(f), info = ext)
        expect_match(file_magic(f), expected[[ext]], info = ext)
    }
})

#----- Vector output converts units to inches
test_that("GenSave converts cm to inches for PDF output", {
    hm <- create_test_heatmap()
    f <- file.path(withr::local_tempdir(), "hm.pdf")
    GenSave(hm, f, width = 2.54 * 4, height = 2.54 * 3, units = "cm")

    # PDF MediaBox is in points (72 per inch): 4 x 3 inches = 288 x 216
    # The pdf device truncates to whole points, so allow 1 point of rounding
    # Search the raw bytes: the binary PDF header is invalid UTF-8, so string
    # regex on rawToChar() output finds no match in UTF-8 locales
    pdf_raw <- readBin(f, "raw", file.size(f))
    box <- rawToChar(
        grepRaw("MediaBox \\[0 0 [0-9.]+ [0-9.]+", pdf_raw, value = TRUE)
    )
    dims <- as.numeric(strsplit(sub("MediaBox \\[0 0 ", "", box), " ")[[1]])
    expect_true(all(abs(dims - c(288, 216)) <= 1))
})

#----- Unsupported extensions error
test_that("GenSave errors on unsupported file formats", {
    hm <- create_test_heatmap()
    expect_error(
        GenSave(hm, file.path(tempdir(), "hm.xyz")),
        "Unsupported file format: xyz"
    )
})
