# RGenEDA Development List

- During `GenEDA` object construction, allow counts in dataframe format, convert to matrix internally
- Change variable names for previous functions to lower case
- Plots to ggplots
- PCA plotting script (arguments should be PCA df, x, y, color_by, shape_by, split_by, return_data which would call ExtractPCA under the hood if TRUE)
- Would be cool to have a differential expression slot that takes your DEGs from DESeq2, adds labels based on numerator/denominator and then you can directly plot volcano, and MA