#!/usr/bin/env Rscript

# 02_download_data.R
# -------------------
# Download GEO series data and save expression + phenotype tables.
#
# This is meant for datasets where we want to pull the *processed* series
# matrix from GEO (not raw CEL files).
#
# It will:
#   - Read a list of GEO IDs from command-line args OR use a default vector.
#   - For each GSE:
#       * Use GEOquery::getGEO() to download the series matrix.
#       * Extract expression matrix and phenotype (sample metadata).
#       * Save them as:
#             data/processed/<GSE_ID>_expr.rds
#             data/processed/<GSE_ID>_pheno.rds
#
# Usage (from project root, with conda env active):
#   Rscript scripts/02_download_data.R GSE12345 GSE67890
#
# If you run without arguments, it will use a small default list defined below.

message("=== 02_download_data.R: Starting GEO download ===")

## 1. Parse command-line arguments -----------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # TODO: replace with the actual GEO IDs you want as defaults
  GSE_IDS <- c("GSE_PLACEHOLDER1", "GSE_PLACEHOLDER2")
  message("[info] No GEO IDs supplied on the command line.")
  message("[info] Using default GEO IDs: ", paste(GSE_IDS, collapse = ", "))
} else {
  GSE_IDS <- args
  message("[info] GEO IDs from command line: ", paste(GSE_IDS, collapse = ", "))
}

## 2. Load required package(s) ---------------------------------------------

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop(
    "Package 'GEOquery' is not installed in this environment.\n",
    "Please install it via conda, e.g.:\n",
    "  conda install bioconductor-geoquery -c conda-forge -c bioconda\n"
  )
}

library(GEOquery)

## 3. Ensure output directories exist --------------------------------------

raw_dir       <- file.path("data", "raw")
processed_dir <- file.path("data", "processed")

if (!dir.exists(raw_dir)) {
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(processed_dir)) {
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
}

## 4. Helper to process a single GSE ---------------------------------------

process_one_gse <- function(gse_id) {
  message("\n--- Processing GEO series: ", gse_id, " ---")

  # Directory where GEOquery will store downloaded files
  gse_raw_dir <- file.path(raw_dir, gse_id)
  if (!dir.exists(gse_raw_dir)) {
    dir.create(gse_raw_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Output paths for processed data
  expr_out_path  <- file.path(processed_dir, sprintf("%s_expr.rds",  gse_id))
  pheno_out_path <- file.path(processed_dir, sprintf("%s_pheno.rds", gse_id))

  # If outputs already exist, you may choose to skip
  if (file.exists(expr_out_path) && file.exists(pheno_out_path)) {
    message("[info] Processed files already exist for ", gse_id, "; skipping.")
    return(invisible(NULL))
  }

  message("[step] Downloading series matrix using GEOquery::getGEO() ...")

  # getGEO usually returns a list if multiple platforms are present
  gse <- NULL
  try({
    gse <- GEOquery::getGEO(
      GSEMatrix = TRUE,
      GEO       = gse_id,
      destdir   = gse_raw_dir
    )
  }, silent = TRUE)

  if (is.null(gse)) {
    warning("[warn] getGEO() returned NULL for ", gse_id, "; skipping.")
    return(invisible(NULL))
  }

  # If multiple platforms, pick the first; you can refine this later
  if (is.list(gse)) {
    message("[info] getGEO() returned ", length(gse), " ExpressionSet objects; using the first one.")
    eset <- gse[[1]]
  } else {
    eset <- gse
  }

  # Extract expression matrix and phenotype
  expr_mat <- Biobase::exprs(eset)
  pheno    <- Biobase::pData(eset)

  message("[info] Expression matrix dimensions (features x samples): ",
          nrow(expr_mat), " x ", ncol(expr_mat))

  # Basic sanity check on sample names
  if (!identical(colnames(expr_mat), rownames(pheno))) {
    warning(
      "[warn] Sample names differ between expression and pheno data for ",
      gse_id, ".\n",
      "  - expr colnames (first few): ", paste(head(colnames(expr_mat)), collapse = ", "), "\n",
      "  - pheno rownames  (first few): ", paste(head(rownames(pheno)), collapse = ", "), "\n",
      "Downstream scripts should carefully match samples by ID."
    )
  } else {
    message("[info] Sample names match between expression and phenotype data.")
  }

  # Save processed files
  message("[step] Saving processed expression and phenotype data for ", gse_id, " ...")
  saveRDS(expr_mat,  file = expr_out_path)
  saveRDS(pheno,     file = pheno_out_path)

  message("[done] Finished processing ", gse_id,
          ". Outputs:\n  - ", expr_out_path, "\n  - ", pheno_out_path)
}

## 5. Loop over all GSE IDs -------------------------------------------------

for (gse_id in GSE_IDS) {
  process_one_gse(gse_id)
}

message("\n=== 02_download_data.R: Finished ===")
