#!/usr/bin/env Rscript

# 02_download_data.R
# -------------------
# Preprocess GEO microarray datasets by using the provided "series matrix"
# as our normalized expression data.
#
# For each GSE ID:
#   - Download series matrix with GEOquery::getGEO()
#   - Extract expression matrix and phenotype (sample metadata)
#   - Do a light check on scale (log2-ish vs raw) and log2-transform if needed
#   - Save to:
#       data/processed/<GSE>_expr.rds
#       data/processed/<GSE>_pheno.rds
#
# Usage (from project root, with conda env active):
#   Rscript scripts/02_download_data.R GSE26910 GSE42568 GSE65194
#
# If no GSE IDs are supplied as arguments, it falls back to a default list.

message("=== 02_download_data.R: Starting GEO preprocessing ===")

## 1. Parse command-line arguments -----------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # You can hard-code your main microarray GSEs here
  GSE_IDS <- c("GSE26910", "GSE42568", "GSE65194")
  message("[info] No GEO IDs supplied on the command line.")
  message("[info] Using default GEO IDs: ", paste(GSE_IDS, collapse = ", "))
} else {
  GSE_IDS <- args
  message("[info] GEO IDs from command line: ", paste(GSE_IDS, collapse = ", "))
}

## 2. Load required packages -----------------------------------------------

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop(
    "Package 'GEOquery' is not installed in this environment.\n",
    "Install via conda:\n",
    "  conda install bioconductor-geoquery -c conda-forge -c bioconda\n"
  )
}
if (!requireNamespace("Biobase", quietly = TRUE)) {
  stop("Package 'Biobase' is not installed in this environment.")
}

library(GEOquery)
library(Biobase)

## 3. Ensure output directories exist --------------------------------------

raw_dir       <- file.path("data", "raw")
processed_dir <- file.path("data", "processed")

if (!dir.exists(raw_dir)) {
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(processed_dir)) {
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
}

## 4. Helper: light "is this log-ish?" check -------------------------------

is_log_like <- function(x) {
  # drop NAs
  x <- x[is.finite(x)]
  if (length(x) == 0) return(TRUE)
  # heuristic: if most values are between ~0 and ~20, assume log-like
  max_val <- quantile(x, 0.99, na.rm = TRUE)
  return(max_val < 100)  # arbitrary but works for separating raw vs log2
}

## 5. Process one GSE ------------------------------------------------------

process_one_gse <- function(gse_id) {
  message("\n--- Processing GEO series: ", gse_id, " ---")

  # Where to store raw GEO files (optional)
  gse_raw_dir <- file.path(raw_dir, gse_id)
  if (!dir.exists(gse_raw_dir)) {
    dir.create(gse_raw_dir, recursive = TRUE, showWarnings = FALSE)
  }

  expr_out_path  <- file.path(processed_dir, sprintf("%s_expr.rds",  gse_id))
  pheno_out_path <- file.path(processed_dir, sprintf("%s_pheno.rds", gse_id))

  # If we already processed this GSE, you can choose to skip
  if (file.exists(expr_out_path) && file.exists(pheno_out_path)) {
    message("[info] Processed files already exist for ", gse_id, "; skipping.")
    return(invisible(NULL))
  }

  message("[step] Downloading series matrix with GEOquery::getGEO() ...")

  gse <- NULL
  try({
    gse <- GEOquery::getGEO(
      GEO       = gse_id,
      GSEMatrix = TRUE,
      destdir   = gse_raw_dir
    )
  }, silent = TRUE)

  if (is.null(gse)) {
    warning("[warn] getGEO() returned NULL for ", gse_id, "; skipping.")
    return(invisible(NULL))
  }

  # If multiple platforms, pick the first for now (can refine later)
  if (is.list(gse)) {
    message("[info] getGEO() returned ", length(gse), " ExpressionSet objects; using the first one.")
    eset <- gse[[1]]
  } else {
    eset <- gse
  }

  expr_mat <- Biobase::exprs(eset)
  pheno    <- Biobase::pData(eset)

  message("[info] Expression matrix dimensions (features x samples): ",
          nrow(expr_mat), " x ", ncol(expr_mat))

  # Light check: log vs linear
  if (!is_log_like(expr_mat)) {
    message("[info] Expression values look non-log; applying log2(x + 1) transform.")
    expr_mat <- log2(expr_mat + 1)
  } else {
    message("[info] Expression values already look log-like; leaving as-is.")
  }

  # Sample name sanity check
  if (!identical(colnames(expr_mat), rownames(pheno))) {
    warning(
      "[warn] Sample names differ between expression and phenotype data for ",
      gse_id, ".\n",
      "  - expr colnames (first few): ", paste(head(colnames(expr_mat)), collapse = ", "), "\n",
      "  - pheno rownames  (first few): ", paste(head(rownames(pheno)), collapse = ", "), "\n",
      "Downstream scripts should match samples carefully."
    )
  } else {
    message("[info] Sample names match between expression and phenotype data.")
  }

  # Save processed objects
  message("[step] Saving processed expression and phenotype data for ", gse_id, " ...")
  saveRDS(expr_mat,  file = expr_out_path)
  saveRDS(pheno,     file = pheno_out_path)

  message("[done] Finished processing ", gse_id,
          ". Outputs:\n  - ", expr_out_path, "\n  - ", pheno_out_path)
}

## 6. Loop over all GSE IDs ------------------------------------------------

for (gse_id in GSE_IDS) {
  process_one_gse(gse_id)
}

message("\n=== 02_download_data.R: Finished GEO preprocessing ===")
