#!/usr/bin/env Rscript

# 00_preprocess_microarray_from_CEL.R
# -----------------------------------
# Preprocess raw microarray .CEL files into a normalized expression matrix.
#
# - Reads all .CEL files from data/raw/<DATASET_ID>/CEL/
# - Runs RMA normalization (background correction + normalization + summarization)
# - Saves:
#     data/processed/<DATASET_ID>_expr.rds   # expression matrix (probes x samples)
#     data/processed/<DATASET_ID>_pheno.rds  # sample metadata
#
# Usage (from project root, after loading R/4.4.1 and activating renv):
#   Rscript scripts/01_preprocess_microarray_from_CEL.R GSE12345
#
# If you don't pass an argument, it will use DATASET_ID <- "GSE_PLACEHOLDER"
# (change that below).

message("=== 01_preprocess_microarray_from_CEL.R: Starting ===")

# ---- 0. Parse command-line argument (dataset ID) ----

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  DATASET_ID <- args[1]
} else {
  # TODO: change this to a real default dataset ID if you want
  DATASET_ID <- "GSE_PLACEHOLDER"
  message(sprintf(
    "[info] No dataset ID supplied as argument. Using default: %s",
    DATASET_ID
  ))
}

message(sprintf("[info] Processing dataset: %s", DATASET_ID))

# ---- 1. Load required packages ----

needed_pkgs <- c("oligo", "Biobase")

for (pkg in needed_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf(
        "Package '%s' is not installed. Add it to 00_setup_environment.R and re-run that script.",
        pkg
      )
    )
  }
}

library(oligo)
library(Biobase)

# ---- 2. Define input/output paths ----

cel_dir <- file.path("data", "raw", DATASET_ID, "CEL")
if (!dir.exists(cel_dir)) {
  stop(sprintf("CEL directory not found: %s", cel_dir))
}

# Make sure processed dir exists
processed_dir <- file.path("data", "processed")
if (!dir.exists(processed_dir)) {
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
}

expr_out_path  <- file.path(processed_dir, sprintf("%s_expr.rds",  DATASET_ID))
pheno_out_path <- file.path(processed_dir, sprintf("%s_pheno.rds", DATASET_ID))

message(sprintf("[info] CEL input directory: %s", cel_dir))
message(sprintf("[info] Expression output:   %s", expr_out_path))
message(sprintf("[info] Pheno output:        %s", pheno_out_path))

# ---- 3. List CEL files ----

cel_files <- list.files(cel_dir, pattern = "\\.[cC][eE][lL]$", full.names = TRUE)

if (length(cel_files) == 0) {
  stop(sprintf("No .CEL files found in %s", cel_dir))
}

message(sprintf("[info] Found %d CEL files.", length(cel_files)))

# ---- 4. Read CEL files and run RMA normalization ----

message("[step] Reading CEL files with oligo::read.celfiles() ...")

# NOTE:
#  For Affymetrix arrays, oligo can usually auto-detect the platform if the
#  appropriate 'pd-*' annotation package is installed.
#  If this fails, you may need to install a platform-specific package, e.g.:
#     BiocManager::install('pd.hugene.1.0.st.v1')
#
#  For now we assume the correct platform package has been installed already
#  via 00_setup_environment.R.

raw_data <- read.celfiles(cel_files)

message("[step] Running RMA normalization (rma()) ...")
norm_data <- rma(raw_data)   # ExpressionSet

# ---- 5. Extract expression matrix and phenotype data ----

expr_mat <- Biobase::exprs(norm_data)  # rows = probes/features, cols = samples
pheno    <- Biobase::pData(norm_data)  # sample metadata

message("[info] Expression matrix dimensions (features x samples):")
message(sprintf("       %d x %d", nrow(expr_mat), ncol(expr_mat)))

# Ensure sample names are consistent
sample_names_expr  <- colnames(expr_mat)
sample_names_pheno <- rownames(pheno)

if (!identical(sample_names_expr, sample_names_pheno)) {
  warning(
    "Sample names in expression matrix and phenotype data are not identical.\n",
    "  - Expression colnames: first few = ", paste(head(sample_names_expr), collapse = ", "), "\n",
    "  - Pheno rownames:      first few = ", paste(head(sample_names_pheno), collapse = ", "), "\n",
    "Proceeding, but downstream scripts should be careful when matching samples."
  )
} else {
  message("[info] Sample names match between expression matrix and phenotype data.")
}

# ---- 6. (Optional) Map probes to gene symbols (platform-specific) ----
#
# At this stage, expr_mat rows are PROBE IDs. That is fine for limma,
# but for biological interpretation you often want gene symbols.
#
# This requires a platform-specific annotation package (e.g. hgu133plus2.db,
# hugene10sttranscriptcluster.db, etc.). Because this differs per array,
# we put the mapping step as a clearly-marked TODO.

# Example structure (commented out):
#
# library(AnnotationDbi)
# library(hgu133plus2.db)  # <-- change to your array-specific package
#
# probe_ids <- rownames(expr_mat)
# ann <- AnnotationDbi::select(
#   hgu133plus2.db,
#   keys    = probe_ids,
#   columns = c("SYMBOL", "GENENAME"),
#   keytype = "PROBEID"
# )
#
# # Merge and collapse to one row per gene symbol, e.g. median across probes
# # (Implementation depends on how you want to handle duplicates; we keep it
# # as probe-level for now to keep the script general.)

# ---- 7. Save outputs ----

message("[step] Saving processed expression and phenotype data as .rds files ...")

saveRDS(expr_mat,  file = expr_out_path)
saveRDS(pheno,     file = pheno_out_path)

message("=== 01_preprocess_microarray_from_CEL.R: Finished successfully ===")
