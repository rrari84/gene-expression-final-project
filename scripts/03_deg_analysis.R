#!/usr/bin/env Rscript

# 03_deg_analysis.R
# -----------------------------------
# Performs DEG analysis using clusterProfiler.
#
# Input:
#   data/processed/<GSE_ID>_expr.rds
#   data/processed/<GSE_ID>_pheno.rds
#
# Output:
#   results/deg/<GSE_ID>_deg_results.csv
#      columns: gene, log2FC, pvalue, padj, dataset, direction
#
# Usage:
#   Rscript scripts/03_deg_analysis.R


message("=== 03_deg_analysis.R: Starting DEG Analysis ===")

# ---- 1. Load required packages ----

needed_pkgs <- c("limma", "edgeR") 

  for (pkg in needed_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        sprintf(
          "Package '%s' is not installed, Add it to 00_setup_environment.R and re-run that script.",
          pkg
        )
      )
    }
  }

library(limma)
library (edgeR)

# ---- 2. Define input/output paths ----

# Make sure in directory exists
rds_dir <- file.path("data", "processed")

if (!dir.exists(rds_dir)) {
  stop(sprintf("RDS directory not found: %s", rds_dir))
}

# Make sure out directory exists
deg_dir <- file.path("results", "deg")

if (!dir.exists(deg_dir)) {
  message(sprintf("[info] Output directory not found, creating now at %s/%s", getwd(), deg_dir))
  dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
}

# ---- 3. List expr.rds files ----

expr_files <- list.files(rds_dir, pattern = "_expr\\.rds$", full.names = TRUE)

expr_gse_ids <- gsub("_expr\\.rds$", "", basename(expr_files))


if (length(expr_files) == 0) {
  stop(sprintf("No expression files found in %s", rds_dir))
}

message(sprintf("[info] Found %d expression file(s)", length(expr_files)))

# ---- 4. Check for matching phenotype files ----
pheno_files <- file.path(rds_dir, paste0(expr_gse_ids, "_pheno.rds"))
pheno_exists <- file.exists(pheno_files)

if (!all(pheno_exists)) {
  missing_ids <- expr_gse_ids[!pheno_exists]
  message(sprintf("[warning] Missing phenotype files for: %s", 
                  paste(missing_ids, collapse = ", ")))
  
  # Keep only datasets with both files
  expr_gse_ids <- expr_gse_ids[pheno_exists]
  expr_files <- expr_files[pheno_exists]
  pheno_files <- pheno_files[pheno_exists]
}

if (length(expr_gse_ids) == 0) {
  stop("No datasets with matching expression and phenotype files found.")
}

message(sprintf("[info] Found %d dataset(s) with matching expr and pheno files", 
                length(expr_gse_ids)))

# ---- 5. Create Analysis function to run across all expression matrices ----
deg_analysis <- function(gse_id, expr_file, pheno_file, out_dir) {
  
  message(sprintf("Processing %s", gse_id))
  
  out_file <- file.path(out_dir, sprintf("%s_deg_results.csv", gse_id))
  
  expr_mat <- readRDS(expr_file)
  
  pheno <- readRDS(pheno_file)
  
  expr_pheno_matches <- intersect(colnames(expr_mat), rownames(pheno))
  if (length(expr_pheno_matches) == 0) {
    stop("No matching samples between expression and phenotype data")
  }
  
  expr_mat <- expr_mat[, expr_pheno_matches]
  pheno <- pheno[expr_pheno_matches, , drop = FALSE]
  
  if (!"group" %in% colnames(pheno)) {
    stop(sprintf("'group' column not found in phenotype data. Available columns: %s", colnames(pheno)))
  }
  
  # check group distribution
  group_dist <- table(pheno$group)
  message("[info] Group distribution:")
  print(group_dist)
  
  if (length(unique(pheno$group)) < 2) {
    stop("Need at least 2 groups for comparison")
  }
  # Create DGEList object
  dge <- DGEList(counts = expr_mat)
  
  keep <- filterByExpr(dge, group = pheno$group)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Normalize DGE
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Create design matrix 
  design <- model.matrix(~ 0 + group, data = pheno)
  colnames(design) <- gsub("^group", "", colnames(design))
  message(sprintf("[info] Design matrix columns: %s", paste0(colnames(design), collapse = ", ")))
  
  # Voom prep for limma
  v <- voom(dge, design, plot = FALSE)
  
  # Fitting linear model
  fit <- lmFit(v, design)
  
  # Create contrast and groupings
  group_names <- colnames(design)
  
  # Group column names into groups by variety of potential indicators
  cancer_group <- grep("cancer|tumor|tumour|malignant", group_names, ignore.case = TRUE, value = TRUE)
  
  normal_group <- grep("normal|control|healthy|benign", group_names, ignore.case = TRUE, value = TRUE)
  
  if (length(cancer_group) > 0 && length(normal_group) > 0) {
    cancer_group <- cancer_group[1]
    normal_group <- normal_group[1]
    contrast_str <- paste0(cancer_group, "-", normal_group)
    message(sprintf("[info] Auto-detected contrast: %s vs %s", cancer_group, normal_group))
  } else {
    # If grouping fails, default to contrast first and second group
    
    contrast_str <- paste0(group_names[1], "-", group_names[2])
    message(sprintf("[warning] Could not auto-group cancer vs normal groups"))
    message(sprintf("[info] Using contrast: %s vs %s", group_names[1], group_names[2]))
  }
  
  contrast_mat <- makeContrasts(contrasts = contrast_str, levels = design)
  
  # Fit contrast matrix
  fit2 <- contrasts.fit(fit, contrast_mat)
  fit2 <- eBayes(fit2)
  
  # Extract results
  results <- topTable(fit2, number = Inf, adjust.method = "BH")
  
  # Format output
  results_cleaned <- data.frame(
    gene = rownames(results),
    log2FC = results$logFC,
    pvalue = results$P.Value,
    padj = results$adj.P.Val,
    dataset = gse_id,
    stringsAsFactors = FALSE
  )
  
  # Add upreg vs downreg
  results_cleaned$direction <- "None"
  results_cleaned$direction[results_cleaned$padj < 0.05 & results_cleaned$log2FC > 0] <- "Up"
  results_cleaned$direction[results_cleaned$padj < 0.05 & results_cleaned$log2FC < 0] <- "Down"
  
  # Sort matrix by padj
  results_cleaned <- results_cleaned[order(results_cleaned$padj), ]
  
  # Summary
  n_sig <- sum(results_cleaned$padj < 0.05, na.rm = TRUE)
  n_up <- sum(results_cleaned$direction == "Up", na.rm = TRUE)
  n_down <- sum(results_cleaned$direction == "Down", na.rm = TRUE)
  
  message(sprintf("[info] Results summary:"))
  message(sprintf("  Total genes tested: %d", nrow(results_cleaned)))
  message(sprintf("  Significant (padj < 0.05): %d", n_sig))
  message(sprintf("  Upregulated: %d", n_up))
  message(sprintf("  Downregulated: %d", n_down))
  
  # Save
  write.csv(results_cleaned, out_file, row.names = FALSE)
  message(sprintf("[info] Saved results to: %s", out_file))
  
  return(results_cleaned)
}

# ---- 6. Apply the analysis function to datasets ----

n_success <- 0
n_total <- length(expr_gse_ids)

for (i in seq_along(expr_gse_ids)) {
  gse_id <- expr_gse_ids[i]
  expr_file <- expr_files[i]
  pheno_file <- pheno_files[i]
  
  result <- tryCatch({
    deg_analysis(gse_id, expr_file, pheno_file, deg_dir)
  }, error = function(e) {
    message(sprintf("[error] Failed to process %s: %s", gse_id, conditionMessage(e)))
    return(NULL)
  })
  
  if (!is.null(result)) {
    n_success <- n_success + 1
  }
  
}

message("\n=== DEG Analysis Complete ===")
message(sprintf("[info] Successfully processed: %d/%d datasets", n_success, n_total))

if (n_success < n_total) {
  message("[warning] Some datasets failed. Check error messages above.")
}