#!/usr/bin/env Rscript

# 03_deg_analysis.R
# -----------------------------------
# Performs per-dataset DEG analysis with limma on
# preprocessed GEO expression matrices.
#
# Input (from 02_download_data.R):
#   data/processed/<GSE_ID>_expr.rds  # normalized expression (genes x samples)
#   data/processed/<GSE_ID>_pheno.rds # sample metadata; must contain 'group'
#
# Output:
#   results/deg/<GSE_ID>_deg_results.csv
#      columns: gene, log2FC, pvalue, padj, dataset, direction
#
# Usage:
#   Rscript scripts/03_deg_analysis.R
#
# Notes:
#   - This script assumes expression values are already on a log scale
#     (which 02_download_data.R enforces).
#   - 'group' in the phenotype table should encode the comparison groups
#     (e.g. "cancer" vs "normal").

message("=== 03_deg_analysis.R: Starting DEG Analysis ===")

# ---- 1. Load required packages ----

needed_pkgs <- c("limma")

for (pkg in needed_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf(
        "Package '%s' is not installed. Install it in the conda env, e.g.:\n  conda install r-%s -c conda-forge",
        pkg, pkg
      )
    )
  }
}

library(limma)

# ---- 2. Define input/output paths ----

rds_dir <- file.path("data", "processed")

if (!dir.exists(rds_dir)) {
  stop(sprintf("Processed RDS directory not found: %s", rds_dir))
}

deg_dir <- file.path("results", "deg")

if (!dir.exists(deg_dir)) {
  message(sprintf("[info] Output directory not found, creating now at %s/%s",
                  getwd(), deg_dir))
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

pheno_files   <- file.path(rds_dir, paste0(expr_gse_ids, "_pheno.rds"))
pheno_exists  <- file.exists(pheno_files)

if (!all(pheno_exists)) {
  missing_ids <- expr_gse_ids[!pheno_exists]
  message(sprintf("[warning] Missing phenotype files for: %s",
                  paste(missing_ids, collapse = ", ")))
  # Keep only datasets with both files
  expr_gse_ids <- expr_gse_ids[pheno_exists]
  expr_files   <- expr_files[pheno_exists]
  pheno_files  <- pheno_files[pheno_exists]
}

if (length(expr_gse_ids) == 0) {
  stop("No datasets with matching expression and phenotype files found.")
}

message(sprintf("[info] Found %d dataset(s) with matching expr and pheno files",
                length(expr_gse_ids)))

# ---- 5. DEG analysis function for one dataset ----

deg_analysis <- function(gse_id, expr_file, pheno_file, out_dir) {
  message(sprintf("\n--- Processing %s ---", gse_id))

  out_file <- file.path(out_dir, sprintf("%s_deg_results.csv", gse_id))

  expr_mat <- readRDS(expr_file)
  pheno    <- readRDS(pheno_file)

  # If someone accidentally saved an ExpressionSet, handle that:
  if (inherits(expr_mat, "ExpressionSet")) {
    expr_mat <- Biobase::exprs(expr_mat)
  }

  # Ensure it's a numeric matrix
  expr_mat <- as.matrix(expr_mat)

  # Match samples between expression and phenotype
  expr_pheno_matches <- intersect(colnames(expr_mat), rownames(pheno))
  if (length(expr_pheno_matches) == 0) {
    stop("No matching samples between expression and phenotype data")
  }

  expr_mat <- expr_mat[, expr_pheno_matches, drop = FALSE]
  pheno    <- pheno[expr_pheno_matches, , drop = FALSE]

  # ---- 5a. Ensure pheno$group exists (dataset-specific mapping) ----
  if (!"group" %in% colnames(pheno)) {
    message("[info] 'group' column not found; attempting to derive it.")
    cols <- colnames(pheno)
    message("[info] Available columns: ", paste(cols, collapse = ", "))

    if (gse_id == "GSE26910" && "normal/tumor:ch1" %in% cols) {
      col_nm <- "normal/tumor:ch1"
      vals <- pheno[[col_nm]]
      message("[info] Using column '", col_nm, "' for group (normal vs tumor)")
      print(table(vals))

      pheno$group <- ifelse(
        grepl("tumor", vals, ignore.case = TRUE),
        "cancer",
        "normal"
      )

    } else if (gse_id == "GSE54002" && "tumor/nontumor:ch1" %in% cols) {
      col_nm <- "tumor/nontumor:ch1"
      vals <- pheno[[col_nm]]
      message("[info] Using column '", col_nm, "' for group (tumor vs nontumor)")
      print(table(vals))

      pheno$group <- ifelse(
        grepl("tumor", vals, ignore.case = TRUE),
        "cancer",
        "normal"
      )

    } else if (gse_id == "GSE65194" && "sample_group:ch1" %in% cols) {
      col_nm <- "sample_group:ch1"
      vals <- pheno[[col_nm]]
      message("[info] Using column '", col_nm, "' for group from sample_group:ch1")
      print(table(vals))

      pheno$group <- dplyr::case_when(
        grepl("normal",  vals, ignore.case = TRUE) ~ "normal",
        grepl("control", vals, ignore.case = TRUE) ~ "normal",
        grepl("tumor|tumour|cancer|carcinoma", vals, ignore.case = TRUE) ~ "cancer",
        TRUE ~ "other"
      )

    } else if (gse_id == "GSE42568" && "tissue:ch1" %in% cols) {
      col_nm <- "tissue:ch1"
      vals <- pheno[[col_nm]]
      message("[info] Using column '", col_nm, "' for group (tumor vs normal tissue)")
      print(table(vals))

      pheno$group <- dplyr::case_when(
        grepl("normal",  vals, ignore.case = TRUE) ~ "normal",
        grepl("control", vals, ignore.case = TRUE) ~ "normal",
        grepl("tumor|tumour|cancer|carcinoma", vals, ignore.case = TRUE) ~ "cancer",
        TRUE ~ "other"
      )

    } else {
      stop(sprintf(
        "'group' column not found and no mapping implemented for %s.\nAvailable columns: %s",
        gse_id,
        paste(cols, collapse = ", ")
      ))
    }
  }

  # Check group distribution
  group_dist <- table(pheno$group)
  message("[info] Group distribution (pheno$group):")
  print(group_dist)

  if (length(unique(pheno$group)) < 2) {
    stop("Need at least 2 groups in 'group' for comparison")
  }

  # ---- limma design and model ----

  # Design matrix with no intercept: one column per group level
  design <- model.matrix(~ 0 + group, data = pheno)
  colnames(design) <- gsub("^group", "", colnames(design))

  group_names <- colnames(design)
  message(sprintf("[info] Design matrix columns: %s",
                  paste(group_names, collapse = ", ")))

  # Try to auto-detect cancer vs normal
  cancer_group <- grep("cancer|tumor|tumour|malignant", group_names,
                       ignore.case = TRUE, value = TRUE)
  normal_group <- grep("normal|control|healthy|benign", group_names,
                       ignore.case = TRUE, value = TRUE)

  if (length(cancer_group) > 0 && length(normal_group) > 0) {
    cancer_group <- cancer_group[1]
    normal_group <- normal_group[1]
    contrast_str <- paste0(cancer_group, " - ", normal_group)
    message(sprintf("[info] Auto-detected contrast: %s vs %s",
                    cancer_group, normal_group))
  } else {
    # If we can't auto-detect, fall back to the first two groups
    if (length(group_names) < 2) {
      stop("Design matrix has fewer than 2 group columns; cannot build contrast.")
    }
    contrast_str <- paste0(group_names[1], " - ", group_names[2])
    message("[warning] Could not auto-group cancer vs normal.")
    message(sprintf("[info] Using contrast: %s vs %s",
                    group_names[1], group_names[2]))
  }

  # Fit model
  fit <- lmFit(expr_mat, design)

  contrast_mat <- makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- contrasts.fit(fit, contrast_mat)

  # Empirical Bayes
  fit2 <- eBayes(fit2)

  # Extract all genes
  results <- topTable(fit2, number = Inf, adjust.method = "BH")

  # Format output
  results_cleaned <- data.frame(
    gene    = rownames(results),
    log2FC  = results$logFC,
    pvalue  = results$P.Value,
    padj    = results$adj.P.Val,
    dataset = gse_id,
    stringsAsFactors = FALSE
  )

  # Label direction
  results_cleaned$direction <- "None"
  results_cleaned$direction[
    results_cleaned$padj < 0.05 & results_cleaned$log2FC > 0
  ] <- "Up"
  results_cleaned$direction[
    results_cleaned$padj < 0.05 & results_cleaned$log2FC < 0
  ] <- "Down"

  # Sort by adjusted p-value
  results_cleaned <- results_cleaned[order(results_cleaned$padj), ]

  # Summary
  n_sig  <- sum(results_cleaned$padj < 0.05, na.rm = TRUE)
  n_up   <- sum(results_cleaned$direction == "Up",   na.rm = TRUE)
  n_down <- sum(results_cleaned$direction == "Down", na.rm = TRUE)

  message("[info] Results summary:")
  message(sprintf("  Total genes tested: %d", nrow(results_cleaned)))
  message(sprintf("  Significant (padj < 0.05): %d", n_sig))
  message(sprintf("  Upregulated: %d", n_up))
  message(sprintf("  Downregulated: %d", n_down))

  # Save to CSV
  write.csv(results_cleaned, out_file, row.names = FALSE)
  message(sprintf("[info] Saved results to: %s", out_file))

  invisible(results_cleaned)
}


# ---- 6. Apply analysis to all datasets ----

n_success <- 0
n_total   <- length(expr_gse_ids)

for (i in seq_along(expr_gse_ids)) {
  gse_id     <- expr_gse_ids[i]
  expr_file  <- expr_files[i]
  pheno_file <- pheno_files[i]

  result <- tryCatch(
    {
      deg_analysis(gse_id, expr_file, pheno_file, deg_dir)
    },
    error = function(e) {
      message(sprintf("[error] Failed to process %s: %s",
                      gse_id, conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(result)) {
    n_success <- n_success + 1
  }
}

message("\n=== DEG Analysis Complete ===")
message(sprintf("[info] Successfully processed: %d/%d datasets",
                n_success, n_total))

if (n_success < n_total) {
  message("[warning] Some datasets failed. Check error messages above.")
}
