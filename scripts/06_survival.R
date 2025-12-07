#!/usr/bin/env Rscript

# 06_survival.R
# -----------------------------------
# Survival analysis for hub genes using TCGA BRCA expression + clinical data.
#
# Inputs:
#   data/processed/hub_genes.csv         (column: gene = HGNC SYMBOL)
#   data/processed/tcga_brca_expr.rds    (TCGA expression, saved once)
#   data/processed/tcga_brca_clinical.rds (TCGA clinical, saved once)
#
# Outputs:
#   results/survival/survival_<GENE>.png
#
# Usage:
#   Rscript scripts/06_survival.R

message("=== 06_survival.R: Starting Survival Analysis ===")

# ---- 1. Load packages ----
pkgs <- c("TCGAbiolinks", "dplyr", "survival", "survminer",
          "SummarizedExperiment", "AnnotationDbi", "org.Hs.eg.db", "readr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' not installed. Install it first.", p))
  }
}

library(TCGAbiolinks)
library(dplyr)
library(survival)
library(survminer)
library(SummarizedExperiment)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readr)

# ---- 2. Load hub genes (SYMBOLS) ----
hub_path <- file.path("data", "processed", "hub_genes.csv")
if (!file.exists(hub_path)) stop("Hub genes file not found: ", hub_path)

hub_genes <- read_csv(hub_path, show_col_types = FALSE)

if (!"gene" %in% colnames(hub_genes)) {
  stop("hub_genes.csv must contain a column named 'gene' with HGNC symbols.")
}

genes <- unique(hub_genes$gene)
message("[info] Hub genes (symbols): ", paste(genes, collapse = ", "))

# ---- 3. Load or download TCGA BRCA expression + clinical ----

expr_rds  <- file.path("data", "processed", "tcga_brca_expr.rds")
clin_rds  <- file.path("data", "processed", "tcga_brca_clinical.rds")

if (file.exists(expr_rds) && file.exists(clin_rds)) {
  message("[info] Loading pre-saved TCGA BRCA data from RDS files ...")
  exp_data <- readRDS(expr_rds)
  clinical <- readRDS(clin_rds)
} else {
  message("[step] Downloading TCGA BRCA expression (HTSeq - FPKM) ...")

  query.exp <- GDCquery(
    project       = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "HTSeq - FPKM"
  )

  GDCdownload(query.exp)
  exp_data <- GDCprepare(query.exp)

  message("[step] Downloading TCGA BRCA clinical data ...")
  clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

  # Save for future runs
  dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
  saveRDS(exp_data, expr_rds)
  saveRDS(clinical, clin_rds)
  message("[info] Saved TCGA BRCA expression + clinical to data/processed/")
}

# ---- 4. Prepare expression at SYMBOL level ----

message("[step] Preparing expression matrix at gene symbol level ...")

# assay(exp_data) is usually genes (rows) x samples (columns), rownames = Ensembl IDs
exp_mat <- as.data.frame(SummarizedExperiment::assay(exp_data))
exp_mat$ensembl <- rownames(exp_mat)

# strip version from Ensembl IDs, e.g. ENSG00000141510.15 -> ENSG00000141510
exp_mat$ensembl_clean <- sub("\\..*$", "", exp_mat$ensembl)

# Map Ensembl -> SYMBOL using org.Hs.eg.db
anno_tcga <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(exp_mat$ensembl_clean),
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)

anno_tcga <- anno_tcga %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::distinct(ENSEMBL, SYMBOL)

# Merge annotation back into expression
exp_mat <- exp_mat %>%
  dplyr::inner_join(anno_tcga, by = c("ensembl_clean" = "ENSEMBL"))

# Now we have multiple rows per SYMBOL sometimes; we summarise by mean expression
sample_cols <- setdiff(colnames(exp_mat), c("ensembl", "ensembl_clean", "SYMBOL"))

message("[info] Collapsing expression to SYMBOL level by averaging Ensembl rows ...")
exp_by_symbol <- exp_mat %>%
  dplyr::select(SYMBOL, dplyr::all_of(sample_cols)) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(sample_cols), mean, na.rm = TRUE), .groups = "drop")

# Ensure SYMBOL as rownames for convenience
exp_mat_sym <- as.data.frame(exp_by_symbol)
rownames(exp_mat_sym) <- exp_mat_sym$SYMBOL
exp_mat_sym$SYMBOL <- NULL

# ---- 5. Prepare clinical OS variables + align IDs ----

message("[step] Preparing clinical overall survival variables ...")

clinical_os <- clinical %>%
  dplyr::mutate(
    # Try different possible names for last follow-up
    days_last_followup = dplyr::coalesce(
      .data$days_to_last_follow_up,
      .data$days_to_last_followup
    ),
    os_time = dplyr::coalesce(.data$days_to_death, days_last_followup),
    os_event = ifelse(.data$vital_status == "Dead", 1, 0)
  ) %>%
  dplyr::filter(!is.na(os_time))

# We assume exp_mat_sym columns are sample barcodes like TCGA-XX-XXXX-01A-...
sample_barcodes <- colnames(exp_mat_sym)
patient_ids <- substr(sample_barcodes, 1, 12)

expr_long <- data.frame(
  sample_barcode = sample_barcodes,
  patient_id     = patient_ids,
  stringsAsFactors = FALSE
)

# ---- 6. Survival analysis for each hub gene ----

surv_dir <- file.path("results", "survival")
dir.create(surv_dir, showWarnings = FALSE, recursive = TRUE)

for (g in genes) {
  message("\n[info] Running survival analysis for gene: ", g)

  if (!(g %in% rownames(exp_mat_sym))) {
    warning(sprintf("Gene %s not found in TCGA expression (after mapping to SYMBOL).", g))
    next
  }

  # Expression vector for that gene across samples
  gene_expr <- as.numeric(exp_mat_sym[g, ])
  names(gene_expr) <- sample_barcodes

  df_expr <- data.frame(
    sample_barcode = names(gene_expr),
    expr           = gene_expr,
    stringsAsFactors = FALSE
  ) %>%
    left_join(expr_long, by = "sample_barcode") %>%
    dplyr::filter(!is.na(patient_id))

  # Merge expression with clinical OS
  merged <- df_expr %>%
    dplyr::inner_join(
      clinical_os,
      by = c("patient_id" = "submitter_id")
    )

  if (nrow(merged) < 20) {
    warning(sprintf("Not enough patients with data for gene %s; skipping.", g))
    next
  }

  median_expr <- median(merged$expr, na.rm = TRUE)
  merged$group <- ifelse(merged$expr > median_expr, "High", "Low")

  s <- Surv(time = merged$os_time, event = merged$os_event)

  fit <- survfit(s ~ group, data = merged)

  p <- ggsurvplot(
    fit,
    data        = merged,
    pval        = TRUE,
    risk.table  = TRUE,
    title       = paste("Overall Survival -", g, "expression"),
    legend.title = "Expression",
    legend.labs  = c("High", "Low")
  )

  out_png <- file.path(surv_dir, sprintf("survival_%s.png", g))

  ggsave(
    filename = out_png,
    plot     = p,
    width    = 8,
    height   = 6
  )

  message("[info] Saved survival plot to: ", out_png)
}

message("\n=== 06_survival.R: Completed Successfully ===")
