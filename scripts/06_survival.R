#!/usr/bin/env Rscript

# 06_survival.R
# -----------------------------------
# Performs survival analysis using TCGA BRCA clinical + expression data.
#
# Input:
#   data/processed/hub_genes.csv  (must contain column: "gene")
#
# Output:
#   plots/survival_<GENE>.png
#
# Usage:
#   Rscript scripts/06_survival.R

message("=== 06_survival.R: Starting Survival Analysis ===")

# ---- 1. Load packages ----
pkgs <- c("TCGAbiolinks", "dplyr", "survival", "survminer")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' not installed. Install it first.", p))
  }
}
library(TCGAbiolinks)
library(dplyr)
library(survival)
library(survminer)

# ---- 2. Load hub genes ----
hub_path <- file.path("data", "processed", "hub_genes.csv")
if (!file.exists(hub_path)) stop("Hub genes file not found.")

hub_genes <- read.csv(hub_path)
if (!"gene" %in% colnames(hub_genes)) stop("hub_genes.csv must contain 'gene' column")

genes <- hub_genes$gene
message("[info] Hub genes: ", paste(genes, collapse=", "))

# ---- 3. Download TCGA BRCA gene expression ----
message("[step] Downloading TCGA BRCA expression data ...")

query.exp <- GDCquery(
  project     = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "HTSeq - FPKM"
)

GDCdownload(query.exp)
exp_data <- GDCprepare(query.exp)

# ---- 4. Download clinical data ----
message("[step] Downloading TCGA BRCA clinical data ...")

clinical <- GDCquery_clinic(project="TCGA-BRCA", type="clinical")

# ---- 5. Prepare expression matrix ----
exp_mat <- as.data.frame(assay(exp_data))
exp_mat$gene <- rownames(exp_mat)

# ---- 6. Loop through each hub gene for survival ----

plots_dir <- "plots"
dir.create(plots_dir, showWarnings = FALSE)

for (g in genes) {
  message("[info] Running survival analysis for gene: ", g)

  if (!(g %in% exp_mat$gene)) {
    warning(sprintf("Gene %s not found in TCGA expression data.", g))
    next
  }

  gene_row <- exp_mat[exp_mat$gene == g, ]

  gene_expr <- t(gene_row[ , -ncol(gene_row)]) # remove gene column

  df <- data.frame(
    sample = colnames(exp_mat)[-ncol(exp_mat)],
    expr   = as.numeric(gene_expr)
  )

  merged <- df %>%
    inner_join(clinical, by = c("sample" = "submitter_id"))

  median_expr <- median(merged$expr, na.rm = TRUE)
  merged$group <- ifelse(merged$expr > median_expr, "High", "Low")

  s <- Surv(
    time = merged$days_to_death,
    event = merged$vital_status == "Dead"
  )

  fit <- survfit(s ~ group, data = merged)

  p <- ggsurvplot(
    fit,
    data = merged,
    pval = TRUE,
    title = paste("Survival Curve for", g),
    legend.title = "Expression",
    legend.labs = c("High", "Low")
  )

  ggsave(
    filename = file.path(plots_dir, sprintf("survival_%s.png", g)),
    plot = p,
    width = 8,
    height = 6
  )
}

message("=== 06_survival.R: Completed Successfully ===")
