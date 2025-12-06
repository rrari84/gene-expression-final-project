#!/usr/bin/env Rscript

# 05_enrichment.R
# -----------------------------------
# Performs GO and KEGG enrichment analysis using clusterProfiler
# on robust up- and down-regulated genes from 04_integration.R.
#
# Inputs (from 04_integration.R):
#   results/integration/robust_up_genes.csv
#   results/integration/robust_down_genes.csv
#
# Outputs:
#   results/enrichment/GO_up.csv
#   results/enrichment/GO_down.csv
#   results/enrichment/KEGG_up.csv
#   results/enrichment/KEGG_down.csv
#   results/enrichment/GO_up_dotplot.png
#   results/enrichment/GO_down_dotplot.png
#   results/enrichment/KEGG_up_dotplot.png
#   results/enrichment/KEGG_down_dotplot.png
#
# Usage:
#   Rscript scripts/05_enrichment.R

message("=== 05_enrichment.R: Starting Enrichment Analysis ===")

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(AnnotationDbi)
})

# Try to load array annotation for HG-U133 Plus 2.0
if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) {
  stop(
    "Package 'hgu133plus2.db' is not installed.\n",
    "Install it with:\n",
    "  BiocManager::install('hgu133plus2.db')\n",
    "or via conda:\n",
    "  conda install -c bioconda bioconductor-hgu133plus2.db"
  )
}
library(hgu133plus2.db)

int_dir    <- file.path("results", "integration")
enrich_dir <- file.path("results", "enrichment")

if (!dir.exists(int_dir)) {
  stop("Integration directory not found: ", int_dir)
}
dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# Helper: map Affy probes -> SYMBOL + ENTREZ and run GO + KEGG
# ------------------------------------------------------------------
run_enrichment_for_probes <- function(probe_ids, label) {
  message("\n[info] Running enrichment for: ", label)
  probe_ids <- unique(probe_ids)
  message("[info] # unique probes: ", length(probe_ids))

  # Map Affy PROBEID -> SYMBOL, ENTREZ using hgu133plus2.db
  message("[step] Mapping Affy probe IDs to SYMBOL and ENTREZID ...")

  anno <- AnnotationDbi::select(
    x       = hgu133plus2.db,
    keys    = probe_ids,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "PROBEID"
  ) %>%
    distinct(PROBEID, SYMBOL, ENTREZID) %>%
    filter(!is.na(ENTREZID))

  if (nrow(anno) == 0) {
    stop("No valid ENTREZ IDs after mapping for ", label)
  }

  message("[info] # probes with mapped ENTREZID: ", nrow(anno))

  entrez_ids <- unique(anno$ENTREZID)

  # ----------------- GO BP enrichment -----------------
  message("[step] Running GO BP enrichment for ", label, " ...")

  ego <- enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )

  # ----------------- KEGG enrichment ------------------
  message("[step] Running KEGG enrichment for ", label, " ...")

  ekegg <- enrichKEGG(
    gene         = entrez_ids,
    organism     = "hsa",
    pvalueCutoff = 0.05
  )

  # ----------------- Save tables ----------------------
  go_out    <- file.path(enrich_dir, paste0("GO_", label, ".csv"))
  kegg_out  <- file.path(enrich_dir, paste0("KEGG_", label, ".csv"))

  readr::write_csv(as.data.frame(ego),   go_out)
  readr::write_csv(as.data.frame(ekegg), kegg_out)

  message("[info] Saved GO results to:   ", go_out)
  message("[info] Saved KEGG results to: ", kegg_out)

  # ----------------- Save dotplots --------------------
  message("[step] Saving dotplots for ", label, " ...")

  png(file.path(enrich_dir, paste0("GO_", label, "_dotplot.png")),
      width = 1200, height = 900)
  print(dotplot(ego, showCategory = 15))
  dev.off()

  png(file.path(enrich_dir, paste0("KEGG_", label, "_dotplot.png")),
      width = 1200, height = 900)
  print(dotplot(ekegg, showCategory = 15))
  dev.off()

  invisible(list(ego = ego, ekegg = ekegg, anno = anno))
}

# ------------------------------------------------------------------
# 1. Load robust up/down gene lists (from 04_integration.R)
# ------------------------------------------------------------------
up_file   <- file.path(int_dir, "robust_up_genes.csv")
down_file <- file.path(int_dir, "robust_down_genes.csv")

if (!file.exists(up_file))   stop("File not found: ", up_file)
if (!file.exists(down_file)) stop("File not found: ", down_file)

robust_up   <- readr::read_csv(up_file, show_col_types = FALSE)
robust_down <- readr::read_csv(down_file, show_col_types = FALSE)

if (!"gene" %in% colnames(robust_up) ||
    !"gene" %in% colnames(robust_down)) {
  stop("robust_up_genes.csv and robust_down_genes.csv must contain a 'gene' column (Affy probeset IDs).")
}

# ------------------------------------------------------------------
# 2. Run enrichment for up- and down-regulated genes separately
# ------------------------------------------------------------------
res_up   <- run_enrichment_for_probes(robust_up$gene,   "up")
res_down <- run_enrichment_for_probes(robust_down$gene, "down")

message("\n=== 05_enrichment.R: Completed Successfully ===")
message("Results written under: ", enrich_dir)
