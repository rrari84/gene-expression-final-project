#!/usr/bin/env Rscript

# 05_enrichment.R
# -----------------------------------
# Performs GO and KEGG enrichment analysis using clusterProfiler.
#
# Input:
#   data/processed/common_DEGs.csv  (must contain "gene" column)
#
# Output:
#   data/results/enrichment_GO.csv
#   data/results/enrichment_KEGG.csv
#   plots/enrichment_GO_dotplot.png
#   plots/enrichment_KEGG_dotplot.png
#
# Usage:
#   Rscript scripts/05_enrichment.R

message("=== 05_enrichment.R: Starting Enrichment Analysis ===")

# ---- 1. Load packages ----
pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' not installed. Install via BiocManager::install().", p))
  }
}
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# ---- 2. Load input DEGs ----
deg_path <- file.path("data", "processed", "common_DEGs.csv")
if (!file.exists(deg_path)) stop("DEG file not found: ", deg_path)

deg_df <- read.csv(deg_path)
if (!"gene" %in% colnames(deg_df)) {
  stop("common_DEGs.csv must contain a column named 'gene'")
}

genes <- unique(deg_df$gene)
message("[info] Number of DEGs used for enrichment: ", length(genes))

# Convert gene symbols → ENTREZ IDs
message("[step] Mapping gene symbols to Entrez IDs ...")
entrez <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

if (nrow(entrez) == 0) stop("No valid Entrez IDs mapped.")

entrez_ids <- entrez$ENTREZID

# ---- 3. GO enrichment ----
message("[step] Running GO enrichment (Biological Process) ...")

ego <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# ---- 4. KEGG enrichment ----
message("[step] Running KEGG enrichment ...")

ekegg <- enrichKEGG(
  gene         = entrez_ids,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

# ---- 5. Export results ----

results_dir <- file.path("data", "results")
plots_dir   <- file.path("plots")

dir.create(results_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(plots_dir,   showWarnings=FALSE, recursive=TRUE)

write.csv(as.data.frame(ego),   file.path(results_dir, "enrichment_GO.csv"),   row.names=FALSE)
write.csv(as.data.frame(ekegg), file.path(results_dir, "enrichment_KEGG.csv"), row.names=FALSE)

# ---- 6. Save plots ----
message("[step] Saving GO + KEGG dotplots ...")

png(file.path(plots_dir, "enrichment_GO_dotplot.png"), width=1200, height=900)
print(dotplot(ego, showCategory=15))
dev.off()

png(file.path(plots_dir, "enrichment_KEGG_dotplot.png"), width=1200, height=900)
print(dotplot(ekegg, showCategory=15))
dev.off()

message("=== 05_enrichment.R: Completed Successfully ===")
