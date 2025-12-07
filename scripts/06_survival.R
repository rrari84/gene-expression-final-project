#!/usr/bin/env Rscript
# 06_survival.R
# -----------------------------------
# Performs survival analysis using TCGA BRCA clinical + expression data.
# - Expects data/processed/hub_genes.csv with column 'gene' (HGNC SYMBOL)
# - Handles Ensembl -> SYMBOL mapping and collapsing expression to SYMBOLs
# - Caches prepared TCGA data in data/processed/*.rds to avoid re-downloading
#
# Usage:
#  Rscript scripts/06_survival.R

message("=== 06_survival.R: Starting Survival Analysis ===")

# ---- 1. Load packages ----
pkgs <- c("TCGAbiolinks", "dplyr", "survival", "survminer",
          "AnnotationDbi", "org.Hs.eg.db", "readr", "SummarizedExperiment")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf("Package '%s' not installed. Install it first.", p))
  }
}
library(TCGAbiolinks)
library(dplyr)
library(survival)
library(survminer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readr)
library(SummarizedExperiment)

# ---- 2. Load hub genes ----
hub_path <- file.path("data", "processed", "hub_genes.csv")
if (!file.exists(hub_path)) stop("Hub genes file not found. Run mapping script first (05_map_probes_to_symbols.R).")

hub_genes <- read_csv(hub_path, show_col_types = FALSE)
if (!"gene" %in% colnames(hub_genes)) stop("hub_genes.csv must contain 'gene' column")
genes <- unique(as.character(hub_genes$gene))
message("[info] Hub genes: ", paste(head(genes, 30), collapse = ", "), ifelse(length(genes) > 30, " ...", ""))

# ---- 3. Load or download TCGA BRCA expression + clinical ----
expr_rds <- file.path("data", "processed", "tcga_brca_expr.rds")
clin_rds <- file.path("data", "processed", "tcga_brca_clinical.rds")

if (file.exists(expr_rds) && file.exists(clin_rds)) {
  message("[info] Loading cached TCGA data from RDS files.")
  exp_data <- readRDS(expr_rds)
  clinical <- readRDS(clin_rds)
} else {
  message("[step] Downloading TCGA BRCA expression data ... (this will take time, saving to RDS)")
  query.exp <- GDCquery(
    project     = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "HTSeq - FPKM"
  )
  GDCdownload(query.exp)
  exp_data <- GDCprepare(query.exp)

  message("[step] Downloading TCGA BRCA clinical data ...")
  clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

  dir.create(dirname(expr_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(exp_data, expr_rds)
  saveRDS(clinical, clin_rds)
  message("[info] Saved prepared TCGA data to RDS.")
}

# ---- 4. Prepare expression matrix: Ensembl -> SYMBOL ----
message("[step] Preparing expression matrix (ENSEMBL -> SYMBOL, collapse duplicates) ...")

assay_mat <- as.data.frame(assay(exp_data))  # rows = Ensembl (with versions), cols = sample barcodes
ensembl_ids <- rownames(assay_mat)
ensembl_clean <- sub("\\..*$", "", ensembl_ids)   # strip version: ENSG00000.1 -> ENSG00000

# map ENSEMBL -> SYMBOL
mapping <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = unique(ensembl_clean),
                                 keytype = "ENSEMBL",
                                 columns = c("SYMBOL"))
# mapping may have multiple lines per ENSEMBL if multiple symbols; keep first SYMBOL per ENSEMBL
mapping <- mapping %>% distinct(ENSEMBL, SYMBOL)

# attach mapping to rows
map_df <- data.frame(ensembl = ensembl_clean, orig_row = seq_along(ensembl_clean), stringsAsFactors = FALSE)
map_df <- left_join(map_df, mapping, by = c("ensembl" = "ENSEMBL"))

# add SYMBOL column to expression rows
assay_mat$ENSEMBL <- ensembl_clean
assay_mat$SYMBOL <- map_df$SYMBOL

# remove rows with NA SYMBOL
assay_mat2 <- assay_mat %>% filter(!is.na(SYMBOL))

if (nrow(assay_mat2) == 0) stop("No Ensembl IDs could be mapped to SYMBOL. Check org.Hs.eg.db mapping.")

# collapse to 1 row per SYMBOL by averaging expression across Ensembl rows mapping to same SYMBOL
expr_by_symbol <- assay_mat2 %>%
  dplyr::select(-ENSEMBL) %>%
  group_by(SYMBOL) %>%
  summarize(across(everything(), ~ mean(as.numeric(.), na.rm = TRUE))) %>%
  ungroup()

# make rownames the SYMBOL and ensure matrix-like
expr_mat_symbol <- as.data.frame(expr_by_symbol)
rownames(expr_mat_symbol) <- expr_mat_symbol$SYMBOL
expr_mat_symbol$SYMBOL <- NULL

message("[info] Expression matrix now has ", nrow(expr_mat_symbol), " genes (SYMBOLs) and ", ncol(expr_mat_symbol), " samples (barcodes).")

# ---- 5. Prepare clinical: compute os_time and os_event ----
message("[step] Preparing clinical metrics (os_time, os_event) ...")

clin <- as.data.frame(clinical, stringsAsFactors = FALSE)

# helper to pick proper days_to_last_follow_up field (common variants)
possible_last <- c("days_to_last_follow_up", "days_to_last_followup", "days_to_last_followup", "days_to_last_contact", "days_to_last_follow_up")
last_field <- intersect(possible_last, colnames(clin))
last_field <- if (length(last_field) > 0) last_field[1] else NULL

if (!"days_to_death" %in% colnames(clin) && is.null(last_field)) {
  # if expected columns missing, still proceed but warn
  warning("Clinical fields 'days_to_death' and candidate last-follow-up fields not found. os_time may be incomplete.")
}

clin$days_to_death <- if ("days_to_death" %in% colnames(clin)) as.numeric(clin$days_to_death) else NA
if (!is.null(last_field)) {
  clin$days_to_last_follow_up <- as.numeric(clin[[last_field]])
} else {
  clin$days_to_last_follow_up <- NA
}
# compute os_time and os_event
clin <- clin %>%
  mutate(
    os_time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    os_event = ifelse(tolower(as.character(vital_status)) == "dead", 1L, 0L)
  )

# some rows may have negative or NA os_time -- convert to NA and remove later
clin$os_time <- as.numeric(clin$os_time)
clin$os_event <- as.integer(clin$os_event)

message("[info] Clinical rows: ", nrow(clin), " ; with non-NA os_time: ", sum(!is.na(clin$os_time)))

# ---- 6. Align sample barcodes -> patient ID and join ----
message("[step] Aligning expression samples and clinical patient IDs ...")

sample_barcodes <- colnames(expr_mat_symbol)
# patient id is first 12 chars of barcode (TCGA-XX-YYYY-ZZ...)
patient_id <- substr(sample_barcodes, 1, 12)
expr_samples_df <- data.frame(sample_barcode = sample_barcodes, patient_id = patient_id, stringsAsFactors = FALSE)

# make sure clinical submitter_id column exists
submitter_col <- if ("submitter_id" %in% colnames(clin)) "submitter_id" else if ("bcr_patient_barcode" %in% colnames(clin)) "bcr_patient_barcode" else NULL
if (is.null(submitter_col)) {
  # try case-insensitive match
  sm <- grep("submitter|patient_barcode", colnames(clin), ignore.case = TRUE, value = TRUE)
  if (length(sm) >= 1) submitter_col <- sm[1]
}
if (is.null(submitter_col)) stop("Could not find clinical patient ID column (submitter_id or bcr_patient_barcode).")

# rename clinical patient id to patient_id for join
clin$patient_id <- as.character(clin[[submitter_col]])

message("[info] Number of expression samples: ", nrow(expr_samples_df))
message("[info] Number of clinical patients: ", length(unique(clin$patient_id)))

# ---- 7. Create output dir ----
out_plots_dir <- file.path("results", "survival")
dir.create(out_plots_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 8. Loop through hub genes and run survival analysis ----
missing_genes <- setdiff(genes, rownames(expr_mat_symbol))
if (length(missing_genes) > 0) {
  message("[warn] ", length(missing_genes), " hub genes not found in expression matrix. Examples: ",
          paste(head(missing_genes, 10), collapse = ", "), ifelse(length(missing_genes) > 10, " ...", ""))
}

for (g in genes) {
  message("[info] Running survival analysis for gene: ", g)
  if (!(g %in% rownames(expr_mat_symbol))) {
    warning(sprintf("Gene %s not found in expression matrix. Skipping.", g))
    next
  }

  expr_values <- as.numeric(expr_mat_symbol[g, ])
  df_expr <- data.frame(sample_barcode = colnames(expr_mat_symbol), expr = expr_values, stringsAsFactors = FALSE)
  df_expr$patient_id <- substr(df_expr$sample_barcode, 1, 12)

  merged <- df_expr %>%
    inner_join(clin, by = "patient_id")

  # if merge yields too few samples, skip
  if (nrow(merged) < 10) {
    warning(sprintf("Too few samples after merging for gene %s (n=%d). Skipping.", g, nrow(merged)))
    next
  }

  # compute os_time and os_event availability
  merged <- merged %>% filter(!is.na(os_time))
  if (nrow(merged) < 10) {
    warning(sprintf("Too few samples with os_time for gene %s after filtering (n=%d). Skipping.", g, nrow(merged)))
    next
  }

  median_expr <- median(merged$expr, na.rm = TRUE)
  merged$group <- ifelse(merged$expr > median_expr, "High", "Low")
  merged$group <- factor(merged$group, levels = c("High", "Low"))

  surv_obj <- Surv(time = merged$os_time, event = merged$os_event)

  # check for at least one event in each group
  tab_events <- table(merged$group, merged$os_event)
  if (ncol(tab_events) < 2 || sum(tab_events[, "1"], na.rm = TRUE) < 3) {
    warning(sprintf("Not enough events for gene %s to compute KM (events by group: %s). Skipping plot.", g,
                    paste(capture.output(print(tab_events)), collapse = " ")))
    next
  }

  fit <- tryCatch(survfit(surv_obj ~ group, data = merged), error = function(e) { warning("survfit failed: ", e); NULL })
  if (is.null(fit)) next

  p <- tryCatch(
    ggsurvplot(
      fit,
      data = merged,
      pval = TRUE,
      title = paste("Overall survival by", g),
      legend.title = "Expression",
      legend.labs = c("High", "Low"),
      risk.table = TRUE,
      risk.table.col = "strata"
    ),
    error = function(e) { warning("ggsurvplot failed: ", e); NULL }
  )
  if (is.null(p)) next

  out_file <- file.path(out_plots_dir, sprintf("survival_%s.png", g))
  # ggsurvplot returns a list; use $plot (ggplot object) for ggsave
  tryCatch({
    ggsave(filename = out_file, plot = p$plot, width = 8, height = 6)
    message("[ok] Saved survival plot to: ", out_file)
  }, error = function(e) {
    warning("Failed to save plot for gene ", g, ": ", e)
  })
}

message("=== 06_survival.R: Completed ===")
